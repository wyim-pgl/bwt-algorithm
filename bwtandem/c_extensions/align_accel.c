/*
 * align_accel.c — Fast C implementation of align_repeat_region loop.
 *
 * Replaces the Python-level loop that calls _align_unit_to_window per copy,
 * updates consensus via Counter, and accumulates stats. Moving the entire
 * loop to C eliminates ~3ms/call Python overhead (string slicing, object
 * creation, Counter updates) down to ~0.1ms/call.
 */
#include <stdlib.h>
#include <string.h>

/* ── Banded DP alignment (one motif copy → one window) ──────────── */

typedef struct {
    int consumed;       /* columns consumed from window */
    int mismatches;
    int ins_len;
    int del_len;
    int edit_dist;
} AlignOneResult;

/*
 * Banded semi-global DP: align motif[0..m-1] to prefix of window[0..w-1].
 * Band width = max_indel + 2.  Returns 0 on failure.
 *
 * `ops_out` receives the traceback in forward order ('M'/'S'/'D'/'I'), and
 * `i_stop_out`/`j_stop_out` where the traceback stopped — normally the origin,
 * but a traceback that walks onto an unwritten cell stops early and its ops then
 * describe the alignment from (i_stop, j_stop) onward, exactly as the reversed
 * aligned_ref/aligned_query lists do in MotifUtils._align_unit_to_window.  The
 * caller turns those into the per-copy variation records the Python loop builds.
 */
static int align_one(
    const unsigned char *motif, int m,
    const unsigned char *window, int w,
    int max_indel, int mismatch_tol,
    AlignOneResult *out,
    unsigned char *obs_bases,   /* obs_bases[motif_pos] = observed base (or 0xFF if gap) */
    int *obs_valid,             /* obs_valid[motif_pos] = 1 if base observed */
    char *ops_out,              /* capacity >= m + w */
    int *ops_len_out,
    int *i_stop_out,
    int *j_stop_out
)
{
    if (m == 0 || w == 0) return 0;

    int lower = m - max_indel;
    if (lower < 0) lower = 0;
    int upper = m + max_indel;
    if (upper > w) upper = w;
    if (lower > upper) return 0;

    int band = max_indel + 2;
    int inf = m + w + 10;

    /* Flatten DP into 1D rolling arrays to avoid malloc for large tables */
    /* We need (m+1) x (w+1) but use banded approach with two rows */
    int cols = w + 1;
    /* The full traceback table is (m+1)*cols bytes and the table is indexed by
     * i*cols + j.  Both quantities MUST be computed in 64-bit: for a large
     * pseudo-motif (e.g. an LCP candidate with period ~1e5, giving m≈1e5 and
     * w≈1.1e5) the product (m+1)*cols ≈ 1e10 overflows 32-bit int, wraps to a
     * small/negative value, malloc returns an undersized buffer, and the DP fill
     * then writes far out of bounds -> heap corruption / SIGSEGV.  Using size_t
     * for the size and every i*cols index keeps the math correct; an honestly
     * huge table simply fails to allocate and we return 0 gracefully. */
    size_t cols_sz = (size_t)cols;

    /* Guard against pathological per-copy alignments.  The traceback table is
     * O(m*w); real tandem motifs are at most a few kb, but a spurious LCP
     * candidate can carry a period of ~1e5, which would demand a ~10 GB table
     * and ~1e10-cell fill.  Refuse such cases up front (return 0 -> the caller
     * simply skips this candidate) rather than thrashing memory.  256 M cells
     * (~256 MB) corresponds to motif/window ~16 kb each — far above any genuine
     * tandem-repeat unit, so legitimate inputs are never affected. */
    if (((size_t)m + 1) * cols_sz > (size_t)256 * 1024 * 1024) {
        return 0;
    }

    /* Stack-allocate for small problems, heap for large */
    int *prev_row, *curr_row;
    char *prev_ptr, *curr_ptr;
    int on_heap = 0;

    if (cols <= 4096) {
        prev_row = (int *)__builtin_alloca(cols * sizeof(int));
        curr_row = (int *)__builtin_alloca(cols * sizeof(int));
        prev_ptr = (char *)__builtin_alloca(cols * sizeof(char));
        curr_ptr = (char *)__builtin_alloca(cols * sizeof(char));
    } else {
        prev_row = (int *)malloc(cols_sz * sizeof(int));
        curr_row = (int *)malloc(cols_sz * sizeof(int));
        prev_ptr = (char *)malloc(cols_sz * sizeof(char));
        curr_ptr = (char *)malloc(cols_sz * sizeof(char));
        on_heap = 1;
        if (!prev_row || !curr_row || !prev_ptr || !curr_ptr) {
            free(prev_row); free(curr_row); free(prev_ptr); free(curr_ptr);
            return 0;
        }
    }

    /* Full ptr table needed for traceback — allocate on heap (64-bit size). */
    /* calloc, not malloc: the fill loop below writes only [j_min..j_max] of each
     * row, so nothing else in this table is ever assigned. 0 is the traceback's
     * Stop code, matching _accelerators.pyx, which zero-fills its whole table. */
    char *ptr_table = (char *)calloc(((size_t)m + 1) * cols_sz, sizeof(char));
    if (!ptr_table) {
        if (on_heap) { free(prev_row); free(curr_row); free(prev_ptr); free(curr_ptr); }
        return 0;
    }

    /* Init row 0 */
    prev_row[0] = 0;
    ptr_table[0] = 0;
    for (int j = 1; j <= w; j++) {
        prev_row[j] = j;
        ptr_table[j] = 'I';
    }

    /* Fill DP */
    for (int i = 1; i <= m; i++) {
        int j_min = i - band;
        if (j_min < 1) j_min = 1;
        int j_max = i + band;
        if (j_max > w) j_max = w;

        /* Init column 0 for this row. The backpointer must go into ptr_table,
         * which is what the traceback reads — curr_ptr is never read. Writing only
         * curr_ptr left column 0 uninitialised, so any alignment that walked off
         * the left edge read heap garbage and the whole caller became
         * nondeterministic. _accelerators.pyx sets ptr[i*cols] = D here. */
        curr_row[0] = i;
        curr_ptr[0] = 'D';
        ptr_table[(size_t)i * cols_sz] = 'D';
        /* Set out-of-band to inf */
        for (int j = 1; j < j_min; j++) curr_row[j] = inf;
        for (int j = j_max + 1; j <= w; j++) curr_row[j] = inf;

        for (int j = j_min; j <= j_max; j++) {
            int sub_cost = prev_row[j - 1] + (motif[i - 1] != window[j - 1]);
            int del_cost = prev_row[j] + 1;
            int ins_cost = curr_row[j - 1] + 1;

            int best = sub_cost;
            char bp = (motif[i - 1] == window[j - 1]) ? 'M' : 'S';

            if (del_cost < best) { best = del_cost; bp = 'D'; }
            if (ins_cost < best) { best = ins_cost; bp = 'I'; }

            curr_row[j] = best;
            ptr_table[(size_t)i * cols_sz + j] = bp;
        }

        /* Swap rows */
        int *tmp_r = prev_row; prev_row = curr_row; curr_row = tmp_r;
        char *tmp_p = prev_ptr; prev_ptr = curr_ptr; curr_ptr = tmp_p;
    }

    /* Find best endpoint in last row (prev_row holds row m) */
    int best_j = -1, best_cost = inf;
    for (int j = lower; j <= upper; j++) {
        if (prev_row[j] < best_cost) {
            best_cost = prev_row[j];
            best_j = j;
        }
    }

    if (best_j <= 0 || best_cost >= inf) {
        free(ptr_table);
        if (on_heap) { free(prev_row); free(curr_row); free(prev_ptr); free(curr_ptr); }
        return 0;
    }

    /* Traceback to compute stats */
    int mm = 0, ins = 0, del_ = 0;
    memset(obs_valid, 0, m * sizeof(int));

    int i = m, j = best_j;
    int n_ops = 0;
    while (i > 0 || j > 0) {
        char op = ptr_table[(size_t)i * cols_sz + j];
        if (op == 'M' || op == 'S') {
            obs_bases[i - 1] = window[j - 1];
            obs_valid[i - 1] = 1;
            if (op == 'S') mm++;
            ops_out[n_ops++] = op;
            i--; j--;
        } else if (op == 'D') {
            del_++;
            ops_out[n_ops++] = op;
            i--;
        } else if (op == 'I') {
            ins++;
            ops_out[n_ops++] = op;
            j--;
        } else {
            break;
        }
    }

    /* The traceback walks backwards; the extraction below reads forwards. */
    for (int a = 0, b = n_ops - 1; a < b; a++, b--) {
        char t = ops_out[a]; ops_out[a] = ops_out[b]; ops_out[b] = t;
    }
    *ops_len_out = n_ops;
    *i_stop_out = i;
    *j_stop_out = j;

    free(ptr_table);
    if (on_heap) { free(prev_row); free(curr_row); free(prev_ptr); free(curr_ptr); }

    if (mm > mismatch_tol || ins > max_indel || del_ > max_indel)
        return 0;

    out->consumed = best_j;
    out->mismatches = mm;
    out->ins_len = ins;
    out->del_len = del_;
    out->edit_dist = best_cost;
    return 1;
}

/* ── Per-copy variation records ─────────────────────────────────
 *
 * The Python loop keeps `variations`, `copy_sequences` and `error_counts`
 * per copy; this loop used to return them empty, so the STRfinder Variations
 * column was populated or blank depending on whether libalign_accel happened to
 * load.  Rather than format strings here, the C reports the structure and the
 * Python side renders it, so the output format still lives in one place.
 */

#define VAR_SUB 0
#define VAR_INS 1
#define VAR_DEL 2

typedef struct {
    int *meta;              /* 5 ints per record: copy, kind, pos, len, char_off */
    int cap;
    int n;
    unsigned char *chars;   /* sub: ref,alt.  ins: the inserted bases.  del: none */
    int chars_cap;
    int chars_n;
    int overflow;
} VarSink;

static void var_emit(VarSink *vs, int copy_idx, int kind, int pos, int len,
                     const unsigned char *chars, int n_chars)
{
    if (vs->n >= vs->cap || vs->chars_n + n_chars > vs->chars_cap) {
        vs->overflow = 1;
        return;
    }
    int *m = vs->meta + (size_t)vs->n * 5;
    m[0] = copy_idx;
    m[1] = kind;
    m[2] = pos;
    m[3] = len;
    m[4] = n_chars ? vs->chars_n : -1;
    for (int k = 0; k < n_chars; k++) vs->chars[vs->chars_n + k] = chars[k];
    vs->chars_n += n_chars;
    vs->n++;
}

/*
 * Transcription of the operation-extraction loop in
 * MotifUtils._align_unit_to_window: walk the aligned columns, coalescing runs of
 * inserted and deleted bases, and record a sub/ins/del per group.  Positions and
 * emission order match that loop exactly — an insertion is reported at the ref
 * position *before* the column advances, a deletion at the position after, and a
 * pending insertion is flushed ahead of a pending deletion at the end.
 */
static void extract_variations(
    const unsigned char *motif, const unsigned char *window,
    const char *ops, int ops_len, int i_stop, int j_stop,
    int copy_idx, VarSink *vs)
{
    int ref_pos = 0;
    int mi = i_stop, wi = j_stop;
    int pend_ins_len = 0, pend_ins_pos = 0, pend_ins_wi = 0;
    int pend_del_len = 0, pend_del_pos = 0;

    for (int t = 0; t < ops_len; t++) {
        char op = ops[t];

        if (op == 'I') {                       /* ref column is a gap */
            if (pend_ins_len == 0) { pend_ins_pos = ref_pos; pend_ins_wi = wi; }
            pend_ins_len++;
            wi++;
            continue;
        }

        if (pend_ins_len) {
            var_emit(vs, copy_idx, VAR_INS, pend_ins_pos, pend_ins_len,
                     window + pend_ins_wi, pend_ins_len);
            pend_ins_len = 0;
        }

        ref_pos++;

        if (op == 'D') {                       /* query column is a gap */
            if (pend_del_len == 0) pend_del_pos = ref_pos;
            pend_del_len++;
            mi++;
            continue;
        }

        if (pend_del_len) {
            var_emit(vs, copy_idx, VAR_DEL, pend_del_pos, pend_del_len, NULL, 0);
            pend_del_len = 0;
        }

        if (op == 'S') {
            unsigned char pair[2];
            pair[0] = motif[mi];
            pair[1] = window[wi];
            var_emit(vs, copy_idx, VAR_SUB, ref_pos, 1, pair, 2);
        }
        mi++; wi++;
    }

    if (pend_ins_len)
        var_emit(vs, copy_idx, VAR_INS, pend_ins_pos, pend_ins_len,
                 window + pend_ins_wi, pend_ins_len);
    if (pend_del_len)
        var_emit(vs, copy_idx, VAR_DEL, pend_del_pos, pend_del_len, NULL, 0);
}

/* ── Full align_repeat_region loop in C ────────────────────────── */

typedef struct {
    int copies;
    int consumed_length;
    int total_mismatches;
    int total_insertions;
    int total_deletions;
    int max_errors_per_copy;
    double mismatch_rate;
    /* consensus is written into caller-provided buffer */
} AlignRegionResult;

/*
 * align_repeat_region_c: replaces the Python loop.
 *
 * text:           full sequence bytes
 * text_len:       length of text
 * start, end:     region to align within
 * motif:          initial motif template
 * motif_len:      length of motif
 * mismatch_frac:  per-position mismatch tolerance fraction
 * max_indel:      max indels per copy alignment
 * min_copies:     minimum copies required
 * out:            result struct
 * consensus_out:  buffer of at least motif_len bytes for final consensus
 * copy_consumed_out / copy_errors_out:
 *                 per-copy bases consumed and error count, capacity copies_cap
 * var_meta_out / var_chars_out:
 *                 per-copy variation records, capacity var_cap / var_chars_cap
 * n_vars_out:     number of variation records written
 *
 * Returns 1 on success, 0 if min_copies not reached, 2 if an output buffer was
 * too small (the caller then re-runs the Python loop, which has no such bound).
 */
int align_repeat_region_c(
    const unsigned char *text, int text_len,
    int start, int end,
    const unsigned char *motif, int motif_len,
    double mismatch_frac, int max_indel, int min_copies,
    AlignRegionResult *out,
    unsigned char *consensus_out,
    int *copy_consumed_out, int *copy_errors_out, int copies_cap,
    int *var_meta_out, unsigned char *var_chars_out,
    int var_cap, int var_chars_cap,
    int *n_vars_out
)
{
    *n_vars_out = 0;
    if (motif_len <= 0 || text_len <= 0) return 0;
    if (start < 0) start = 0;
    if (end <= start) end = text_len;
    if (end > text_len) end = text_len;

    int tolerance = (int)(motif_len * mismatch_frac);
    if (tolerance < 1) tolerance = 1;

    /* Per-position base counts for consensus: counts[pos][base], base: A=0,C=1,G=2,T=3 */
    int *counts = (int *)calloc(motif_len * 4, sizeof(int));
    unsigned char *current_motif = (unsigned char *)malloc(motif_len);
    unsigned char *obs_bases = (unsigned char *)malloc(motif_len);
    int *obs_valid = (int *)malloc(motif_len * sizeof(int));
    /* traceback is at most m + w steps, and w <= motif_len + max_indel */
    char *ops = (char *)malloc((size_t)2 * motif_len + max_indel + 2);
    if (!counts || !current_motif || !obs_bases || !obs_valid || !ops) {
        free(counts); free(current_motif); free(obs_bases); free(obs_valid); free(ops);
        return 0;
    }

    VarSink vs;
    vs.meta = var_meta_out; vs.cap = var_cap; vs.n = 0;
    vs.chars = var_chars_out; vs.chars_cap = var_chars_cap; vs.chars_n = 0;
    vs.overflow = 0;

    memcpy(current_motif, motif, motif_len);

    int copies = 0;
    int total_mm = 0, total_ins = 0, total_del = 0;
    int max_err = 0;

    /* Same extension bound as MotifUtils.align_repeat_region:
     *   min(seq_len, max(end, start + motif_len*min_copies) + max(motif_len*3, max_indel*4))
     * The old form, end + motif_len*3 + max_indel*4, let the C loop run past where
     * the Python one stopped, so the two implementations placed array boundaries
     * differently. */
    int pos = start;
    int base = end;
    if (start + motif_len * min_copies > base) base = start + motif_len * min_copies;
    int slack = motif_len * 3;
    if (max_indel * 4 > slack) slack = max_indel * 4;
    int safety = base + slack;
    if (safety > text_len) safety = text_len;

    while (pos < safety) {
        int window_end = pos + motif_len + max_indel;
        if (window_end > text_len) window_end = text_len;
        int w_len = window_end - pos;
        if (w_len < motif_len - max_indel) break;

        AlignOneResult ar;
        int ops_len = 0, i_stop = 0, j_stop = 0;
        int ok = align_one(current_motif, motif_len,
                           text + pos, w_len,
                           max_indel, tolerance,
                           &ar, obs_bases, obs_valid,
                           ops, &ops_len, &i_stop, &j_stop);
        if (!ok || ar.consumed == 0) break;

        if (copies >= copies_cap) { vs.overflow = 1; break; }
        copy_consumed_out[copies] = ar.consumed;
        copy_errors_out[copies] = ar.mismatches + ar.ins_len + ar.del_len;
        extract_variations(current_motif, text + pos, ops, ops_len,
                           i_stop, j_stop, copies + 1, &vs);

        copies++;
        total_mm += ar.mismatches;
        total_ins += ar.ins_len;
        total_del += ar.del_len;
        int err = ar.mismatches + ar.ins_len + ar.del_len;
        if (err > max_err) max_err = err;

        /* Update per-position counts */
        for (int k = 0; k < motif_len; k++) {
            if (obs_valid[k]) {
                unsigned char b = obs_bases[k];
                int bi = -1;
                switch (b) {
                    case 'A': case 'a': bi = 0; break;
                    case 'C': case 'c': bi = 1; break;
                    case 'G': case 'g': bi = 2; break;
                    case 'T': case 't': bi = 3; break;
                }
                if (bi >= 0) counts[k * 4 + bi]++;
            }
        }

        /* Update consensus motif */
        for (int k = 0; k < motif_len; k++) {
            int best_cnt = 0, best_base = current_motif[k];
            static const unsigned char bases[4] = {'A', 'C', 'G', 'T'};
            for (int b = 0; b < 4; b++) {
                if (counts[k * 4 + b] > best_cnt) {
                    best_cnt = counts[k * 4 + b];
                    best_base = bases[b];
                }
            }
            current_motif[k] = best_base;
        }

        pos += ar.consumed;
    }

    /* A truncated result would disagree with the Python loop, which has no
     * capacity bound; say so and let the caller fall back to it instead. */
    if (vs.overflow) {
        free(counts); free(current_motif); free(obs_bases); free(obs_valid); free(ops);
        return 2;
    }

    if (copies < min_copies) {
        free(counts); free(current_motif); free(obs_bases); free(obs_valid); free(ops);
        return 0;
    }

    int consumed_len = pos - start;
    if (consumed_len <= 0) {
        free(counts); free(current_motif); free(obs_bases); free(obs_valid); free(ops);
        return 0;
    }

    out->copies = copies;
    out->consumed_length = consumed_len;
    out->total_mismatches = total_mm;
    out->total_insertions = total_ins;
    out->total_deletions = total_del;
    out->max_errors_per_copy = max_err;

    int denom = copies * motif_len;
    out->mismatch_rate = denom > 0 ? (double)total_mm / denom : 0.0;

    memcpy(consensus_out, current_motif, motif_len);
    *n_vars_out = vs.n;

    free(counts); free(current_motif); free(obs_bases); free(obs_valid); free(ops);
    return 1;
}
