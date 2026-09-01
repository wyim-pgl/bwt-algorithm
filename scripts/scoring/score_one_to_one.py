#!/usr/bin/env python3
"""One-to-one, boundary-aware accuracy scoring (issue #15).

The historical matcher allows one prediction to satisfy several truth records
and accepts integer-multiple/±20% periods, which is fine for regression
testing but inflates publication sensitivity. This scorer reports the strict
counterpart alongside it:

  * one-to-one assignment -- maximum-cardinality bipartite matching
    (Hopcroft-Karp); each record participates in at most one match;
  * boundary error -- |start offset| and |end offset| per matched pair,
    reported as median/p90;
  * period-length agreement -- exact, integer-multiple, and ±20% rates
    reported separately instead of pooled (length agreement, NOT motif
    sequence identity);
  * copy-number error -- median relative error over matched pairs where both
    sides carry a copy count;
  * stratification by TRUTH period band (1-6, 7-20, 21-100, 101-2000).

Inputs are BED-like TSVs (gzip accepted): truth needs chrom/start/end and,
when available, motif (col 4) and copies (col 5). Column-5 semantics differ
between tools -- BWTandem BEDs carry a copy count there, but the converted
ULTRA/tantan/TRF baselines carry the PERIOD (see convert_to_bed.py) -- so
--pred-col5 must say which; 'period' disables the copy-error metric rather
than silently comparing a period against a copy count. --pred-motif-is-
sequence marks BEDs whose column 4 is the full array sequence (TRF), which
disables prediction-period metrics while truth-based strata still fill.
--truth-chroms-only drops predictions on sequences absent from the truth
(e.g. unplaced scaffolds a chr1-22XY catalog can never match), so precision
denominators are comparable across tools run on different FASTA scopes.

Usage: score_one_to_one.py TRUTH_BED PRED_BED [--min-overlap 0.5]
       [--pred-col5 {copies,period}] [--pred-motif-is-sequence]
       [--truth-chroms-only] [--json OUT]
"""
import argparse
import gzip
import json
import math
import statistics
import sys
from collections import defaultdict

STRATA = ((1, 6), (7, 20), (21, 100), (101, 2000))


def parse_strata(spec):
    """'1-6,7-20,101-500' -> ((1,6),(7,20),(101,500)); validated."""
    bands = []
    for part in spec.split(","):
        try:
            lo, hi = part.strip().split("-")
            lo, hi = int(lo), int(hi)
        except ValueError:
            sys.exit(f"FATAL: bad stratum '{part}' (expected LO-HI)")
        if lo <= 0 or hi < lo:
            sys.exit(f"FATAL: bad stratum '{part}' (need 0 < LO <= HI)")
        if bands and lo <= bands[-1][1]:
            sys.exit(f"FATAL: strata must be ascending and non-overlapping ('{part}')")
        bands.append((lo, hi))
    if not bands:
        sys.exit("FATAL: --strata is empty")
    return tuple(bands)


def load(path, need_cols=3, col5="copies", motif_is_sequence=False):
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for i, line in enumerate(f, 1):
            p = line.rstrip("\n").split("\t")
            if len(p) < need_cols or line.startswith(("#", "track")):
                continue
            try:
                chrom, s, e = p[0], int(p[1]), int(p[2])
            except ValueError:
                sys.exit(f"FATAL {path}:{i}: bad coordinates")
            if e <= s:
                sys.exit(f"FATAL {path}:{i}: empty/inverted interval")
            motif = p[3] if len(p) > 3 and p[3] and not motif_is_sequence else None
            copies = None
            if len(p) > 4 and col5 == "copies":
                try:
                    copies = float(p[4])
                except ValueError:
                    pass
            rows.append({"chrom": chrom, "start": s, "end": e,
                         "motif": motif, "copies": copies, "line": i})
    if not rows:
        sys.exit(f"FATAL: {path} contains no usable records")
    return rows


def overlap(a, b):
    return max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]))


def _eligible_edges(truth, preds, min_frac):
    """Sweep-line candidate generation: only overlapping pairs are materialized.

    Both sides are sorted by start per chromosome and predictions are held in
    an active window, so memory is O(active) and time O((T+P) log + E) instead
    of the Cartesian O(T*P) within each chromosome.
    """
    by_chrom_t = defaultdict(list)
    by_chrom_p = defaultdict(list)
    for i, t in enumerate(truth):
        by_chrom_t[t["chrom"]].append(i)
    for j, p in enumerate(preds):
        by_chrom_p[p["chrom"]].append(j)
    edges = defaultdict(dict)   # i -> {j: overlap}
    for chrom, tis in by_chrom_t.items():
        pjs = sorted(by_chrom_p.get(chrom, ()), key=lambda j: preds[j]["start"])
        tis = sorted(tis, key=lambda i: truth[i]["start"])
        import heapq
        active = []            # (end, j) heap
        k = 0
        for i in tis:
            t = truth[i]
            while k < len(pjs) and preds[pjs[k]]["start"] < t["end"]:
                pj = pjs[k]
                heapq.heappush(active, (preds[pj]["end"], pj))
                k += 1
            while active and active[0][0] <= t["start"]:
                heapq.heappop(active)
            for _, j in active:
                pcall = preds[j]
                ov = overlap(t, pcall)
                if ov <= 0:
                    continue
                # reciprocal: the overlap must cover min_frac of BOTH records
                if ov < min_frac * (t["end"] - t["start"]):
                    continue
                if ov < min_frac * (pcall["end"] - pcall["start"]):
                    continue
                edges[i][j] = ov
    return edges


def one_to_one(truth, preds, min_frac):
    """Maximum-cardinality one-to-one assignment (Hopcroft-Karp).

    A greedy descending-overlap pass is biased: one long candidate can consume
    the only partner of a truth record when a slightly smaller edge would have
    permitted two matches, so sensitivity/precision depended on a heuristic.
    Maximum cardinality removes that bias. Overlap is NOT a secondary
    objective: among equal-cardinality matchings the pairing (and therefore
    the boundary/period/copy statistics) is arbitrary, which is disclosed in
    the output header.
    """
    edges = _eligible_edges(truth, preds, min_frac)
    INF = float("inf")
    match_t = {}   # i -> j
    match_p = {}   # j -> i

    def bfs():
        dist = {}
        q = []
        for i in edges:
            if i not in match_t:
                dist[i] = 0
                q.append(i)
        found = False
        head = 0
        while head < len(q):
            i = q[head]; head += 1
            for j in edges[i]:
                pi = match_p.get(j)
                if pi is None:
                    found = True
                elif pi not in dist:
                    dist[pi] = dist[i] + 1
                    q.append(pi)
        return dist, found

    def dfs(root, dist):
        # Iterative: an augmenting path in a dense contested pileup can span
        # thousands of layers, and the recursive form dies at Python's
        # ~1000-frame default (reproduced on a 3000-interval half-shifted
        # chain). Each stack frame is (vertex, iterator over its edges).
        stack = [(root, iter(edges[root]))]
        path = []  # (i, j) edges taken down the current branch
        while stack:
            i, it = stack[-1]
            advanced = False
            for j in it:
                pi = match_p.get(j)
                if pi is None:
                    # augmenting path found: flip every edge on the path
                    path.append((i, j))
                    for pi_, pj_ in path:
                        match_t[pi_] = pj_
                        match_p[pj_] = pi_
                    return True
                if dist.get(pi) == dist[i] + 1:
                    path.append((i, j))
                    stack.append((pi, iter(edges[pi])))
                    advanced = True
                    break
            if not advanced:
                dist[i] = INF
                stack.pop()
                if path:
                    path.pop()
        return False

    while True:
        dist, found = bfs()
        if not found:
            break
        for i in list(edges):
            if i not in match_t:
                dfs(i, dist)

    pairs = [(i, j, edges[i][j]) for i, j in match_t.items()]
    return pairs, set(match_t), set(match_p)


def period_of(rec):
    return len(rec["motif"]) if rec["motif"] else None


def stratum(p, strata=STRATA):
    for lo, hi in strata:
        if p is not None and lo <= p <= hi:
            return f"{lo}-{hi}"
    return "other"


def pct(n, d):
    return 100.0 * n / d if d else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("truth_bed")
    ap.add_argument("pred_bed")
    ap.add_argument("--min-overlap", type=float, default=0.5,
                    help="reciprocal overlap fraction required (default 0.5)")
    ap.add_argument("--pred-col5", choices=("copies", "period"), default="copies",
                    help="meaning of the prediction BED's 5th column; 'period' "
                         "disables the copy-error metric (default: copies)")
    ap.add_argument("--truth-col5", choices=("copies", "period"), default="copies",
                    help="meaning of the truth BED's 5th column (symmetric with "
                         "--pred-col5; 'period' disables the copy-error metric)")
    ap.add_argument("--pred-motif-is-sequence", action="store_true",
                    help="prediction column 4 is the full array sequence, not a "
                         "motif (TRF-style); disables prediction-period metrics")
    ap.add_argument("--strata", default="1-6,7-20,21-100,101-2000",
                    help="comma list of ascending truth period bands LO-HI "
                         "(default: 1-6,7-20,21-100,101-2000)")
    ap.add_argument("--truth-chroms-only", action="store_true",
                    help="drop predictions on sequences absent from the truth "
                         "before computing precision")
    ap.add_argument("--json", help="also write the full result as JSON")
    args = ap.parse_args()
    if not math.isfinite(args.min_overlap) or not 0 < args.min_overlap <= 1:
        ap.error("--min-overlap must be a finite value in (0, 1]")
    strata_bands = parse_strata(args.strata)

    truth = load(args.truth_bed, col5=args.truth_col5)
    preds = load(args.pred_bed, col5=args.pred_col5,
                 motif_is_sequence=args.pred_motif_is_sequence)
    pred_records_loaded = len(preds)
    dropped_off_truth_chroms = 0
    if args.truth_chroms_only:
        truth_chroms = {t["chrom"] for t in truth}
        preds = [p for p in preds if p["chrom"] in truth_chroms]
        dropped_off_truth_chroms = pred_records_loaded - len(preds)
        if not preds:
            sys.exit("FATAL: no predictions left on the truth's sequences")
    pairs, t_used, p_used = one_to_one(truth, preds, args.min_overlap)

    res = {
        "truth_records": len(truth), "pred_records": len(preds),
        "pred_records_loaded": pred_records_loaded,
        "dropped_off_truth_chroms": dropped_off_truth_chroms,
        "pred_col5": args.pred_col5,
        "pred_motif_is_sequence": args.pred_motif_is_sequence,
        "matched": len(pairs),
        "sensitivity_1to1_maxcard": pct(len(pairs), len(truth)),
        "precision_1to1_maxcard": pct(len(pairs), len(preds)),
        "min_overlap": args.min_overlap,
        "strata_spec": args.strata,
        "strata": {},
    }

    start_off, end_off, copy_err = [], [], []
    per_exact = per_mult = per_20 = per_scored = 0
    strat_counts = defaultdict(lambda: {"truth": 0, "matched": 0})
    for t in truth:
        strat_counts[stratum(period_of(t), strata_bands)]["truth"] += 1
    for i, j, ov in pairs:
        t, p = truth[i], preds[j]
        strat_counts[stratum(period_of(t), strata_bands)]["matched"] += 1
        start_off.append(abs(t["start"] - p["start"]))
        end_off.append(abs(t["end"] - p["end"]))
        tp, pp = period_of(t), period_of(p)
        if tp and pp:
            per_scored += 1
            if tp == pp:
                per_exact += 1
            elif pp % tp == 0 or tp % pp == 0:
                per_mult += 1
            elif (max(tp, pp) - min(tp, pp)) / min(tp, pp) <= 0.2:
                # smaller-period-relative, matching the permissive matcher's
                # periods_compatible() so the strict and permissive rates
                # decompose the same predicate
                per_20 += 1
        if (t["copies"] is not None and t["copies"] > 0
                and p["copies"] is not None):
            copy_err.append(abs(p["copies"] - t["copies"]) / t["copies"])

    def q(v, f):
        return round(f(v), 2) if v else None
    res["boundary"] = {
        "start_offset_median": q(start_off, statistics.median),
        "start_offset_p90": q(start_off, lambda v: statistics.quantiles(v, n=10)[8]) if len(start_off) >= 10 else None,
        "end_offset_median": q(end_off, statistics.median),
        "end_offset_p90": q(end_off, lambda v: statistics.quantiles(v, n=10)[8]) if len(end_off) >= 10 else None,
        "exact_start_pct": round(pct(sum(1 for v in start_off if v == 0), len(pairs)), 2),
        "exact_end_pct": round(pct(sum(1 for v in end_off if v == 0), len(pairs)), 2),
    }
    res["period"] = {
        "scored_pairs": per_scored,
        "exact_pct": round(pct(per_exact, per_scored), 2),
        "integer_multiple_pct": round(pct(per_mult, per_scored), 2),
        "within_20pct_pct": round(pct(per_20, per_scored), 2),
        "outside_pct": round(pct(per_scored - per_exact - per_mult - per_20, per_scored), 2),
    }
    res["copies"] = {
        "scored_pairs": len(copy_err),
        "rel_error_median_pct": q([100 * e for e in copy_err], statistics.median),
    }
    for k in sorted(strat_counts):
        c = strat_counts[k]
        res["strata"][k] = {**c, "sensitivity_pct": round(pct(c["matched"], c["truth"]), 2)}

    print(f"max-cardinality 1:1 (Hopcroft-Karp) @ reciprocal {args.min_overlap:.0%}")
    print("  note: overlap is not a secondary objective; among equal-cardinality "
          "matchings the pairing is arbitrary, so boundary/period/copy stats "
          "carry that indeterminacy")
    if dropped_off_truth_chroms:
        print(f"  note: {dropped_off_truth_chroms:,} of {pred_records_loaded:,} "
              f"predictions dropped (sequences absent from the truth)")
    print(f"  truth {res['truth_records']:,}  preds {res['pred_records']:,}  "
          f"matched {res['matched']:,}")
    print(f"  sensitivity {res['sensitivity_1to1_maxcard']:.2f}%   "
          f"precision {res['precision_1to1_maxcard']:.2f}%")
    b = res["boundary"]
    print(f"  boundary |Δstart| median {b['start_offset_median']} bp "
          f"(exact {b['exact_start_pct']}%), |Δend| median {b['end_offset_median']} bp "
          f"(exact {b['exact_end_pct']}%)")
    pr = res["period"]
    print(f"  period: exact {pr['exact_pct']}%  int-multiple {pr['integer_multiple_pct']}%  "
          f"±20% {pr['within_20pct_pct']}%  outside {pr['outside_pct']}%")
    print(f"  copies: median rel err {res['copies']['rel_error_median_pct']}% "
          f"over {res['copies']['scored_pairs']:,} pairs")
    for k, v in res["strata"].items():
        print(f"  stratum {k:>9}: {v['matched']:,}/{v['truth']:,} = {v['sensitivity_pct']}%")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(res, f, indent=1)
        print(f"json -> {args.json}")


if __name__ == "__main__":
    main()
