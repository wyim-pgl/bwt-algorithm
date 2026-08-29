#!/usr/bin/env python3
"""One-to-one, boundary-aware accuracy scoring (issue #15).

The historical matcher allows one prediction to satisfy several truth records
and accepts integer-multiple/±20% periods, which is fine for regression
testing but inflates publication sensitivity. This scorer reports the strict
counterpart alongside it:

  * one-to-one assignment -- predictions and truths are matched greedily by
    descending overlap; each record participates in at most one match;
  * boundary error -- |start offset| and |end offset| per matched pair,
    reported as median/p90;
  * period accuracy -- exact, integer-multiple, and ±20% agreement rates
    reported separately instead of pooled;
  * copy-number error -- median relative error over matched pairs where both
    sides carry a copy count;
  * stratification by truth period (1-6, 7-20, 21-100, 101-2000).

Inputs are BED-like TSVs: truth needs chrom/start/end and, when available,
motif (col 4) and copies (col 5); predictions are 8-column BWTandem BED.

Usage: score_one_to_one.py TRUTH_BED PRED_BED [--min-overlap 0.5] [--json OUT]
"""
import argparse
import json
import math
import statistics
import sys
from collections import defaultdict

STRATA = ((1, 6), (7, 20), (21, 100), (101, 2000))


def load(path, need_cols=3):
    rows = []
    with open(path) as f:
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
            motif = p[3] if len(p) > 3 and p[3] else None
            copies = None
            if len(p) > 4:
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

    def dfs(i, dist):
        for j in edges[i]:
            pi = match_p.get(j)
            if pi is None or (dist.get(pi) == dist[i] + 1 and dfs(pi, dist)):
                match_t[i] = j
                match_p[j] = i
                return True
        dist[i] = INF
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


def stratum(p):
    for lo, hi in STRATA:
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
    ap.add_argument("--json", help="also write the full result as JSON")
    args = ap.parse_args()
    if not math.isfinite(args.min_overlap) or not 0 < args.min_overlap <= 1:
        ap.error("--min-overlap must be a finite value in (0, 1]")

    truth = load(args.truth_bed)
    preds = load(args.pred_bed)
    pairs, t_used, p_used = one_to_one(truth, preds, args.min_overlap)

    res = {
        "truth_records": len(truth), "pred_records": len(preds),
        "matched": len(pairs),
        "sensitivity_1to1_maxcard": pct(len(pairs), len(truth)),
        "precision_1to1_maxcard": pct(len(pairs), len(preds)),
        "min_overlap": args.min_overlap,
        "strata": {},
    }

    start_off, end_off, copy_err = [], [], []
    per_exact = per_mult = per_20 = per_scored = 0
    strat_counts = defaultdict(lambda: {"truth": 0, "matched": 0})
    for t in truth:
        strat_counts[stratum(period_of(t))]["truth"] += 1
    for i, j, ov in pairs:
        t, p = truth[i], preds[j]
        strat_counts[stratum(period_of(t))]["matched"] += 1
        start_off.append(abs(t["start"] - p["start"]))
        end_off.append(abs(t["end"] - p["end"]))
        tp, pp = period_of(t), period_of(p)
        if tp and pp:
            per_scored += 1
            if tp == pp:
                per_exact += 1
            elif pp % tp == 0 or tp % pp == 0:
                per_mult += 1
            elif abs(pp - tp) <= 0.2 * tp:
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
