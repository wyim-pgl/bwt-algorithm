#!/usr/bin/env python3
"""Input tables for two figures the manuscript does not yet have.

FIGURE T — stacked tool tracks at representative loci.
  One row per call, long format, so a plotting layer can lane the tools vertically
  and draw each call as a segment. Shows fragmentation and boundary disagreement
  directly, which no summary statistic in the paper does.
  Output: figure_inputs/tracks_<locus>.tsv  +  tracks_index.tsv

FIGURE U — UpSet over tool agreement.
  `bedtools multiinter` segments the genome wherever the set of overlapping tools
  changes, and reports that set per segment. Counting segments and bp per distinct
  set is exactly UpSet's input.
  Output: figure_inputs/upset_<genome>_combinations.tsv  +  upset_<genome>_segments.bed.gz

Both are derived only from BEDs already on disk; no detection compute.

Usage:  python3 make_figure_inputs.py [--genome colcen|maize|human|all]
"""
import argparse
import gzip
import os
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
EXP1 = os.path.dirname(HERE)
FILIP = "/data/gpfs/assoc/pgl/filip/bwtandem_results"
BEDS = f"{FILIP}/beds"
GT = f"{FILIP}/ground_truth"
OUT = os.path.join(HERE, "figure_inputs")
BEDTOOLS = os.environ.get(  # same default the scorers use
    "BEDTOOLS", "/data/gpfs/assoc/pgl/bin/bedtools2/bin/bedtools")

# Which BED represents each tool. Where a tool was re-run because its published
# invocation could not reach the periods being scored (WP0-B), the re-run is used
# and the published file is not, so the figure and the tables agree.
GENOMES = {
    "colcen": {
        "reference": (f"{GT}/colcen_cen180.bed", "CEN180 monomer"),
        "tools": [
            ("BWTandem", f"{HERE}/out/remeas_colcen.bed"),
            ("TRF",      f"{BEDS}/trf/Col-CEN_v1.2_output.bed"),
            ("ULTRA",    f"{HERE}/beds/ultra_colcen_p500.bed"),      # -p 500 re-run
            ("mreps",    f"{HERE}/beds/mreps_colcen_150_400.bed"),   # 150-400 re-run
            ("tantan",   f"{BEDS}/tantan/Col-CEN_v1.2_output.bed"),
            ("TRASH",    f"{BEDS}/trash/Col-CEN_v1.2_trash.bed"),
            ("NCRF",     f"{BEDS}/ncrf/Col-CEN_v1.2_output.bed"),
        ],
    },
    "maize": {
        "reference": (f"{GT}/mo17_knob180_arrays.bed", "knob180 array"),
        "tools": [
            ("BWTandem", f"{HERE}/out/remeas_maize.bed"),
            ("TRF",      f"{BEDS}/trf/GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic_trf_exp3B_satellite.bed"),
            ("ULTRA",    f"{HERE}/beds/ultra_maize_exp3B.bed"),
            ("tantan",   f"{BEDS}/tantan/GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic_tantan_exp3B_satellite.bed"),
            ("TRASH",    f"{BEDS}/trash/GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic_trash_denovo_exp3B_satellite.bed"),
        ],
    },
    "human": {
        "reference": (f"{EXP1}/data/adotto_primary.bed", "adotto region"),
        "tools": [
            ("BWTandem", f"{HERE}/out/remeas_human.bed"),
            ("TRF",      f"{BEDS}/trf/GCA_000001405.15_GRCh38_genomic_output.bed"),
            ("ULTRA",    f"{HERE}/beds/ultra_human_GCA.bed"),
            ("tantan",   f"{BEDS}/tantan/GCA_000001405.15_GRCh38_genomic_output.bed"),
            ("TRASH",    f"{BEDS}/trash/GCA_000001405.15_GRCh38_genomic_trash.bed"),
            # mreps is excluded: every one of its 18,071,983 intervals carries
            # chrom chr4 (WP0-B finding 5), so including it would put a fictitious
            # tool on one chromosome and none on the other 23.
        ],
    },
}

# Windows chosen to be wide enough to show array structure and narrow enough that
# individual calls remain distinguishable. Full array extents are in tracks_index.
LOCI = [
    ("arabidopsis_cen1", "colcen", "Chr1", 15_600_000, 15_680_000,
     "Chr1 centromere interior, inside the 14.38-17.87 Mb CEN180 array"),
    ("maize_knob180",    "maize",  "CM039155.1", 4_700_000, 4_780_000,
     "knob180 array on CM039155.1, inside the 4.64-6.19 Mb curated extent"),
    ("human_gap",        "human",  None, None, None,
     "densest 80 kb of BWTandem calls with no adotto region (chosen at runtime)"),
]

# Added after the first pass: the 80 kb windows above sit inside the arrays, so every
# tool covers them and only fragmentation is visible. These two shapes add what that
# view cannot show — individual call boundaries, and what each tool does where the
# curated array actually stops.
LOCI += [
    ("arabidopsis_cen1_zoom", "colcen", "Chr1", 15_640_000, 15_655_000,
     "15 kb inside the Chr1 CEN180 array: individual monomer-level calls resolvable"),
    ("arabidopsis_cen1_edge", "colcen", "Chr1", 14_362_000, 14_392_000,
     "30 kb spanning the Chr1 centromere start at 14,377,192: boundary disagreement"),
    ("maize_knob180_zoom",    "maize",  "CM039155.1", 4_720_000, 4_735_000,
     "15 kb inside the knob180 array"),
    ("maize_knob180_edge",    "maize",  "CM039155.1", 4_629_000, 4_659_000,
     "30 kb spanning the knob180 array start at 4,643,853: boundary disagreement"),
    ("human_gap_zoom",        "human",  "chr21", 10_672_000, 10_687_000,
     "15 kb inside the chr21 window BWTandem calls and adotto does not"),
]

# Tighter edge windows: the 30 kb pair above resolves the boundary but not individual
# calls. 16 kb centred on the annotated array start keeps both readable.
LOCI += [
    ("arabidopsis_cen1_edge16", "colcen", "Chr1", 14_369_200, 14_385_200,
     "16 kb centred on the Chr1 centromere start, 14,377,192"),
    ("maize_knob180_edge16",    "maize",  "CM039155.1", 4_635_850, 4_651_850,
     "16 kb centred on the knob180 array start, 4,643,853"),
    ("human_adotto_edge16",     "human",  "chr21", 10_730_700, 10_746_700,
     "16 kb centred on the first adotto region after the uncovered window, 10,738,710"),
]


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"FAILED: {' '.join(cmd[:4])}...\n{r.stderr[:600]}")
    return r.stdout


def _short(motif, keep=24):
    """Plots never render a 500 bp consensus; carrying it made the Chr1 track file
    4.4 MB. The full length is kept in its own column so nothing is lost."""
    return motif if len(motif) <= keep else motif[:keep] + "..."


def period_of(parts):
    """int(), never int(float()): column 5 is a period in converted competitor BEDs
    but a copy count like '46.3' in a BED BWTandem wrote, where the period is
    len(motif). int(float('46.3')) silently yields 46 and turns a period-banded
    metric into a copy-count one."""
    if len(parts) >= 5:
        try:
            return int(parts[4])
        except ValueError:
            pass
    return len(parts[3]) if len(parts) >= 4 else None


def slice_bed(path, chrom, start, end):
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or p[0] != chrom:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            if e <= start or s >= end:
                continue
            out.append((s, e, p[3] if len(p) > 3 else "", period_of(p)))
    return sorted(out)


def pick_human_window(width=80_000):
    """Densest `width` window of BWTandem calls that no adotto region overlaps."""
    tool = dict(GENOMES["human"]["tools"])["BWTandem"]
    ref = GENOMES["human"]["reference"][0]
    tmp = os.path.join(OUT, "_tmp_uncovered.bed")
    with open(tmp, "w") as fh:
        subprocess.run([BEDTOOLS, "intersect", "-a", tool, "-b", ref, "-v"],
                       stdout=fh, stderr=subprocess.DEVNULL, check=False)
    best, counts = None, defaultdict(int)
    span = defaultdict(lambda: [10**12, 0])
    with open(tmp) as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) < 3:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            k = (p[0], s // width)
            counts[k] += 1
            span[k][0] = min(span[k][0], s)
            span[k][1] = max(span[k][1], e)
    os.remove(tmp)
    if not counts:
        return None
    (chrom, w), n = max(counts.items(), key=lambda kv: kv[1])
    return chrom, w * width, (w + 1) * width, n


def build_tracks():
    os.makedirs(OUT, exist_ok=True)
    index = [("locus", "genome", "chrom", "start", "end", "n_calls", "description")]
    for locus, genome, chrom, start, end, desc in LOCI:
        if chrom is None:
            picked = pick_human_window()
            if picked is None:
                print(f"  {locus}: no uncovered calls, skipped")
                continue
            chrom, start, end, n = picked
            desc = f"{desc} -> {chrom}:{start}-{end}, {n} uncovered BWTandem calls"
        cfg = GENOMES[genome]
        rows = [("tool", "chrom", "start", "end", "motif", "motif_len", "period")]
        ref_path, ref_label = cfg["reference"]
        for s, e, motif, per in slice_bed(ref_path, chrom, start, end):
            rows.append((ref_label, chrom, s, e, _short(motif), len(motif), per if per else ""))
        total = 0
        for name, path in cfg["tools"]:
            if not os.path.exists(path):
                print(f"  {locus}: MISSING {name} -> {path}")
                continue
            calls = slice_bed(path, chrom, start, end)
            total += len(calls)
            for s, e, motif, per in calls:
                rows.append((name, chrom, s, e, _short(motif), len(motif), per if per else ""))
        dst = os.path.join(OUT, f"tracks_{locus}.tsv")
        with open(dst, "w") as fh:
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        index.append((locus, genome, chrom, start, end, total, desc))
        print(f"  {locus:18s} {chrom}:{start}-{end}  {total} tool calls  -> {os.path.basename(dst)}")
    with open(os.path.join(OUT, "tracks_index.tsv"), "w") as fh:
        for r in index:
            fh.write("\t".join(str(x) for x in r) + "\n")


def build_upset(genome):
    os.makedirs(OUT, exist_ok=True)
    cfg = GENOMES[genome]
    names, paths, tmp = [], [], []
    for name, path in cfg["tools"]:
        if not os.path.exists(path):
            print(f"  UpSet {genome}: MISSING {name}, skipped")
            continue
        srt = os.path.join(OUT, f"_tmp_{genome}_{name}.bed")
        with open(srt, "w") as fh:
            subprocess.run(["sort", "-k1,1", "-k2,2n", path], stdout=fh, check=True)
        names.append(name)
        paths.append(srt)
        tmp.append(srt)
    seg = os.path.join(OUT, f"upset_{genome}_segments.bed.gz")
    combos, combo_bp = defaultdict(int), defaultdict(int)
    proc = subprocess.Popen([BEDTOOLS, "multiinter", "-header",
                             "-names", *names, "-i", *paths],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    with gzip.open(seg, "wt") as gz:
        for line in proc.stdout:
            gz.write(line)
            if line.startswith("chrom"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            key = tuple(sorted(f[4].split(",")))
            combos[key] += 1
            combo_bp[key] += int(f[2]) - int(f[1])
    proc.wait()
    for f in tmp:
        os.remove(f)
    if proc.returncode != 0:
        sys.exit(f"bedtools multiinter failed:\n{proc.stderr.read()[:600]}")
    dst = os.path.join(OUT, f"upset_{genome}_combinations.tsv")
    with open(dst, "w") as fh:
        fh.write("degree\ttools\tsegments\tbp\n")
        for key in sorted(combos, key=lambda k: (-combo_bp[k], k)):
            fh.write(f"{len(key)}\t{','.join(key)}\t{combos[key]}\t{combo_bp[key]}\n")
    print(f"  UpSet {genome:8s} {len(combos)} distinct tool sets over {sum(combos.values())} "
          f"segments -> {os.path.basename(dst)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", default="colcen",
                    help="colcen | maize | human | all (human needs a batch job)")
    ap.add_argument("--tracks-only", action="store_true")
    a = ap.parse_args()
    print("FIGURE T — tool tracks at representative loci")
    build_tracks()
    if not a.tracks_only:
        print("\nFIGURE U — UpSet combinations")
        for g in (["colcen", "maize", "human"] if a.genome == "all" else [a.genome]):
            build_upset(g)
