#!/usr/bin/env python3
"""fig5_array_structure — self-similarity of one satellite array per species.

Left of each row: the pairwise-identity matrix of the array against itself, drawn
as the upper triangle rotated 45 degrees, so the x axis stays in genome
coordinates and height is the separation between the two windows being compared.
Higher-order structure reads as bands parallel to the baseline: a band at height
h means the array repeats itself every h bases.

Right of each row: a self dot plot over the same interval, which shows the same
structure as lines rather than as colour and keeps the monomer period visible.

Windows are compared by shared canonical k-mers (containment), not by alignment,
so an indel between two copies costs a few k-mers instead of shifting everything
after it.

Regions come from BWTandem's own deposited calls; see REGIONS below for the
coordinates and the reported motif length of each.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import figstyle

figstyle.setup()

K_IDENT = 13          # k-mer for the identity matrix
K_DOT = 21            # k-mer for the dot plot; longer, so the plot stays sparse
WIN = 1000            # identity window
STEP = 250            # identity step
DOT_TARGET = 25_000   # aim for this many plotted pairs; stride adapts per region
PROBE_STRIDE = 32     # cheap first pass used only to calibrate the stride

REGIONS = [
    # label, fasta, sequence name, start, end, reported motif length
    ("Human GRCh38 — chr1 α-satellite",
     "/data/gpfs/assoc/pgl/devel/exp1_human/data/hg38_primary.fa",
     "chr1", 122_257_803, 122_317_803, 1866),
    ("Arabidopsis Col-CEN — Chr4 CEN180",
     "/data/gpfs/assoc/pgl/filip/bwt/colcen/colcen.fa",
     "Chr4", 4_985_644, 5_045_644, 178),
    ("Maize Mo17 — CentC array",
     "/data/gpfs/assoc/pgl/devel/exp1_human/tools2026/fasta/zmays.fna",
     "CM039157.1", 53_652_014, 53_678_133, 155),
]


def read_region(path, name, start, end):
    """Pull one interval out of a FASTA without loading the whole file."""
    seq, grabbing, kept = [], False, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if grabbing:
                    break
                grabbing = line[1:].split()[0] == name
                pos = 0
                continue
            if not grabbing:
                continue
            n = len(line) - 1
            if pos + n > start and pos < end:
                lo, hi = max(0, start - pos), min(n, end - pos)
                seq.append(line[lo:hi])
                kept += hi - lo
            pos += n
            if pos >= end:
                break
    return "".join(seq).upper()


def canonical_kmers(s, k):
    """Set of canonical k-mers, as 2-bit packed ints. Windows containing N are
    simply shorter; they are not padded, so identity stays a real fraction."""
    code = {"A": 0, "C": 1, "G": 2, "T": 3}
    out, fwd, rev, valid = set(), 0, 0, 0
    mask = (1 << (2 * k)) - 1
    for ch in s:
        c = code.get(ch)
        if c is None:
            fwd = rev = valid = 0
            continue
        fwd = ((fwd << 2) | c) & mask
        rev = (rev >> 2) | ((3 - c) << (2 * (k - 1)))
        valid += 1
        if valid >= k:
            out.add(min(fwd, rev))
    return out


def identity_matrix(seq, win, step, k):
    starts = list(range(0, max(1, len(seq) - win + 1), step))
    sets = [canonical_kmers(seq[s:s + win], k) for s in starts]
    n = len(sets)
    m = np.full((n, n), np.nan)
    for i in range(n):
        a = sets[i]
        if not a:
            continue
        for j in range(i, n):
            b = sets[j]
            if not b:
                continue
            # containment against the smaller set: robust to one window
            # carrying a partial copy at the array edge
            m[i, j] = m[j, i] = len(a & b) / min(len(a), len(b))
    return np.array(starts), m


def triangle(ax, starts, m, origin, title):
    """Upper triangle rotated 45 degrees: x is the midpoint of the compared pair
    in genome coordinates, y is half their separation."""
    verts, vals = [], []
    step = starts[1] - starts[0] if len(starts) > 1 else 1
    h = step / 2.0
    for i in range(len(starts)):
        for j in range(i, len(starts)):
            v = m[i, j]
            if not np.isfinite(v):
                continue
            cx = ((starts[i] + starts[j]) / 2.0 + origin) / 1e6
            cy = (starts[j] - starts[i]) / 2.0
            hx = h / 1e6
            verts.append([(cx - hx, cy), (cx, cy + h), (cx + hx, cy), (cx, cy - h)])
            vals.append(v)
    pc = PolyCollection(verts, array=np.array(vals), cmap="magma",
                        edgecolors="none", rasterized=True)
    pc.set_clim(0.0, 1.0)
    ax.add_collection(pc)
    ax.set_xlim(origin / 1e6, (origin + starts[-1] + step) / 1e6)
    ax.set_ylim(0, (starts[-1] / 2.0) * 1.02)
    ax.set_title(title, fontweight="bold", loc="left")
    ax.set_ylabel("separation / 2 (bp)")
    return pc


def dotplot(ax, seq, origin, k, stride, motif_len, mb):
    pos = {}
    code = {"A": 0, "C": 1, "G": 2, "T": 3}
    fwd, valid, mask = 0, 0, (1 << (2 * k)) - 1
    xs, ys = [], []
    for idx, ch in enumerate(seq):
        c = code.get(ch)
        if c is None:
            fwd = valid = 0
            continue
        fwd = ((fwd << 2) | c) & mask
        valid += 1
        if valid < k or (idx % stride):
            continue
        start = idx - k + 1
        for prev in pos.get(fwd, ()):
            xs.append((start + origin) / 1e6)
            ys.append((prev + origin) / 1e6)
        pos.setdefault(fwd, []).append(start)
    ax.scatter(xs, ys, s=0.5, c=figstyle.BWT, alpha=0.30, linewidths=0,
               rasterized=True)
    ax.set_title(f"self dot plot — k={k}, every {stride} bp\nreported motif {motif_len:,} bp",
                 fontweight="bold", loc="left", fontsize=11)
    ax.set_aspect("equal", adjustable="box")
    return len(xs)


def _null_ax():
    """Throwaway axes for the density probe; never drawn."""
    f = plt.figure()
    a = f.add_subplot(111)
    plt.close(f)
    return a


fig, axes = plt.subplots(len(REGIONS), 2, figsize=(15, 4.2 * len(REGIONS)),
                         gridspec_kw={"width_ratios": [1.35, 1]})

mappable = None
for row, (label, fasta, name, start, end, motif) in enumerate(REGIONS):
    seq = read_region(fasta, name, start, end)
    starts, m = identity_matrix(seq, WIN, STEP, K_IDENT)
    mappable = triangle(axes[row, 0], starts, m, start,
                        f"{label} — {len(seq)/1000:.0f} kb")
    axes[row, 0].set_xlabel(f"{name} position (Mb)")
    # one cheap probe at a coarse stride, then scale so every panel carries a
    # comparable number of points instead of one being a solid block
    probe = dotplot(_null_ax(), seq, start, K_DOT, PROBE_STRIDE, motif, True)
    # pairs scale as (1/stride)^2, so scale the stride by the square root of the
    # overshoot rather than linearly
    stride = int(round(PROBE_STRIDE * (probe / max(1, DOT_TARGET)) ** 0.5))
    stride = max(2, min(256, stride))
    n_dots = dotplot(axes[row, 1], seq, start, K_DOT, stride, motif, True)
    axes[row, 1].set_xlabel(f"{name} position (Mb)")
    axes[row, 1].set_ylabel(f"{name} position (Mb)")
    print(f"{label}: {len(seq):,} bp, {len(starts)} windows, {n_dots:,} dot pairs")

# A dedicated axes on the far right: attaching the bar to the axes list put it
# in the whitespace left by the square dot plots, on top of them.
fig.subplots_adjust(right=0.88)
cax = fig.add_axes([0.935, 0.34, 0.011, 0.32])
cb = fig.colorbar(mappable, cax=cax)
cb.set_label(f"shared canonical {K_IDENT}-mers (containment)")

caption = (
    r"$\mathbf{Figure\ 5.}$ "
    "Internal structure of one satellite array per species, taken from BWTandem's own calls. "
    "Left: pairwise window identity drawn as the upper triangle rotated 45 degrees, so x stays in genome "
    f"coordinates and height is half the separation between the compared windows ({WIN:,} bp windows, "
    f"{STEP} bp step, shared canonical {K_IDENT}-mers). Bands parallel to the baseline are higher-order "
    "repeats: a band at height h means the array repeats every 2h bases. Right: self dot plot over the same "
    f"interval at k={K_DOT}; the sampling stride is set per panel so each carries a comparable number of "
    "points. Identity is containment against the smaller "
    "window, so a window straddling an array edge is not penalised for being half empty."
)
plt.figtext(0.009, 0.02, caption, ha="left", va="top", fontsize=13, wrap=True)
plt.suptitle("Higher-order structure differs by species: human α-satellite is built from a multi-monomer "
             "unit, CEN180 and CentC repeat closer to the monomer.", y=0.995, fontsize=14)
fig.subplots_adjust(hspace=0.36, bottom=0.10, right=0.88)
figstyle.save(fig, "fig5_array_structure")
