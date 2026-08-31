# Figure outputs

Figure 1 is pending native period-100 job 6129408 at commit `0363d8b`.
`figure_curve_data.csv` is deliberately a non-numeric pending placeholder; it
must be replaced with the scored P/B/F/H values before a new active PNG/PDF is
generated.  Full-range identity-sweep job 6129410 separately supplies
Supplementary Table S3 and does not fill Figure 1.

After replacing the placeholder with the scored CSV, render the active figure:

```bash
python3 results/figures/make_curve_figure.py
```

The CSV must contain the columns `series,label,region_recall,region_precision,`
`bp_recall,bp_precision`, with BWTandem rows P/B/F/H and competitor rows ULTRA,
tantan, TRF and TRASH. The renderer validates the complete schema, labels and
0--100 numeric range before importing matplotlib or creating output files. With
the checked-in pending placeholder it exits non-zero with a clear message and
does not create `figure_curve.png` or `figure_curve.pdf`.

Files containing `superseded` preserve the withdrawn earlier-build curve for
audit history only.  They are not manuscript inputs.  The superseded renderer
writes only superseded filenames, so running it cannot recreate an apparently
active Figure 1.
