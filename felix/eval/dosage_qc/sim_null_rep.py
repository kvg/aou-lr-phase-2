#!/usr/bin/env python3
"""2,000 null repeat loci (REF-relative length dosage, scale s) for calibration."""

import numpy as np
import pandas as pd

rng = np.random.default_rng(7)
ph = pd.read_csv("pheno.tsv", sep="\t")
ids, q = ph["ID"].tolist(), ph["q"].to_numpy()
N = len(ids)
support = np.arange(-3, 6)
hdr = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000000>\n"
    '##FORMAT=<ID=DS1,Number=1,Type=Float,Description="Dosage of ANC1">\n'
    '##FORMAT=<ID=DS2,Number=1,Type=Float,Description="Dosage of ANC2">\n'
    '##FORMAT=<ID=ANC1,Number=1,Type=Float,Description="Ancestry count 1">\n'
    '##FORMAT=<ID=ANC2,Number=1,Type=Float,Description="Ancestry count 2">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n"
)
with open("rep_null.vcf", "w") as fh:
    fh.write(hdr)
    for j in range(2000):
        a1 = (rng.random(N) > q).astype(int)
        a2 = (rng.random(N) > q).astype(int)
        # Mostly-REF loci with ancestry-specific allele spectra; some skewed/rare.
        conc = 0.3 if j % 4 == 0 else 1.0
        p0 = rng.dirichlet(np.full(9, conc))
        p1 = rng.dirichlet(np.full(9, conc))
        p0[3] += 2.0
        p1[3] += 2.0  # index 3 = x 0 (REF)
        p0, p1 = p0 / p0.sum(), p1 / p1.sum()
        x1 = np.where(a1 == 0, rng.choice(support, N, p=p0), rng.choice(support, N, p=p1))
        x2 = np.where(a2 == 0, rng.choice(support, N, p=p0), rng.choice(support, N, p=p1))
        s = float(max(abs(x1).max(), abs(x2).max(), 1))
        ds = np.zeros((N, 2))
        anc = np.zeros((N, 2), dtype=int)
        for x, a in ((x1, a1), (x2, a2)):
            ds[np.arange(N), a] += x / s
            anc[np.arange(N), a] += 1
        cells = "\t".join(
            f"{ds[i, 0]:.6g}:{ds[i, 1]:.6g}:{anc[i, 0]}:{anc[i, 1]}" for i in range(N)
        )
        fh.write(f"chr1\t{1_000_000 + j * 100}\tn{j}\tA\t<RPT>\t.\tPASS\t.\tDS1:DS2:ANC1:ANC2\t{cells}\n")
print("ok")
