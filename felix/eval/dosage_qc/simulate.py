#!/usr/bin/env python3
"""Synthetic two-ancestry cohort for the M6 Option A spike.

Writes (in the current directory):
  grm.vcf             unphased markers for the step-1 GRM
  pheno.tsv           ID, x1, q, y_bin10, y_bin2, y_q
  joint.vcf           SNVs: phased GT + AN1/AN2 (FELIXla / FLARE layout)
  snv_admix.vcf       same SNVs as admixed dosage VCF (DS1, DS2, ANC1, ANC2)
  rep_s1.vcf          repeat loci: DS_k = sum of x_h over ancestry-k haps / s
  rep_s2.vcf          same loci with scale 2s
  rep_shift.vcf       same loci shifted so the shortest allele is 0, scale s
  rep_truth.tsv       per-locus scale/shift
"""

import numpy as np

rng = np.random.default_rng(20261002)
N, K = 2000, 2
ids = [f"s{i:04d}" for i in range(N)]
q = rng.beta(2, 2, N)  # proportion ancestry 0
x1 = rng.normal(size=N)


def lin(eta):
    return 1 / (1 + np.exp(-eta))


def thresh(liab, prev):
    return (liab > np.quantile(liab, 1 - prev)).astype(int)


liab = 1.0 * q + 0.5 * x1 + rng.logistic(size=N)
y_bin10 = thresh(liab, 0.10)
y_bin2 = thresh(liab, 0.02)
y_q = 0.5 * q + 0.3 * x1 + rng.normal(size=N)
with open("pheno.tsv", "w") as fh:
    fh.write("ID\tx1\tq\ty_bin10\ty_bin2\ty_q\n")
    for i in range(N):
        fh.write(f"{ids[i]}\t{x1[i]:.5f}\t{q[i]:.5f}\t{y_bin10[i]}\t{y_bin2[i]}\t{y_q[i]:.5f}\n")

hdr_samples = "\t".join(ids)

# GRM markers (chr2, unphased)
with open("grm.vcf", "w") as fh:
    fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr2,length=100000000>\n")
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{hdr_samples}\n")
    for j in range(1500):
        af = rng.uniform(0.05, 0.5)
        g = rng.binomial(2, af, N)
        gts = ["0/0", "0/1", "1/1"]
        fh.write(f"chr2\t{1000 + j * 1000}\tg{j}\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(gts[x] for x in g) + "\n")

# Local ancestry per haplotype per site: independent draws from q (no tracts needed here).
def draw_anc():
    a1 = (rng.random(N) > q).astype(int)  # 0 with prob q
    a2 = (rng.random(N) > q).astype(int)
    return a1, a2


snv_lines, adm_lines = [], []
for j in range(300):
    a1, a2 = draw_anc()
    if j < 60:  # rare: exercise SPA
        afs = rng.uniform(0.002, 0.01, K)
    else:
        afs = rng.uniform(0.02, 0.5, K)
    h1 = (rng.random(N) < afs[a1]).astype(int)
    h2 = (rng.random(N) < afs[a2]).astype(int)
    pos = 10_000 + j * 1000
    snv_lines.append(
        f"chr1\t{pos}\tv{j}\tA\tG\t.\tPASS\t.\tGT:AN1:AN2\t"
        + "\t".join(f"{h1[i]}|{h2[i]}:{a1[i]}:{a2[i]}" for i in range(N))
    )
    cells = []
    for i in range(N):
        ds = [0, 0]
        anc = [0, 0]
        for h, a in ((h1[i], a1[i]), (h2[i], a2[i])):
            ds[a] += h
            anc[a] += 1
        cells.append(f"{ds[0]}:{ds[1]}:{anc[0]}:{anc[1]}")
    adm_lines.append(f"chr1\t{pos}\tv{j}\tA\tG\t.\tPASS\t.\tDS1:DS2:ANC1:ANC2\t" + "\t".join(cells))

joint_hdr = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="Ancestry hap1">\n'
    '##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Ancestry hap2">\n'
    "##ANCESTRY=<anc0=0,anc1=1>\n"
    f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{hdr_samples}\n"
)
adm_hdr = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000000>\n"
    '##FORMAT=<ID=DS1,Number=1,Type=Float,Description="Dosage of ANC1">\n'
    '##FORMAT=<ID=DS2,Number=1,Type=Float,Description="Dosage of ANC2">\n'
    '##FORMAT=<ID=ANC1,Number=1,Type=Float,Description="Ancestry count 1">\n'
    '##FORMAT=<ID=ANC2,Number=1,Type=Float,Description="Ancestry count 2">\n'
    f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{hdr_samples}\n"
)
open("joint.vcf", "w").write(joint_hdr + "\n".join(snv_lines) + "\n")
open("snv_admix.vcf", "w").write(adm_hdr + "\n".join(adm_lines) + "\n")

# Repeat loci: x_h = repeat-unit deviation from REF (0 = REF), ancestry-specific spread,
# contractions included. All 40 loci are null.
rep = {"s1": [], "s2": [], "shift": []}
truth = ["locus\tscale\tshift\tmin_x\tmax_x"]
for j in range(40):
    a1, a2 = draw_anc()
    p0 = rng.dirichlet(np.ones(9))
    p1 = rng.dirichlet(np.ones(9))
    support = np.arange(-3, 6)
    x_h1 = np.where(a1 == 0, rng.choice(support, N, p=p0), rng.choice(support, N, p=p1))
    x_h2 = np.where(a2 == 0, rng.choice(support, N, p=p0), rng.choice(support, N, p=p1))
    s = float(max(abs(x_h1).max(), abs(x_h2).max()))
    mn = float(min(x_h1.min(), x_h2.min()))
    truth.append(f"r{j}\t{s}\t{mn}\t{mn}\t{max(x_h1.max(), x_h2.max())}")
    pos = 500_000 + j * 1000
    for key, scale, shift in (("s1", s, 0.0), ("s2", 2 * s, 0.0), ("shift", s - mn, mn)):
        # shift: x' = (x - min) / (max - min) keeps every per-hap value in [0, 1]
        sc = (max(x_h1.max(), x_h2.max()) - mn) if key == "shift" else scale
        cells = []
        for i in range(N):
            ds = [0.0, 0.0]
            anc = [0, 0]
            for x, a in ((x_h1[i], a1[i]), (x_h2[i], a2[i])):
                ds[a] += (x - shift) / sc
                anc[a] += 1
            cells.append(f"{ds[0]:.6g}:{ds[1]:.6g}:{anc[0]}:{anc[1]}")
        rep[key].append(f"chr1\t{pos}\tr{j}\tA\t<RPT>\t.\tPASS\t.\tDS1:DS2:ANC1:ANC2\t" + "\t".join(cells))
for key, lines in rep.items():
    open(f"rep_{key}.vcf", "w").write(adm_hdr + "\n".join(lines) + "\n")
open("rep_truth.tsv", "w").write("\n".join(truth) + "\n")
print("ok")
