#!/usr/bin/env python3
"""
HyDe hybridization test, target as the putative hybrid of the parent1 and parent2 species.

DNA matrix from the panel VCF (heterozygous target sites as IUPAC codes). Sites kept where the
target and every parent2-species sample are called, >= 2 parent1-species samples are called
(all of them if fewer than 2) and >= 1 sample of an outgroup species is called; thinned to
<= 200,000 sites. run_hyde.py is then run once per non-parent species as the outgroup.
Outputs (outdir): seqs.txt, map.txt, hyde_<species>-out.txt / -out-filtered.txt, hyde_summary.tsv
"""
import argparse, os, subprocess, sys
import pysam
from common import read_references

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--vcf", required=True, help="panel VCF (bgzipped)")
ap.add_argument("--references", required=True)
ap.add_argument("--target_name", required=True)
ap.add_argument("--outdir", required=True)
a = ap.parse_args()

refs = read_references(a.references)
T = a.target_name
POP = {T: T, **{r["name"]: r["species"] for r in refs}}
SAMP = list(POP)
sp1 = next(r["species"] for r in refs if r["role"] == "parent1")
sp2 = next(r["species"] for r in refs if r["role"] == "parent2")
og_species = {r["species"] for r in refs if r["role"] == "outgroup"}
p1_samples = [s for s in SAMP if POP[s] == sp1]
p2_samples = [s for s in SAMP if POP[s] == sp2]
og_samples = [s for s in SAMP if POP[s] in og_species]
min_p1 = min(2, len(p1_samples))
IUPAC = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("AC"): "M",
         frozenset("GT"): "K", frozenset("AT"): "W", frozenset("CG"): "S"}

seqs = {s: [] for s in SAMP}
n = 0
for rec in pysam.VariantFile(a.vcf):
    if rec.alts is None or len(rec.alts) != 1:
        continue
    ref, alt = rec.ref.upper(), rec.alts[0].upper()
    if ref not in "ACGT" or alt not in "ACGT":
        continue
    row = {}
    for s in SAMP:
        gt = rec.samples[s]["GT"]
        if gt is None or gt[0] is None:
            row[s] = "N"
        elif gt[0] == gt[1]:
            row[s] = ref if gt[0] == 0 else alt
        else:
            row[s] = IUPAC[frozenset(ref + alt)]
    if row[T] == "N" or any(row[s] == "N" for s in p2_samples):
        continue
    if sum(row[s] != "N" for s in p1_samples) < min_p1:
        continue
    if all(row[s] == "N" for s in og_samples):
        continue
    for s in SAMP:
        seqs[s].append(row[s])
    n += 1
step = max(1, n // 200000)
nsub = len(range(0, n, step))
print(f"{n} sites, {nsub} after thinning", file=sys.stderr)

od = os.path.abspath(a.outdir)
with open(f"{od}/seqs.txt", "w") as o:
    for s in SAMP:
        o.write(s + "\t" + "".join(seqs[s][i] for i in range(0, n, step)) + "\n")
with open(f"{od}/map.txt", "w") as o:
    for s in SAMP:
        o.write(f"{s}\t{POP[s]}\n")

summary = []
for og in dict.fromkeys(POP[s] for s in SAMP if POP[s] not in (T, sp1, sp2)):
    pref = f"{od}/hyde_{og}"
    cmd = ["run_hyde.py", "-i", f"{od}/seqs.txt", "-m", f"{od}/map.txt", "-o", og,
           "-n", str(len(SAMP)), "-t", str(len(set(POP.values()))), "-s", str(nsub), "--prefix", pref]
    print(">>", " ".join(cmd), file=sys.stderr)
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=od)
    sys.stderr.write(r.stdout[-2000:] + r.stderr[-2000:])
    with open(f"{pref}-out.txt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if f["Hybrid"] == T and {f["P1"], f["P2"]} == {sp1, sp2}:
                summary.append((og, f["P1"], f["Hybrid"], f["P2"], f["Zscore"], f["Pvalue"], f["Gamma"]))
with open(f"{od}/hyde_summary.tsv", "w") as o:
    o.write("outgroup\tP1\tHybrid\tP2\tZscore\tPvalue\tGamma\n")
    for r in summary:
        o.write("\t".join(r) + "\n")
print(open(f"{od}/hyde_summary.tsv").read())
