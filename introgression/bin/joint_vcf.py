#!/usr/bin/env python3
"""
Joint biallelic-SNP VCF in target coordinates.

Genotypes:
 - references: paftools projections (1/1 = carries the ALT, 0/0 = callable and reference,
   ./. = not callable or a different ALT)
 - target: Clair3 pileup calls (0/1 or 1/1 = PASS SNV, 0/0 = RefCall with GQ >= 8, ./. otherwise)
Only sites with one ALT across all samples are kept, and sites invariant among called samples
are dropped.

--set quartet: parent1, parent2, target, outgroups (for Dsuite and Twisst); parent1, parent2
               and the target must be called, plus at least one outgroup.
--set panel:   target then every reference in references.tsv order (for HyDe, PCA/ADMIXTURE);
               at least 4 samples called.
"""
import argparse, collections, sys
import pysam
from common import read_references, by_role, load_callable, contains, load_snv

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--set", choices=["quartet", "panel"], required=True)
ap.add_argument("--references", required=True)
ap.add_argument("--aln_dir", required=True)
ap.add_argument("--clair3_vcf", required=True, help="Clair3 pileup.vcf.gz of the target")
ap.add_argument("--target", required=True, help="target FASTA (indexed)")
ap.add_argument("--target_name", required=True)
ap.add_argument("--out", required=True, help="uncompressed VCF")
a = ap.parse_args()
GQ_MIN = 8

refs = read_references(a.references)
T = a.target_name
if a.set == "quartet":
    samples = by_role(refs, "parent1") + by_role(refs, "parent2") + [T] + by_role(refs, "outgroup")
else:
    samples = [T] + [r["name"] for r in refs]
asm = [s for s in samples if s != T]

asm_snv, asm_call = {}, {}
for name in asm:
    asm_snv[name] = load_snv(f"{a.aln_dir}/{name}.vcf.gz")
    asm_call[name] = load_callable(f"{a.aln_dir}/{name}.callable.bed")
    print(f"  {name}: {sum(len(x) for x in asm_snv[name].values())} SNV", file=sys.stderr)

t_snv = collections.defaultdict(dict)   # chrom -> pos -> (alt, gt)
t_ref = collections.defaultdict(set)    # confident hom-ref positions
for r in pysam.VariantFile(a.clair3_vcf):
    s = r.samples[0]
    gt = s.get("GT")
    flt = list(r.filter)
    if flt == ["PASS"] and r.alts and len(r.alts) == 1 and len(r.ref) == 1 and len(r.alts[0]) == 1:
        if gt in ((0, 1), (1, 0)):
            t_snv[r.chrom][r.pos] = (r.alts[0].upper(), "0/1")
        elif gt == (1, 1):
            t_snv[r.chrom][r.pos] = (r.alts[0].upper(), "1/1")
    elif "RefCall" in flt and (s.get("GQ") or 0) >= GQ_MIN:
        t_ref[r.chrom].add(r.pos)
print(f"  {T}: {sum(len(x) for x in t_snv.values())} SNV, "
      f"{sum(len(x) for x in t_ref.values())} confident reference", file=sys.stderr)


def target_gt(c, p, alt):
    v = t_snv.get(c, {}).get(p)
    if v:
        return v[1] if v[0] == alt else "./."
    return "0/0" if p in t_ref.get(c, ()) else "./."


def asm_gt(n, c, p, alt):
    b = asm_snv[n].get(c, {}).get(p)
    if b is not None:
        return "1/1" if b == alt else "./."
    iv = asm_call[n].get(c)
    return "0/0" if (iv and contains(iv, p - 1)) else "./."


def keep(gts):
    called = [g for g in gts if g != "./."]
    if len(set(called)) <= 1:
        return False
    if a.set == "quartet":
        return all(g != "./." for g in gts[:3]) and any(g != "./." for g in gts[3:])
    return len(called) >= 4


fa = pysam.FastaFile(a.target)
n_out = 0
with open(a.out, "w") as o:
    o.write("##fileformat=VCFv4.2\n")
    for c, ln in zip(fa.references, fa.lengths):
        o.write(f"##contig=<ID={c},length={ln}>\n")
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
    for c in fa.references:
        cand = set(t_snv.get(c, {}))
        for n in asm:
            cand |= set(asm_snv[n].get(c, {}))
        for p in sorted(cand):
            alts = {asm_snv[n][c][p] for n in asm if p in asm_snv[n].get(c, {})}
            if p in t_snv.get(c, {}):
                alts.add(t_snv[c][p][0])
            if len(alts) != 1:
                continue
            alt = alts.pop()
            ref = fa.fetch(c, p - 1, p).upper()
            if ref not in "ACGT" or ref == alt:
                continue
            gts = [target_gt(c, p, alt) if s == T else asm_gt(s, c, p, alt) for s in samples]
            if not keep(gts):
                continue
            o.write(f"{c}\t{p}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t" + "\t".join(gts) + "\n")
            n_out += 1
print(f"wrote {n_out} sites ({a.set}: {', '.join(samples)}) -> {a.out}", file=sys.stderr)
