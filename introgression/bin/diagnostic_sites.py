#!/usr/bin/env python3
"""
Sites where the two putative parents carry different alleles, in target coordinates.

A diagnostic site is a target position inside the paftools-callable regions of both parents
where the parents' alleles (alt from their projected SNV, else the target reference base)
differ. Biallelic SNVs only.

Output: chrom pos ref p1 p2 (p1 = parent1 allele, p2 = parent2 allele), gzipped TSV.
"""
import argparse, gzip, sys
import pysam
from common import read_references, by_role, load_callable, contains, load_snv

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--references", required=True)
ap.add_argument("--aln_dir", required=True, help="reference_alignment output (<name>.vcf.gz, <name>.callable.bed)")
ap.add_argument("--target", required=True, help="target FASTA (indexed)")
ap.add_argument("--out", required=True)
a = ap.parse_args()

refs = read_references(a.references)
p1, p2 = by_role(refs, "parent1")[0], by_role(refs, "parent2")[0]
print(f"parent1={p1} parent2={p2}", file=sys.stderr)
p1_call = load_callable(f"{a.aln_dir}/{p1}.callable.bed")
p2_call = load_callable(f"{a.aln_dir}/{p2}.callable.bed")
p1_snv = load_snv(f"{a.aln_dir}/{p1}.vcf.gz")
p2_snv = load_snv(f"{a.aln_dir}/{p2}.vcf.gz")
fa = pysam.FastaFile(a.target)

n = 0
with gzip.open(a.out, "wt") as out:
    out.write("chrom\tpos\tref\tp1\tp2\n")
    for chrom in sorted(set(p1_call) | set(p2_call)):
        c1, c2 = p1_call.get(chrom, []), p2_call.get(chrom, [])
        if not c1 or not c2:
            continue
        s1, s2 = p1_snv.get(chrom, {}), p2_snv.get(chrom, {})
        for pos in sorted(set(s1) | set(s2)):
            if not contains(c1, pos - 1) or not contains(c2, pos - 1):
                continue
            ref = fa.fetch(chrom, pos - 1, pos).upper()
            if ref not in "ACGT":
                continue
            b1, b2 = s1.get(pos, ref), s2.get(pos, ref)   # no SNV in a callable region = reference base
            if b1 == b2:
                continue
            out.write(f"{chrom}\t{pos}\t{ref}\t{b1}\t{b2}\n")
            n += 1
print(f"diagnostic sites: {n}", file=sys.stderr)
