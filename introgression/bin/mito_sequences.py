#!/usr/bin/env python3
"""
One mitogenome sequence per sample, named by sample, for the mitochondrial tree.

Target: --target_mito (given, e.g. the pipeline's MitoHiFi mitogenome, or the Flye assembly
from mito-mapped reads, of which the contig with the best bait hit is taken).
References: the `mito` column of references.tsv when filled, else extracted from the
assembly: the contig with the best minimap2 -x asm20 hit of the bait; if that contig is
0.8-1.2x the bait length and >= 50% of the bait aligns to it, it is taken whole, otherwise the span of the hits on the same copy as
the best hit (same alignment diagonal) +- 300 bp (contigs can hold tandem copies or sit inside
nuclear scaffolds). A reference whose copy covers < 80% of the bait is skipped (no complete
mitogenome in the assembly).
Complete sequences (given, assembled or whole contigs) are oriented to the bait and rotated to
the bait's start; cut segments are oriented only. Outputs: mito_seqs.fasta, mito_sources.tsv
"""
import argparse, os, subprocess, sys
from multiprocessing import Pool
from common import read_references

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--references", required=True)
ap.add_argument("--target_name", required=True)
ap.add_argument("--target_mito", required=True)
ap.add_argument("--target_source", default="given")
ap.add_argument("--bait", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()
MARGIN = 300
MIN_BAIT_COV = 0.8
MIN_COV_WHOLE = 0.5  # a mitogenome-sized contig of a divergent species aligns only in part
DIAG_SLACK = 2000   # bp of indels allowed between hits on one copy
COMP = str.maketrans("ACGTRYKMSWNacgtrykmswn", "TGCAYRMKSWNtgcayrmkswn")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path, want=None):
    """{name: seq}; only `want` if given (streams large assemblies without an index)"""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (want is None or name == want):
                    seqs[name] = "".join(buf)
                    if want:
                        return seqs
                name, buf = line[1:].split()[0], []
            elif want is None or name == want:
                buf.append(line.strip())
    if name is not None and (want is None or name == want):
        seqs[name] = "".join(buf)
    return seqs


def hits(target_fa, threads=1):
    """PAF hits of the bait (query) on target_fa, best (longest block) first"""
    r = subprocess.run(["minimap2", "-x", "asm20", "-t", str(threads), target_fa, a.bait],
                       capture_output=True, text=True, check=True)
    h = [l.split("\t") for l in r.stdout.splitlines()]
    return sorted(h, key=lambda f: -int(f[10]))


def orient_rotate(seq, h):
    """rotate a complete mitogenome so it starts where the bait starts, on the bait's strand"""
    f = min((x for x in h if x[5] == h[0][5]), key=lambda x: int(x[2]))   # hit nearest bait start
    L, qs, ts, te = len(seq), int(f[2]), int(f[7]), int(f[8])
    if f[4] == "+":
        k = (ts - qs) % L
        return seq[k:] + seq[:k]
    p = (te - 1 + qs) % L
    r, k = rc(seq), L - 1 - p
    return r[k:] + r[:k]


def diag(x):
    """target offset of bait position 0 for a hit; equal for hits on the same copy"""
    return int(x[7]) - int(x[2]) if x[4] == "+" else int(x[8]) + int(x[2])


def from_assembly(name, fna):
    h = hits(fna)
    if not h:
        return name, None, "no bait hit"
    best = h[0]
    ctg, clen, bait_len = best[5], int(best[6]), int(best[1])
    copy = [x for x in h if x[5] == ctg and x[4] == best[4] and abs(diag(x) - diag(best)) <= DIAG_SLACK]
    covered = set()
    for x in copy:
        covered.update(range(int(x[2]), int(x[3])))
    cov = len(covered) / bait_len
    whole = 0.8 * bait_len <= clen <= 1.2 * bait_len   # mitogenome-sized contig
    if cov < (MIN_COV_WHOLE if whole else MIN_BAIT_COV):
        return name, None, f"skipped: hits on {ctg} ({clen} bp) cover {cov:.0%} of the bait (partial copy or NUMT)"
    seq = read_fasta(fna, ctg)[ctg]
    if whole:
        return name, orient_rotate(seq, [x for x in h if x[5] == ctg]), f"extracted {ctg} (whole contig, {cov:.0%} of bait)"
    s = max(0, min(int(x[7]) for x in copy) - MARGIN)
    e = min(clen, max(int(x[8]) for x in copy) + MARGIN)
    seq = seq[s:e]
    if best[4] == "-":
        seq = rc(seq)
    return name, seq, f"extracted {ctg}:{s+1}-{e} ({cov:.0%} of bait)"


def complete(name, path, source):
    seqs = read_fasta(path)
    tmp = f"{a.outdir}/.{name}.fa"
    with open(tmp, "w") as o:
        o.write("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    h = hits(tmp)
    os.remove(tmp)
    if not h:
        return name, None, f"{source}: no bait hit"
    seq = seqs[h[0][5]]
    return name, orient_rotate(seq, [x for x in h if x[5] == h[0][5]]), f"{source} {os.path.basename(path)}:{h[0][5]}"


refs = read_references(a.references)
jobs = [(complete, (a.target_name, a.target_mito, a.target_source))]
for r in refs:
    jobs.append((complete, (r["name"], r["mito"], "given")) if r["mito"] else (from_assembly, (r["name"], r["fasta"])))
with Pool(min(a.threads, len(jobs))) as pool:
    res = [pool.apply_async(f, args) for f, args in jobs]
    res = [x.get() for x in res]
with open(f"{a.outdir}/mito_seqs.fasta", "w") as o, open(f"{a.outdir}/mito_sources.tsv", "w") as t:
    t.write("sample\tsource\tlength\n")
    for name, seq, src in res:
        t.write(f"{name}\t{src}\t{len(seq) if seq else 0}\n")
        if seq:
            o.write(f">{name}\n{seq}\n")
        else:
            print(f"WARN {name}: {src}", file=sys.stderr)
print(open(f"{a.outdir}/mito_sources.tsv").read())
