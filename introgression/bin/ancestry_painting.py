#!/usr/bin/env python3
"""
Per-haplotype ancestry painting from the haplotagged ONT reads.

At every diagnostic site, pile up the reads of haplotype 1 and haplotype 2 (HP tag) and call
the base each haplotype carries (>= MIN_READS reads, >= MIN_FRAC agreement, base quality >= 7).
Each haplotype allele is classed P1 (parent1 allele), P2 (parent2 allele) or other, and the
site pattern is HOM_P1, HOM_P2, HET_opposite (one haplotype each: the F1 pattern), half_* (one
haplotype resolved) or other. Aggregated per phase block (>= 10 sites) and per 100 kb window.

Outputs (outdir): site_ancestry.tsv.gz, block_ancestry.tsv, window_ancestry.tsv,
ancestry_summary.txt, ancestry_painting.png
"""
import argparse, collections, gzip, sys
from multiprocessing import Pool
import numpy as np
import pysam
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import read_references, by_role

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--diag", required=True, help="diagnostic_sites.tsv.gz")
ap.add_argument("--bam", required=True, help="haplotagged BAM")
ap.add_argument("--fai", required=True)
ap.add_argument("--references", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--threads", type=int, default=32)
a = ap.parse_args()

WIN = 100_000
MIN_READS = 3      # per haplotype to make a call
MIN_FRAC = 0.75    # consensus fraction
MIN_BLOCK_SITES = 10
CHUNK = 1_000_000   # pileup jobs of 1 Mb, so all threads are used
refs = read_references(a.references)
P1, P2 = by_role(refs, "parent1")[0], by_role(refs, "parent2")[0]

diag = collections.defaultdict(dict)
with gzip.open(a.diag, "rt") as fh:
    next(fh)
    for line in fh:
        c, p, ref, b1, b2 = line.rstrip("\n").split("\t")
        diag[c][int(p)] = (ref, b1, b2)


def hap_base(counts):
    tot = sum(counts.values())
    if tot < MIN_READS:
        return None
    b, n = max(counts.items(), key=lambda x: x[1])
    return b if n / tot >= MIN_FRAC else None


def process_chunk(job):
    chrom, start, end = job
    sites = {p: v for p, v in diag[chrom].items() if start < p <= end}
    want = set(sites)
    bam = pysam.AlignmentFile(a.bam)
    h1c = collections.defaultdict(collections.Counter)
    h2c = collections.defaultdict(collections.Counter)
    psv = collections.defaultdict(collections.Counter)
    for col in bam.pileup(chrom, start, end, min_base_quality=7, truncate=True, stepper="samtools",
                          ignore_orphans=False, max_depth=200):
        pos1 = col.reference_pos + 1
        if pos1 not in want:
            continue
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            aln = pr.alignment
            if not aln.has_tag("HP"):
                continue
            hp = aln.get_tag("HP")
            base = aln.query_sequence[pr.query_position].upper()
            if base not in "ACGT":
                continue
            if hp == 1:
                h1c[pos1][base] += 1
            elif hp == 2:
                h2c[pos1][base] += 1
            if aln.has_tag("PS"):
                psv[pos1][aln.get_tag("PS")] += 1
    out = []
    for pos1, (ref, b1, b2) in sites.items():
        h1, h2 = hap_base(h1c.get(pos1, {})), hap_base(h2c.get(pos1, {}))
        if h1 is None and h2 is None:
            continue
        ps = psv[pos1].most_common(1)[0][0] if psv.get(pos1) else -1

        def cls(b):
            if b is None: return "?"
            if b == b1 and b != b2: return "P1"
            if b == b2 and b != b1: return "P2"
            return "o"
        c1, c2 = cls(h1), cls(h2)
        if c1 in ("P1", "P2") and c2 in ("P1", "P2"):
            pattern = "HET_opposite" if c1 != c2 else f"HOM_{c1}"
        elif "?" in (c1, c2):
            pattern = "half_" + (c1 if c2 == "?" else c2)
        else:
            pattern = "other"
        out.append((chrom, pos1, ref, b1, b2, h1 or ".", h2 or ".", c1, c2, ps, pattern))
    return out


def side(f):
    return "?" if f != f else ("P1" if f >= .8 else "P2" if f <= .2 else "mixed")


if __name__ == "__main__":
    jobs = [(c, s, s + CHUNK) for c in diag for s in range(0, max(diag[c]) + 1, CHUNK)]
    with Pool(a.threads) as pool:
        rows = sorted(r for sub in pool.map(process_chunk, jobs, chunksize=1) for r in sub)
    print(f"{len(rows)} diagnostic sites with a haplotype-resolved call", file=sys.stderr)
    od = a.outdir
    with gzip.open(f"{od}/site_ancestry.tsv.gz", "wt") as o:
        o.write("chrom\tpos\tref\tp1\tp2\th1\th2\th1_class\th2_class\tps\tpattern\n")
        for r in rows:
            o.write("\t".join(map(str, r)) + "\n")

    pat = collections.Counter(r[10] for r in rows)
    h1P1 = sum(r[7] == "P1" for r in rows); h1P2 = sum(r[7] == "P2" for r in rows)
    h2P1 = sum(r[8] == "P1" for r in rows); h2P2 = sum(r[8] == "P2" for r in rows)
    n1, n2 = h1P1 + h2P1, h1P2 + h2P2

    blk = collections.defaultdict(collections.Counter)
    for chrom, pos, ref, b1, b2, h1, h2, c1, c2, ps, pattern in rows:
        if ps == -1: continue
        b = blk[(chrom, ps)]
        b["n"] += 1
        b["start"] = pos if b["start"] == 0 else min(b["start"], pos)
        b["end"] = max(b["end"], pos)
        b["h1" + c1] += 1; b["h2" + c2] += 1; b[pattern] += 1
    bc = collections.Counter()
    with open(f"{od}/block_ancestry.tsv", "w") as o:
        o.write("chrom\tps\tstart\tend\tn\th1_P1\th1_P2\th2_P1\th2_P2\thet_opposite\thom_P1\thom_P2\th1_P1frac\th2_P1frac\tcall\n")
        for (chrom, ps), b in sorted(blk.items()):
            if b["n"] < MIN_BLOCK_SITES: continue
            a1, a2 = b["h1P1"] + b["h1P2"], b["h2P1"] + b["h2P2"]
            f1 = b["h1P1"] / a1 if a1 else float("nan"); f2 = b["h2P1"] / a2 if a2 else float("nan")
            s1, s2 = side(f1), side(f2)
            call = ("opposite" if {s1, s2} == {"P1", "P2"} else f"both_{s1}" if s1 == s2 and s1 in ("P1", "P2")
                    else f"{s1}/{s2}")
            bc[call] += 1
            o.write(f"{chrom}\t{ps}\t{b['start']}\t{b['end']}\t{b['n']}\t{b['h1P1']}\t{b['h1P2']}\t{b['h2P1']}\t{b['h2P2']}\t"
                    f"{b['HET_opposite']}\t{b['HOM_P1']}\t{b['HOM_P2']}\t{f1:.3f}\t{f2:.3f}\t{call}\n")

    win = collections.defaultdict(collections.Counter)
    for chrom, pos, ref, b1, b2, h1, h2, c1, c2, ps, pattern in rows:
        w = win[(chrom, pos // WIN)]
        for tag, c in (("h1", c1), ("h2", c2)):
            if c in ("P1", "P2"): w[tag + c] += 1
        w[pattern] += 1; w["n"] += 1
    with open(f"{od}/window_ancestry.tsv", "w") as o:
        o.write("chrom\twin_start\tn\th1_P1\th1_P2\th2_P1\th2_P2\thet_opposite\thom_P1\thom_P2\th1_P1frac\th2_P1frac\n")
        for (chrom, wi), w in sorted(win.items()):
            a1, a2 = w["h1P1"] + w["h1P2"], w["h2P1"] + w["h2P2"]
            f1 = w["h1P1"] / a1 if a1 else float("nan"); f2 = w["h2P1"] / a2 if a2 else float("nan")
            o.write(f"{chrom}\t{wi*WIN}\t{w['n']}\t{w['h1P1']}\t{w['h1P2']}\t{w['h2P1']}\t{w['h2P2']}\t"
                    f"{w['HET_opposite']}\t{w['HOM_P1']}\t{w['HOM_P2']}\t{f1:.3f}\t{f2:.3f}\n")

    with open(f"{od}/ancestry_summary.txt", "w") as o:
        o.write(f"parent1 (P1) = {P1}; parent2 (P2) = {P2}\n")
        o.write(f"diagnostic sites: {sum(len(v) for v in diag.values())}\n")
        o.write(f"sites with >= 1 haplotype-resolved call: {len(rows)}\n\nsite patterns:\n")
        for k, v in pat.most_common():
            o.write(f"  {k:14s} {v:9d}  ({100*v/len(rows):.1f}%)\n")
        o.write(f"\nallele matches genome-wide: P1={n1}  P2={n2}  P1_fraction={n1/(n1+n2):.4f}\n")
        o.write(f"  H1: P1={h1P1} P2={h1P2}  P1_fraction={h1P1/(h1P1+h1P2):.4f}\n")
        o.write(f"  H2: P1={h2P1} P2={h2P2}  P1_fraction={h2P1/(h2P1+h2P2):.4f}\n\n")
        o.write(f"phase blocks with >= {MIN_BLOCK_SITES} resolved sites: {sum(bc.values())}\n")
        for k, v in bc.most_common():
            o.write(f"  block call {k:12s}: {v}\n")
    print(open(f"{od}/ancestry_summary.txt").read())

    # stacked fraction of HOM_P1 / HOM_P2 / HET_opposite per window, largest scaffolds
    fai = [(l.split()[0], int(l.split()[1])) for l in open(a.fai)]
    present = {c for c, _ in win}
    big = [c for c, ln in sorted(fai, key=lambda x: -x[1]) if ln >= 3_000_000 and c in present][:30]
    wd = collections.defaultdict(list)
    for (chrom, wi), w in win.items():
        wd[chrom].append((wi * WIN, w["n"], w["HOM_P1"], w["HOM_P2"], w["HET_opposite"]))
    fig, axes = plt.subplots(len(big), 1, figsize=(15, 1.15 * len(big)), squeeze=False)
    for ax, chrom in zip(axes[:, 0], big):
        d = np.array(sorted(wd[chrom]), float)
        x, nn = d[:, 0] / 1e6, d[:, 1].clip(1)
        hp1, hp2, het = d[:, 2] / nn, d[:, 3] / nn, d[:, 4] / nn
        ax.fill_between(x, 0, hp1, color="#1f77b4", alpha=.7, step="mid", label=f"HOM {P1}")
        ax.fill_between(x, hp1, hp1 + hp2, color="#d62728", alpha=.7, step="mid", label=f"HOM {P2}")
        ax.fill_between(x, hp1 + hp2, hp1 + hp2 + het, color="#7f7f7f", alpha=.7, step="mid", label="HET opposite")
        ax.set_ylim(0, 1); ax.set_yticks([0, 1])
        ax.set_ylabel(chrom.replace("_RagTag", ""), fontsize=6, rotation=0, ha="right", va="center")
        ax.tick_params(labelsize=6)
    axes[0, 0].legend(fontsize=7, ncol=3, loc="upper right")
    axes[-1, 0].set_xlabel("position (Mb)")
    fig.tight_layout()
    fig.savefig(f"{od}/ancestry_painting.png", dpi=130)
