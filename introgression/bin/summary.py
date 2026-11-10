#!/usr/bin/env python3
"""One row per test (summary/results.tsv), from the step folders of one run."""
import argparse, csv, glob, os, re
from common import read_references, by_role

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--outdir", required=True, help="results/<name>")
ap.add_argument("--references", required=True)
ap.add_argument("--target_name", required=True)
a = ap.parse_args()
O, T = a.outdir, a.target_name
refs = read_references(a.references)
P1, P2, OGS = by_role(refs, "parent1")[0], by_role(refs, "parent2")[0], by_role(refs, "outgroup")
rows = []


def add(test, result, folder):
    rows.append((test, result, folder))


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def safe(fn, test, folder):
    try:
        add(test, fn(), folder)
    except (OSError, KeyError, ValueError, IndexError, StopIteration) as e:
        add(test, f"not available ({type(e).__name__}: {e})", folder)


def het():
    d = dict(l.rstrip("\n").split("\t") for l in open(f"{O}/01_heterozygosity/heterozygosity.tsv"))
    return (f"{int(d['het_snvs']):,} heterozygous SNVs, {100*float(d['het_per_bp']):.3f}% of "
            f"{int(d['scaffold_bp'])/1e6:.1f} Mb; mean depth {float(d['mean_depth']):.2f}x")


def phasing():
    r = next(x for x in tsv(f"{O}/02_phasing/phase_block_stats.tsv") if x["chromosome"] == "ALL")
    return f"phase-block NG50 {int(r['block_n50'])/1e3:.1f} kb; {int(r['blocks']):,} blocks; {100*float(r['phased_fraction']):.1f}% of het SNVs phased"


def painting():
    s = open(f"{O}/04_ancestry_painting/ancestry_summary.txt").read()
    pct = dict(re.findall(r"^\s+(\S+)\s+\d+\s+\(([\d.]+)%\)", s, re.M))
    frac = re.search(r"P1_fraction=([\d.]+)", s).group(1)
    nres = re.search(r"haplotype-resolved call: (\d+)", s).group(1)
    blocks = re.search(r"resolved sites: (\d+)", s).group(1)
    opp = re.search(r"block call opposite\s*: (\d+)", s)
    return (f"of {int(nres):,} resolved diagnostic sites: {pct.get('HOM_P2', '0')}% homozygous {P2}, "
            f"{pct.get('HOM_P1', '0')}% homozygous {P1}, {pct.get('HET_opposite', '0')}% heterozygous opposite (F1 pattern); "
            f"{P1}-allele fraction {frac}; {opp.group(1) if opp else 0} of {int(blocks):,} blocks opposite-haplotype")


def dstats():
    parts = []
    for og in OGS:
        r = tsv(f"{O}/05_d_statistics/{og}/SETS_BBAA.txt")[0]
        parts.append(f"{og}: D={float(r['Dstatistic']):.3f} Z={float(r['Z-score']):.2f} f4-ratio={float(r['f4-ratio']):.4f}")
    return "; ".join(parts)


def twisst():
    return "; ".join(f"{r['outgroup']}: T3={float(r[f'T3_{P1}_{P2}']):.2f} T1={float(r[f'T1_{T}_{P1}']):.2f} "
                     f"T2={float(r[f'T2_{T}_{P2}']):.2f} (plurality {r['plurality']})"
                     for r in tsv(f"{O}/06_twisst/twisst_summary.tsv"))


def robustness():
    g = tsv(f"{O}/07_twisst_robustness/robustness_summary.tsv")
    n3 = sum(r["plurality"] == "T3" for r in g)
    s = open(f"{O}/07_twisst_robustness/context_summary.txt").read()
    corr = re.findall(r"^\S+\t([+-][\d.]+)\t([+-][\d.]+)\t([+-][\d.]+)$", s, re.M)
    rs = [abs(float(r)) for row in corr for r in row]
    return f"T3 plurality in {n3} of {len(g)} grid cells; max |r| context vs weight {max(rs):.3f}"


def hyde():
    return "; ".join(f"O={r['outgroup']}: Z={float(r['Zscore']):.1f} p={float(r['Pvalue']):.3g} gamma={float(r['Gamma']):.3f}"
                     for r in tsv(f"{O}/08_hyde/hyde_summary.tsv"))


def pca():
    s = open(f"{O}/09_pca_admixture/pca_distances.txt").read()
    near = re.findall(r"^\s+(\S+)\s+([\d.]+)$", s, re.M)[1:5]
    cv = tsv(f"{O}/09_pca_admixture/admixture_cv.tsv")
    best = min(cv, key=lambda r: float(r["cv_error"]))["K"]
    inside = re.search(r"is (inside|outside)", s).group(1)
    return (f"nearest centroids (PC1-4): " + ", ".join(f"{sp} {d}" for sp, d in near)
            + f"; {inside} the {next(r['species'] for r in refs if r['role'] == 'parent1')} cluster; ADMIXTURE best K = {best}")


def mito():
    src = tsv(f"{O}/10_mitogenome/mito_sources.tsv")
    used = sum(int(r["length"]) > 0 for r in src)
    last = open(f"{O}/10_mitogenome/mito_tree_summary.txt").read().strip().splitlines()[-1]
    tl = next(r["length"] for r in src if r["sample"] == T)
    return f"{T} mitogenome {int(tl):,} bp; {used} of {len(src)} samples in the tree; {last}"


for test, fn, folder in [("heterozygosity", het, "01_heterozygosity"), ("phasing", phasing, "02_phasing"),
                         ("ancestry painting", painting, "04_ancestry_painting"), ("D-statistics", dstats, "05_d_statistics"),
                         ("Twisst", twisst, "06_twisst"), ("Twisst robustness", robustness, "07_twisst_robustness"),
                         ("HyDe", hyde, "08_hyde"), ("PCA / ADMIXTURE", pca, "09_pca_admixture"),
                         ("mitogenome", mito, "10_mitogenome")]:
    safe(fn, test, folder)
os.makedirs(f"{O}/summary", exist_ok=True)
with open(f"{O}/summary/results.tsv", "w") as o:
    o.write("test\tresult\tfolder\n")
    for r in rows:
        o.write("\t".join(r) + "\n")
print(open(f"{O}/summary/results.tsv").read())
