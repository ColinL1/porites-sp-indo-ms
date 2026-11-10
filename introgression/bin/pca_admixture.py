#!/usr/bin/env python3
"""
PCA and ADMIXTURE on the filtered panel VCF (biallelic SNPs, <= 10% missing, MAF >= 0.05).

PCA is computed directly (mean-imputed, standardised genotype matrix, SVD): plink2 --pca
refuses fewer than 50 samples. The target's Euclidean distance (PC1-4) to every species
centroid and to the parent1 cluster radius is reported. ADMIXTURE K = 2..5 with
cross-validation on all filtered SNPs (plink2 bed, chromosome codes set to 1).
Outputs (outdir): pca.tsv, pca.png, pca_distances.txt, ad.<K>.Q, admixture_K<K>.log,
admixture_cv.tsv
"""
import argparse, os, subprocess
import numpy as np
import pysam
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import read_references

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--vcf", required=True, help="filtered panel VCF")
ap.add_argument("--references", required=True)
ap.add_argument("--target_name", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--pca_only", action="store_true", help="skip ADMIXTURE (Nextflow runs it per K)")
a = ap.parse_args()

refs = read_references(a.references)
T = a.target_name
SP = {T: T, **{r["name"]: r["species"] for r in refs}}
sp1 = next(r["species"] for r in refs if r["role"] == "parent1")
vcf = os.path.abspath(a.vcf)
os.chdir(a.outdir)

vf = pysam.VariantFile(vcf)
samples = list(vf.header.samples)
G = []
for rec in vf:
    if rec.alts is None or len(rec.alts) != 1:
        continue
    row = []
    for s in samples:
        gt = rec.samples[s]["GT"]
        row.append(np.nan if gt is None or gt[0] is None else gt[0] + (gt[1] if len(gt) > 1 and gt[1] is not None else gt[0]))
    G.append(row)
G = np.array(G, float).T            # samples x sites, 0/1/2 dosage (assemblies are 0 or 2)
nan = np.where(np.isnan(G))
G[nan] = np.take(np.nanmean(G, axis=0), nan[1])
G = G[:, G.std(0) > 0]
Gs = (G - G.mean(0)) / G.std(0)
U, S, Vt = np.linalg.svd(Gs / np.sqrt(Gs.shape[1]), full_matrices=False)
pcs = U * S
varexp = S ** 2 / np.sum(S ** 2)
npc = min(6, pcs.shape[1])
with open("pca.tsv", "w") as o:
    o.write("sample\tspecies\t" + "\t".join(f"PC{i+1}" for i in range(npc)) + "\n")
    for i, s in enumerate(samples):
        o.write(f"{s}\t{SP[s]}\t" + "\t".join(f"{pcs[i, j]:.4f}" for j in range(npc)) + "\n")
    o.write("variance_explained\t\t" + "\t".join(f"{varexp[j]:.4f}" for j in range(npc)) + "\n")

species = list(dict.fromkeys(SP[s] for s in samples if s != T))
cmap = plt.get_cmap("tab10")
col = {sp: cmap(i % 10) for i, sp in enumerate(species)}
col[T] = "black"
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
for ax, (x, y) in zip(axes, [(0, 1), (0, 2)]):
    for i, s in enumerate(samples):
        sp = SP[s]
        ax.scatter(pcs[i, x], pcs[i, y], color=col[sp], s=140 if s == T else 90, marker="*" if s == T else "o",
                   edgecolor="k", zorder=3, label=sp if sp not in ax.get_legend_handles_labels()[1] else None)
        ax.annotate(s, (pcs[i, x], pcs[i, y]), fontsize=6, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel(f"PC{x+1} ({varexp[x]*100:.1f}%)"); ax.set_ylabel(f"PC{y+1} ({varexp[y]*100:.1f}%)")
    ax.axhline(0, lw=.3, c="grey"); ax.axvline(0, lw=.3, c="grey")
axes[1].legend(fontsize=8, loc="best")
fig.tight_layout(); fig.savefig("pca.png", dpi=140)

ti = samples.index(T)
cent = {sp: pcs[[i for i, s in enumerate(samples) if SP[s] == sp], :4].mean(0) for sp in species}
li = [i for i, s in enumerate(samples) if SP[s] == sp1]
radius = max(np.linalg.norm(pcs[i, :4] - cent[sp1]) for i in li)
d1 = np.linalg.norm(pcs[ti, :4] - cent[sp1])
with open("pca_distances.txt", "w") as o:
    o.write(f"{pcs.shape[0]} samples x {G.shape[1]} variable sites; variance explained PC1-5: "
            + " ".join(f"{v:.3f}" for v in varexp[:5]) + "\n")
    o.write(f"{T} Euclidean distance (PC1-4) to each species centroid:\n")
    for sp, c in sorted(cent.items(), key=lambda kv: np.linalg.norm(pcs[ti, :4] - kv[1])):
        o.write(f"  {sp:16s} {np.linalg.norm(pcs[ti, :4] - c):.3f}\n")
    o.write(f"{sp1} cluster radius (max distance to its centroid): {radius:.3f}; "
            f"{T} is {'inside' if d1 <= radius else 'outside'} it ({d1:.3f})\n")
print(open("pca_distances.txt").read())
if a.pca_only:
    raise SystemExit(0)

subprocess.run(f"plink2 --vcf {vcf} --set-all-var-ids @:# --allow-extra-chr --double-id --max-alleles 2 "
               f"--threads {a.threads} --make-bed --out ad", shell=True, check=True, capture_output=True)
with open("ad.bim") as f:
    lines = f.readlines()
with open("ad.bim", "w") as f:   # ADMIXTURE needs integer chromosome codes
    for l in lines:
        p = l.split("\t"); p[0] = "1"; f.write("\t".join(p))
with open("admixture_cv.tsv", "w") as cv:
    cv.write("K\tcv_error\n")
    for K in (2, 3, 4, 5):
        r = subprocess.run(f"admixture --cv ad.bed {K}", shell=True, capture_output=True, text=True)
        open(f"admixture_K{K}.log", "w").write(r.stdout + r.stderr)
        for line in r.stdout.splitlines():
            if "CV error" in line:
                cv.write(f"{K}\t{float(line.split(':')[1]):.5f}\n")
with open("ad.fam") as f:
    order = [l.split()[1] for l in f]
with open("admixture_samples.txt", "w") as f:
    f.write("\n".join(order) + "\n")
print(open("admixture_cv.tsv").read())
