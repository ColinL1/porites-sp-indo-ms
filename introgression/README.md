# Introgression: is the new assembly a *P. lobata* x *P. harrisoni* hybrid?

Hybridization / introgression tests for the final assembly (GenomeMaker `results_v1.2`),
as a Nextflow pipeline with containers. The target assembly is a parameter and the reference
panel is a table (`references.tsv`). What each step tests and how to read it: `docs/`
(start with `docs/overview.md`). The 2026-09-30 run on `canu_lut_rag` (the earlier
assembly, no evidence of hybridization) is kept in `archive/01_canu_lut_rag/`.

## What is needed

| Input | Where | Param |
|---|---|---|
| Final assembly | `../data/results_v1.2/final/genome/Porites_sp.fa` | `new_fasta` |
| Mitogenome | `../data/results_v1.2/final/mitogenome/mitogenome.fasta` | `new_mito` |
| ONT reads of the same individual | link as `../data/ont_reads.fastq.gz` | `reads` |
| Reference genomes (NCBI FASTA, names in `references.tsv`) | folder `../data/introgression_refs/` | `ref_dir` |
| AGP of the target (join test only) | must match the final scaffold names, otherwise leave `null` | `agp` |
| Containers | Apptainer; custom image built from `containers/introgression_tools.def` | `tools_container` |

`../data/` is not tracked (see `../data/README.md`). `@REF_DIR@` and `@MITO_DIR@` in
`references.tsv` are replaced by `ref_dir` and `mitogenomes/` when the run starts.

```bash
cd introgression
nextflow run main.nf -params-file params/Porites_sp_v1_2.yaml                       # real run
nextflow run main.nf -params-file params/Porites_sp_v1_2.yaml -profile stub -stub   # wiring test
```

Parameters are in `nextflow.config` and `params/*.yaml` (relative paths are taken from the
launch directory, so launch from `introgression/`). `-resume` reuses finished tasks.

Results go to `params.outdir` (default `../results/<new_name>_nf/`). The processes publish
nothing themselves: the `publish:` section of the workflow names what is kept, and the
`output {}` block at the end of `main.nf` says where each goes, in one place. To move a folder,
change its line there. Changing the layout doesn't change any task, so `-resume` republishes
without rerunning anything. Each step folder has the result tables at the top and the raw
or per-run files in subfolders:

```
00_inputs/               scaffolds.bed, references.tsv, params.json
01_heterozygosity/       heterozygosity.tsv, snvs.vcf.gz, reads.flagstat.txt, reads.meandepth.txt
  depth/                 mosdepth, 5 kb bins
02_phasing/              phase_stats.txt, phase_block_stats.tsv, phase_block_list.tsv, phased.vcf.gz
03_reference_alignment/  alignment_stats.tsv
  <reference>/           callable BED and VCF in target coordinates
04_ancestry_painting/    ancestry_summary.txt, tables, ancestry_painting.png
05_d_statistics/<og>/    Dsuite Dtrios per outgroup
06_twisst/               twisst_summary.tsv, plots
  <og>/                  trees, weights, windows per outgroup
07_twisst_robustness/    robustness_summary.tsv, context_summary.txt, context.tsv
  grid/<rep>_W<n>_<NJ|ML>/  trees, weights, windows per grid cell; grid/<rep>_groups.tsv
08_hyde/                 hyde_summary.tsv
  raw/                   HyDe input and output per outgroup species
09_pca_admixture/        pca_distances.txt, pca.tsv, pca.png, admixture_cv.tsv, admixture_samples.txt
  admixture/             Q matrix and log per K
10_mitogenome/           mito_tree_summary.txt, mito_sources.tsv, mito_seqs.fasta, mito_aln.fasta
  iqtree/                IQ-TREE output (mito_tree.*), mafft.log
summary/results.tsv      one row per test
_intermediate/           joint VCFs of 05 to 09 (BAMs stay in the work directory)
pipeline_info/           trace, report, timeline
```

Containers (Apptainer, images cached in `~/SINGULARITY_CACHE/`): BioContainers images at the
exact versions below for Clair3, WhatsHap, mosdepth, plink2, ADMIXTURE, Flye, MAFFT, IQ-TREE
and seqkit. Everything else runs in one custom image, `containers/introgression_tools.def`.
It has the tools without a BioContainers image: Dsuite, HyDe (phyde), Twisst and
genomics_general, pinned to the commits of the canu run. Twisst also gets the one local edit of the canu run's
checkout (`np.NaN` to `np.nan`, needed with NumPy 2). It also has the Python packages of
`bin/`, FastTree, and minimap2 / samtools / bcftools for the steps that pipe them together.
Build it once:

```bash
cd containers && apptainer build --fakeroot ~/SINGULARITY_CACHE/introgression_tools_1.0.sif introgression_tools.def
```

The parallel work is in Nextflow's hands: one task per reference
alignment, per scaffold (WhatsHap), per outgroup (Dsuite, Twisst), per robustness-grid cell
and per ADMIXTURE K. Apptainer mounts
`$HOME` automatically, so inputs elsewhere need a bind (`apptainer.runOptions`).

## Config of the earlier shell runner (`archive/shell_runner/`)

`run_all.sh` (conda envs, vendored Dsuite/Twisst checkouts that are not in this repo) is kept for
reference only; `main.nf` is the supported route. Its variables map to the lower-case params.

| Variable | Meaning | Default |
|---|---|---|
| `NEW_NAME` | sample label of the target (`[A-Za-z0-9_]+`) | required |
| `NEW_FASTA` | target assembly | required |
| `READS` | ONT reads of the same individual | required |
| `OUTDIR` | output folder | required |
| `AGP` | RagTag AGP of the target (join test in `07_twisst_robustness`) | none: join test skipped |
| `NEW_MITO` | target mitogenome | none: assembled from the reads with Flye |
| `SCAFFOLDS` / `MIN_SCAFFOLD_LEN` | scaffold set as a BED, or all sequences of at least this length | 10 Mb |
| `REFERENCES` | reference panel | `references.tsv` |
| `THREADS` | threads | 256 |
| `MITO_BAIT`, `MITO_ROOT` | mitogenome used to find mito contigs; root of the mito tree | *P. harrisoni* CM137701.1; `Pdivaricata2` |
| `ROBUSTNESS_OUTGROUP` | outgroup of the robustness grid | first outgroup |
| `ENV`, `CLAIR3_ENV`, `FLYE_CMD` | conda envs | `introgress`, `clair3`, Flye from `VTK2026_Assembly` |
| `DSUITE`, `TWISST`, `PARSEVCF` | the vendored checkouts in `assembly/Introgression/scripts/` | |

`archive/shell_runner/Porites_sp_v1_2.sh`: the final assembly. `archive/canu_lut_rag_rerun/`: the
canu run again, for comparison (its inputs are not part of this repo).

## Reference panel (`references.tsv`)

`name, fasta, preset, role, species, mito`. `preset` is the minimap2 preset of the
projection (`asm20` for divergent genomes). `role`: exactly one `parent1` and one `parent2`,
at least one `outgroup` (D-statistics and Twisst are run per outgroup), and `extra`
(replicate genomes, used only in HyDe, PCA/ADMIXTURE and the mito tree). `species` groups
samples for HyDe and the PCA centroids. `mito`: an optional mitogenome FASTA for that
reference; when empty, it is extracted from the assembly (`docs/10_mitogenome.md`). The
default panel leaves out the two alternate-haplotype assemblies of the canu run
(`docs/overview.md`); `archive/canu_lut_rag_rerun/references_canu_run.tsv` is the canu run's panel.

## Steps

| Step (results folder) | Script(s) | Tools | Settings |
|---|---|---|---|
| `00_inputs` | `main.nf` (INPUTS) | samtools 1.24 | scaffold set, copies of config and panel |
| `01_heterozygosity` | `map_reads.sh`, `call_snvs.sh` | minimap2 2.31, Clair3 2.0.3, bcftools 1.24, mosdepth 0.3.14 | `-ax map-ont`; model `r1041_e82_400bps_sup_v520`, `--min_coverage=4`; biallelic PASS SNVs, QUAL ≥ 10; depth in 5 kb bins |
| `02_phasing` | `phase.sh` | WhatsHap 2.8 | `phase` per scaffold, `stats`, `haplotag` |
| `03_reference_alignment` | `align_reference.sh` | minimap2 2.31 (`paftools.js`) | `-cx asm10` or `asm20`; `call -L 10000 -l 5000` |
| `04_ancestry_painting` | `diagnostic_sites.py`, `ancestry_painting.py` | pysam 0.24.1 | ≥ 3 reads, ≥ 75% agreement, base quality ≥ 7 per haplotype; blocks ≥ 10 sites; 100 kb windows |
| `05_d_statistics` | `joint_vcf.py --set quartet` | Dsuite 0.5 r58 `Dtrios` | one run per outgroup; RefCall GQ ≥ 8 = hom-ref |
| `06_twisst` | `twisst_trees.py`, `twisst_plot.py` | genomics_general 3924344 `parseVCF.py`, Twisst a5e1003 | 50-SNP windows, NJ with ties split evenly, target as two random pseudo-haplotypes (seed 1), `--method complete` |
| `07_twisst_robustness` | `twisst_robustness.py`, `robustness_summary.py`, `context.py` | FastTree 2.2.0, Twisst | W 50/100/200/500 x NJ/ML (`-nt -gtr`) x rand2/phased2/cons1; context of the rand2 W200 NJ windows |
| `08_hyde` | `joint_vcf.py --set panel`, `hyde.py` | phyde 1.0.2 `run_hyde.py` | ≤ 200,000 thinned sites; one run per non-parent species as outgroup |
| `09_pca_admixture` | `pca_admixture.py` | bcftools, numpy, plink2 2.0.0a6.9, ADMIXTURE 1.3.0 | biallelic, ≤ 10% missing, MAF ≥ 0.05; PCA by SVD; K = 2 to 5, `--cv` |
| `10_mitogenome` | `mito_sequences.py`, `mito_from_reads.sh`, `mito_tree_summary.py` | minimap2 `-x asm20`, Flye 2.9.6, MAFFT 7.525, IQ-TREE 3.1.3 | `--auto --adjustdirection`; `-m MFP -B 1000 -seed 1` |
| `summary` | `summary.py` | | `summary/results.tsv`, one row per test |

Tools come from the containers listed above (the shell runner used the conda envs `introgress` and `clair3`); Dsuite,
Twisst and genomics_general from the container (pinned commits).

## Changes from the canu run's scripts

- Paths, sample names and roles come from the config and `references.tsv`. The two
  joint-VCF scripts (`A4a`, `B1a`) are now one with `--set quartet|panel`.
- Twisst topologies are classified from the weights-file header (`common.twisst_weights`),
  so T1/T2/T3 cannot be rotated by column order or Newick spelling (the canu run fixed two
  bugs of that kind by hand).
- The steps run by hand in the canu run are scripted: scaffold set, *P. divaricata* `asm20`,
  `parseVCF.py`, the `twisst.py` runs, mosdepth.
- `A3b_paint.py` (painting from the phased VCF) is dropped: the reported painting was
  `A3c_paint_pileup.py`, kept here as `ancestry_painting.py`.
- PCA: the failing plink2 PCA on LD-pruned SNPs is dropped; the Python PCA and ADMIXTURE on
  all filtered SNPs (what the canu run reported) are kept.
- Mitogenome: the target can be given (`NEW_MITO`) and references can have one (`mito`
  column). Extracted reference contigs need ≥ 80% bait coverage (≥ 50% for a
  mitogenome-sized contig), and complete sequences are rotated to a common start.
- NJ: own implementation with exact distances that keeps every tied topology, weighted
  equally per window (`common.nj_all`, `twisst_average.py`). DendroPy's NJ broke ties by
  memory layout, so the canu run's Twisst weights were not reproducible (`docs/06_twisst.md`).
- Parallel work: the painting runs in 1 Mb chunks, the Twisst outgroups run at the same time,
  and NJ trees are built in a process pool.
- Figures have no titles.

## Tests

The canu run's read-based intermediates (Clair3 calls, phased VCF, haplotagged BAM,
reference projections) were linked into a scratch output folder, and every later step was
run with the shell runner (`archive/shell_runner/run_all.sh`, 2026-11-02, 256 threads, ~40 min) and compared with
`archive/01_canu_lut_rag/`, using the canu run's panel (`archive/canu_lut_rag_rerun/references_canu_run.tsv`):

| Test | Canu run | This code |
|---|---|---|
| heterozygosity | 2,892,125 het SNVs, 0.647% | identical |
| phasing | NG50 21.7 kb | identical |
| ancestry painting | 1,705,448 diagnostic sites, 1,498,807 resolved; 37.4 / 37.0 / 7.8%; 14 of 5,771 blocks | identical |
| D-statistics | D 0.024 / 0.021 / 0.009, Z 1.81 / 1.29 / 0.62 | identical |
| Twisst T3 (cyl / div / aus) | 0.431 / 0.713 / 0.440 (not reproducible, see `docs/06_twisst.md`) | 0.439 / 0.711 / 0.448 with ties split; T3 the plurality, T1 ≈ T2 |
| robustness | T3 the plurality; \|r\| ≤ 0.037 | T3 the plurality in 24 of 24 cells; \|r\| ≤ 0.032 |
| HyDe | p = 1.0, Z −23.5 / −152.1 / −33.9 / −113.5 | identical |
| PCA / ADMIXTURE | 1.318 to the *P. lobata* centroid, best K = 3 (CV 0.676) | identical |
| mitogenome | partial and nuclear-insert reference sequences | 10 complete NCBI mitogenomes (`mitogenomes/`) |

Tree building was run twice and gave identical files.

Shell vs Nextflow (2026-11-02): both versions were run end to end on a small real input, two
canu chromosomes (37.5 Mb) and the 462,740 reads that map to them, with the canu panel.
Every result table was identical. The only differences were file paths printed in WhatsHap's
table and the context header, and the row order of `alignment_stats.tsv` (now sorted the same
way). Then both versions on the full v1.2 input (2026-11-02; shell 2h40 including a restart, Nextflow
2h18 alongside it): `summary/results.tsv` identical, and every result table byte-identical,
except the mito tree. Before `-seed 1` was added, IQ-TREE's bootstrap differed between runs.
With it, topology and support are identical, and two branch lengths differ in the tenth decimal
(multithreaded floating point). Not tested: the Flye fallback for the target mitogenome (both
configs give `new_mito`).
