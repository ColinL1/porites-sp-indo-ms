# Species tree: placement of the new genome among *Porites*

Four trees: two taxon sets (`small`, `large`) x two approaches (concatenated ML, ASTRAL-IV),
plus quartet-support and tanglegram plots.

```bash
cd species_tree
nextflow run main.nf                                   # BUSCO on ../data/results_v1.2/final/genome/Porites_sp.fa
nextflow run main.nf --new_fasta <assembly .fa>        # another assembly
nextflow run main.nf --new_busco_dir <BUSCO genome-mode dir>    # or reuse a BUSCO run
nextflow run main.nf -stub --outdir /scratch/x -w /scratch/x/work    # wiring test, no tools
```

## What is needed

| Input | Param | Default |
|---|---|---|
| BUSCO runs of the reference genomes (`busco_<name>/`, BUSCO 6.1.0, anthozoa_odb12), names in `references.tsv` (`@BUSCO_ROOT@`) | `busco_root` | `../data/busco_refs` |
| BUSCO lineage datasets (offline) | `busco_download_path` | `../data/busco_downloads` |
| Conda env `phylo`: BUSCO, MAFFT, trimAl, AMAS, IQ-TREE 3.1.3, ASTER, R (ape, phytools, phangorn) | `phylo_env` | `~/miniconda3/envs/phylo` |
| New assembly (BUSCO is run on it) or its BUSCO output | `new_fasta` / `new_busco_dir` | final assembly of `results_v1.2` |

`--new_name` sets the tip label (default `Porites_sp_nov_indonesia`). Output goes to
`results/{small,large}/` (not tracked).

## Taxon sets (`references.tsv`)

Both are the sets of the earlier runs, rooted on *Acropora muricata* GCA_036669905.1. The new
assembly is added to both.

- `small` (10 taxa, as `assembly/genomes/single_ref/`): one genome per *Porites* species,
  *P. australiensis*, *P. cylindrica*, *P. divaricata*, *P. evermanni*, *P. harrisoni*,
  *P. lobata*, *P. lutea*, *P. rus*.
- `large` (17 taxa, as `assembly/genomes/astral_expanded/`): `small` plus *P. cylindrica*
  GCA_964035495.1, *P. divaricata* GCA_038088405.1 and GCA_964035725.1, *P. lobata*
  GCA_942486035.1, *P. lutea* GCA_958299795.1 and GCA_958299805.1, *P. rus* GCA_964035695.1.

The reference BUSCO runs are the existing ones in `assembly/genomes/busco_*` (BUSCO 6.1.0,
anthozoa_odb12).

## Steps

| Process | Tool | Settings |
|---|---|---|
| `BUSCO_NEW` | BUSCO 6.1.0 | genome mode, anthozoa_odb12, offline |
| `SHARED_LOCI` | `bin/shared_loci.py` | loci complete and single-copy in every genome of the set; protein sequences |
| `ALIGN_TRIM` | MAFFT 7.526, trimAl 1.5.1 | `--auto`; `-automated1` |
| `CONCAT_TREE` | AMAS 1.0, IQ-TREE 3.1.3 | one partition per locus; `-m MFP -bb 1000 -alrt 1000` |
| `GENE_TREES` | IQ-TREE 3.1.3 | `-S` (one tree per locus), `-m MFP` |
| `ASTRAL` | ASTRAL-IV (ASTER 1.25) | `-u 2 --root Amuricata_GCA_036669905.1`, unconstrained |
| `PLOT` | R (ape, phytools, phangorn), `bin/plot_trees.R` | quartet pies, tanglegram, normalised RF |

The earlier ASTRAL-IV runs passed the concatenated tree with `-c`, which in ASTRAL-IV is a
topology constraint; this pipeline does not.

All tools run from the local `phylo` conda env (`nextflow.config`). Tested with `-stub`, and
the loci, alignment, concatenation and plot steps on the existing BUSCO runs (small set with
`canu_lut_rag`: 1,765 loci).

## Earlier runs (`archive/`)

`archive/runs/01_concat_aa_9taxa_canu` and `02_concat_aa_10taxa_canu_flye` (concatenated ML on the
Canu and Flye trial assemblies) and `archive/context_flye/` (ASTRAL-IV and concatenated trees
with the Flye assembly, with the 17-taxon extension): trees, partitions and commands only; the
alignments are not tracked. Their overview is in `archive/README_earlier_runs.md` and
`archive/README_phylogeny_overview_2026-10-30.md`. `tables/busco_inputs_9taxa.tsv` lists the BUSCO
completeness of the inputs of the first run.
