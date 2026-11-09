# porites-sp-indo-ms

Analysis and manuscript repository for the *Porites* sp. Indonesia genome assembly and annotation.

The assembly and annotation pipeline is a separate repository (**GenomeMaker**) and stays independent. This
repository holds the assembly-specific analyses of the manuscript; they read the pipeline results
(`results_v1.2`) through untracked links in `data/`.

| Folder | Question | Method | Run |
|---|---|---|---|
| `comparison/` | how does the assembly and annotation compare with other *Porites*? | re-annotation of the NCBI genomes with the pipeline's stage 08, `funannotate compare` | `comparison/run_compare.sh` (Nextflow `main.nf`) |
| `species_tree/` | where does the new genome sit among *Porites*? | BUSCO single-copy proteins: concatenated ML (IQ-TREE) and coalescent (ASTRAL-IV) | `nextflow run species_tree/main.nf` |
| `introgression/` | is it a *P. lobata* x *P. harrisoni* hybrid? | ONT phasing and ancestry painting, D-statistics, Twisst, HyDe, PCA/ADMIXTURE, mitogenome tree | `nextflow run introgression/main.nf -params-file ...` |
| `genome-assembly-annotation/` | the pipeline | link to GenomeMaker `results_v1.2` | |
| `data/` | inputs outside the repo | untracked links, see `data/README.md` | |
| `docs/` | repository notes | `docs/CLEANUP.md`: what is redundant | |

Each analysis folder has a `README.md` with what it needs, how to run it, its parameters and its results.
The introgression tests take (*P. lobata*, *P. harrisoni*) + the new genome as the species tree, from the
outgroup-rooted BUSCO trees of `species_tree/` (outgroups *P. cylindrica*, *P. australiensis*,
*P. divaricata*). `comparison/` is independent of the other two.

## Setup

1. Link the pipeline results and the other inputs as in `data/README.md`.
2. Tools: Nextflow, Apptainer (introgression, comparison), conda env `phylo` (species_tree), conda env
   `ncbi_datasets` (comparison download).
3. Large outputs (`results/`, `refs/`, `compare_*/`, `work/`, genomes, alignments) are written inside the
   analysis folders and are git-ignored; small trees, tables and reports under `archive/` are tracked.

## Layout rules

- Pipeline code is never edited here; `comparison/` imports it from `data/GenomeMaker`.
- Earlier runs on superseded assemblies (Canu, Flye) are in each folder's `archive/`.
- No host-absolute paths in code: paths are params relative to `data/`.
