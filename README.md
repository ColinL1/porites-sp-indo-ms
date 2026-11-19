# Porites sp. Indonesia: genome assembly and annotation

Analyses for the manuscript:

> [Authors]. [Title: de novo genome assembly and annotation of *Porites* sp. Indonesia]. [Journal / preprint server], [year]. [DOI]

| | |
| --- | --- |
| Organism | *Porites* sp. Indonesia (NCBI taxid 3903096), isolate ID_GA_Por_15 |
| BioProject / BioSample | PRJNA1548328 / SAMN49911241 |
| Genome accession | [GenBank assembly accession, pending] |
| Assembly release | v1.2 |

The genome was assembled and annotated from Oxford Nanopore reads using the
[GenomeMaker](https://github.com/ColinL1/GenomeMaker) pipeline; run settings and figure code are in
[genome-assembly-annotation](genome-assembly-annotation/). The analyses below place the assembly in a
comparative and evolutionary context.

## Analyses

| Folder | Aim | Methods |
| --- | --- | --- |
| [comparison](comparison/) | Comparative assessment of assembly and annotation quality against published *Porites* genomes | Annotated NCBI *Porites* genomes, with *Acropora millepora* as outgroup, re-annotated with the same functional annotation workflow and compared using `funannotate compare` |
| [species_tree](species_tree/) | Phylogenetic placement within *Porites* | Single-copy BUSCO orthologs analysed by concatenated maximum likelihood (IQ-TREE) and multispecies coalescent (ASTRAL-IV) inference, rooted on *Acropora muricata* |
| [introgression](introgression/) | Test for hybrid origin from *P. lobata* and *P. harrisoni* | Phased heterozygous variants with ancestry painting, D-statistics, topology weighting (Twisst), HyDe, PCA/ADMIXTURE and mitogenome phylogeny, using *P. cylindrica*, *P. australiensis* and *P. divaricata* as outgroups |

The introgression analyses use the *P. lobata* and *P. harrisoni* genomes from the species tree
dataset together with the new assembly. The comparison is independent of both.

Each folder contains its own README describing inputs, parameters and execution.

## Running

Requirements: Nextflow ≥ 25.04 and Apptainer (each process runs in a pinned container; the
introgression pipeline also needs its custom image, see its README), and the input links in
[data/](data/README.md). Each analysis is a Nextflow pipeline with its settings in a params file:

```bash
cd species_tree  && nextflow run main.nf -params-file params/Porites_sp_v1_2.yaml
cd introgression && nextflow run main.nf -params-file params/Porites_sp_v1_2.yaml
cd comparison    && nextflow run main.nf -params-file params/Porites_sp_v1_2.yaml
```

Each pipeline folder has the same layout:

| Path | Content |
| --- | --- |
| `main.nf` | workflow and results layout (`output {}` block) |
| `modules/`, `bin/` | processes; their scripts |
| `params/*.yaml` | run settings (inputs, output folder) |
| `nextflow.config` | defaults, profiles, error handling, run reports |
| `process_resources.config` | cpus, memory and time per process |
| `plots/` | manuscript figures and tables from the results (`plots/run_all.sh`) |
| `archive/` | earlier runs (not tracked) |

`shared/plots/setup.R` holds what all figure scripts share: names, colours and plot functions. Results,
`data/` and the `docs/` folders are not tracked.
