# Porites sp. Indonesia: genome assembly and annotation

Analyses for the manuscript:

> [Authors]. [Title: de novo genome assembly and annotation of *Porites* sp. Indonesia]. [Journal / preprint server], [year]. [DOI]

| | |
|---|---|
| Organism | *Porites* sp. Indonesia (NCBI taxid 3903096), isolate ID_GA_Por_15 |
| BioProject / BioSample | PRJNA1548328 / SAMN49911241 |
| Genome accession | [GenBank assembly accession, pending] |
| Assembly release | v1.2 |

The genome was assembled and annotated from Oxford Nanopore reads with the
[GenomeMaker](https://github.com/ColinL1/GenomeMaker) pipeline. This repository holds the analyses that
place the genome in context.

## Analyses

| Folder | Question | Approach |
|---|---|---|
| [comparison](comparison/) | How do the assembly and annotation compare with other *Porites*? | The annotated NCBI *Porites* genomes (and *Acropora millepora* as outgroup) are re-annotated with the same functional annotation as ours, then compared with `funannotate compare` |
| [species_tree](species_tree/) | Where does the new genome sit among *Porites* species? | Single-copy BUSCO proteins: concatenated maximum likelihood (IQ-TREE) and coalescent (ASTRAL-IV) trees, rooted on *Acropora muricata* |
| [introgression](introgression/) | Is it a hybrid of *P. lobata* and *P. harrisoni*? | Phased heterozygous variants and ancestry painting, D-statistics, topology weighting (Twisst), HyDe, PCA/ADMIXTURE and a mitogenome tree, with *P. cylindrica*, *P. australiensis* and *P. divaricata* as outgroups |

The introgression tests use the species tree (*P. lobata*, *P. harrisoni*) plus the new genome. The comparison is independent of the other two.

Each folder has its own README with inputs, how to run it and its parameters.
