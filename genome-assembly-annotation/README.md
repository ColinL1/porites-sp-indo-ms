# Genome assembly and annotation

The genome was assembled and annotated with [GenomeMaker](https://github.com/ColinL1/GenomeMaker) (a
separate repository). This folder holds the run config (`porites_binac2.config`) and the manuscript figures
and tables of the assembly and annotation; `results_v1.2` links to the pipeline's results (not tracked).

## How to run

```bash
ln -s <path>/GenomeMaker/results_v1.2 genome-assembly-annotation/results_v1.2
genome-assembly-annotation/plots/run_all.sh
```

Conda envs and containers used by the figure scripts are listed in the header of `plots/run_all.sh`.

## Outputs

```text
plots/
├── composite/       assembly figure: k-mer spectrum (GenomeScope2) and snail plot (BlobToolKit)
├── gene_support/    annotation figure: gene-model support, OMArk
├── supplement/      blob plot before contamination filtering
├── mitogenome/      circular map of the mitogenome
└── genome_circos/   circular overview of the genome
```
