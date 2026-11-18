# Comparison: the new genome against the annotated NCBI *Porites*

Assembly and annotation comparison with other *Porites* (and *Acropora millepora* as an outgroup).
Every genome in the params file, ours included, is reduced to bare gene models and re-annotated with
the same functional annotation (eggNOG-mapper, InterProScan, Phobius, `funannotate annotate`), so each
carries the same evidence; then `funannotate compare` runs once over all GenBank files. One Nextflow
run does all of it. The stage 07/08 modules are copies of GenomeMaker's (`modules/`), with the genome
name as an input; nothing is imported from GenomeMaker at run time.

## Run

```bash
cd comparison
nextflow run main.nf -profile local  -params-file params/Porites_sp_v1_2.yaml                # workstation (256 threads)
nextflow run main.nf -profile binac2 -params-file params/Porites_sp_v1_2_binac2.yaml         # BinAC 2 (SLURM + Apptainer)
nextflow run main.nf -profile stub -stub -params-file params/Porites_sp_v1_2.yaml            # wiring test
```

`-resume` reuses finished tasks. `-profile conda` runs the tools from the pinned conda packages
(`conda.useMamba`) instead of Apptainer. Per-process resources: `process_resource.config`.

## What is needed

| Input | Param (params file) | Where |
|---|---|---|
| Genomes | `genomes:` list of `{name, fasta, gff}` or `{name, accession}` | new genome: `../data/results_v1.2/release_curation/Porites_sp_v1.2/{genome,annotation}/Porites_sp._Indonesia_ID_GA_Por_15.{scaffolds.fa,gff3}`; NCBI genomes by accession (downloaded) |
| funannotate / eggNOG databases | `funannotate_db`, `eggnog_db` | local: GenomeMaker's `config/` (via `../data/GenomeMaker`, see `../data/README.md`) |
| Phobius (licensed) | `phobius_dir`, `run_phobius` | same |
| InterProScan data | `interproscan_data` (container) or `interproscan_dir` (local install) | data must match `interproscan_container` (5.73-104.0); with neither set, InterProScan is skipped with a warning. Container only: no conda package matches the 5.73 data, so `-profile conda` needs a local install in `interproscan_dir` |
| ProteinOrtho | `compare_synteny` (default `true`), `compare_proteinortho` | `compare_synteny: false` runs ProteinOrtho without `-synteny` as its own step (`PROTEINORTHO_NOSYNTENY`); `compare_proteinortho` hands in a ready-made table instead |
| Extra GenBank files | `compare_gbks` | optional, added to the compare |
| BinAC 2 | `project`, `apptainer_cache_dir` | `params/Porites_sp_v1_2_binac2.yaml` |

Genome names are letters, digits and `_` (they become the `funannotate -s` name and the label in the
compare). The new genome goes through the same gene-model preparation as the NCBI ones
(`bin/ncbi_ref_prep.pl`: sequence names of at most 16 characters, organelle sequences dropped,
gene → mRNA/tRNA/rRNA → exon/CDS only).

## BinAC 2

Nothing here is tied to a location. `params/Porites_sp_v1_2_binac2.yaml` uses paths relative to the
launch directory: copy `comparison/` to BinAC 2, put the new genome's `genome/*.scaffolds.fa` and
`annotation/*.gff3` in `inputs/` and the databases in `db/` (`funannotate_DB`, `eggnog_DB`, `phobius`,
`interproscan-5.73-104.0/data`, as in GenomeMaker's BinAC 2 setup), set `project` if an account is
needed. If the compute nodes cannot reach NCBI, download the references on a login node and list them
as `fasta` + `gff`.

## Results (`<outdir>`, default `results/<name>`)

```
00_inputs/            params.json, ncbi/ (downloaded FASTA + GFF3)
01_prepared/<sp>/     gene models (<sp>.models.gff3), seqmap, prep summary
02_proteins/<sp>/     proteins extracted from the models
03_functional_annotation/<sp>/{eggnog,interproscan,phobius,funannotate}/   funannotate/final_anno/annotate_results/<sp>.gbk
04_compare/           funannotate_compare/ (HTML report and tables)
pipeline_info/        trace, report, timeline
```

## Steps

| Process | Tool | Notes |
|---|---|---|
| `NCBI_DOWNLOAD` | ncbi-datasets-cli 18.39.0 (staphb container, conda) | only genomes given by accession |
| `PREP_REF_ANNOTATION` | AGAT 1.7.0, `bin/ncbi_ref_prep.pl` | normalise, shorten names, keep gene models |
| `DROP_UNTRANSLATABLE` | gffread 0.12.7 | mRNAs with an in-frame stop (`drop_untranslatable`) |
| `EXTRACT_PROTEINS` | gffread 0.12.7 | |
| `FUNANNOTATE_GENE_IDS` | awk | gene IDs without `_`, so funannotate keeps every annotation |
| `EGGNOG_MAPPER`, `INTERPROSCAN`, `PHOBIUS_RUN` | eggnog-mapper 2.1.15, InterProScan 5.73-104.0, Phobius 1.01 | |
| `FUNANNOTATE_ANNOTATE` | funannotate 1.8.17 | `--busco_db metazoa` |
| `FUNANNOTATE_COMPARE` | funannotate 1.8.17 | `--ml_method raxml` (image has no IQ-TREE), no dN/dS (no codeml), no outgroup |

## Synteny

`funannotate compare` runs ProteinOrtho with `-synteny` unless a table is handed in. Earlier runs: the
Porites-only synteny run finished (2026-10-13); with *A. millepora* added it did not (memory, 340 GB at
12 jobs). Contiguity of the genomes (from the prepared FASTAs and the release genome):

| Genome | Sequences | N50 | ≥ 10 Mb |
|---|---|---|---|
| *Porites* sp. (new) | 16,154 | 26.9 Mb | 14 |
| *P. lutea* | 106 | 36.0 Mb | 14 |
| *A. millepora* | 854 | 19.8 Mb | 14 |
| *P. lobata* | 1,098 | 2.2 Mb | 0 |
| *P. harrisoni* | 1,882 | 0.8 Mb | 0 |
| *P. evermanni* | 8,186 | 0.2 Mb | 0 |

The new assembly has 14 chromosome-scale scaffolds, so it is not what limits synteny; the three
fragmented references are. `compare_synteny: false` is the fallback.

## Earlier runs (`archive/`)

`archive/pre_nextflow_cleanup/` holds the previous layout as it was (shell runner `run_compare.sh`,
`genomes*.tsv`, old `main.nf`/`nextflow.config`, run logs, `refs/<species>/`, `compare_*/`, and the
unattended helpers of the 2026-10-08 to 10-13 runs in its `archive/`). Their Nextflow work
directories cannot be resumed from the new location.
