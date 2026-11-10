# Comparison: the new genome against the annotated NCBI *Porites*

Assembly and annotation comparison with other *Porites* (and *Acropora millepora* as an
outgroup). Each NCBI genome is reduced to bare gene models and re-annotated with the pipeline's
own stage 08 (Phobius, InterProScan, eggNOG-mapper, `funannotate annotate`), so that every genome
carries the same evidence as ours; then `funannotate compare` runs over all GenBank files.
The pipeline (GenomeMaker) is only imported, never copied or edited.

## What is needed

| Input | Where |
|---|---|
| GenomeMaker clone, linked as `../data/GenomeMaker` (stage-08 subworkflow, `config/` databases, containers/profiles) | see `../data/README.md` |
| Annotation of the new genome | `../data/results_v1.2/08_functional_annotation/funannotate/final_anno/annotate_results/Porites_sp.gbk` (`OWN_GBK`) |
| Genome lists | `genomes.tsv` (4 *Porites*), `genomes_outgroup.tsv` (+ *A. millepora*): `accession<TAB>species_name` |
| NCBI download | conda env `ncbi_datasets` (`DATASETS_ENV`); every other step runs in the pipeline's containers (`-profile local` by default) |
| Databases | `funannotate_DB`, `eggnog_DB`, `phobius` in `GenomeMaker/config/` (set in `nextflow.config`) |

## Run

```bash
cd comparison
./run_compare.sh annotate [species_name ...]      # fetch + re-annotate (one Nextflow run per genome in refs/<species>/)
./run_compare.sh compare                          # Porites only        -> compare_porites/
GENOMES_TSV=genomes_outgroup.tsv COMPARE_NAME=outgroup ./run_compare.sh compare   # with outgroup
./run_compare.sh all
```

`main.nf` (`--step annotate_ref | compare`) is what `run_compare.sh` launches; `modules/10_comparison.nf`
holds the comparison-specific processes and `bin/ncbi_ref_prep.pl` the NCBI-to-funannotate
preparation. Each Nextflow run starts from its own folder, so `.nextflow/` and `work/` stay apart;
re-running a step resumes it. Extra Nextflow arguments: `EXTRA_ARGS`.

## Results (not tracked; kept in place)

`refs/<species>/` (download, prepared models, annotation, work), `compare_porites/`,
`compare_outgroup/` (ProteinOrtho with synteny) and `compare_*_nosynteny/` (ProteinOrtho without
`-synteny`, via `--compare_proteinortho`). The synteny step was not viable for the outgroup set
(memory; that run was stopped), the Porites-only synteny run finished on 2026-10-13. Pfam coverage check of each annotation: 57% (*P. lobata*) to 73% (*A. millepora*) of mRNAs.

## Run history (`archive/`)

`archive/` has the unattended helpers of the 2026-10-08 to 10-13 runs (`orchestrate.sh`,
`hmmscan_gate.sh` which caps funannotate's dbCAN `hmmscan` jobs, `proteinortho_nosynteny.sh`,
`finish_nosynteny.sh`) and `STATUS.log`. They assume the old repository layout; kept as a record.
