# data/ (not tracked)

Links to inputs that live outside this repo. Everything here is ignored by git; recreate the links on a new machine.

| Link | Target | Used by |
| --- | --- | --- |
| `results_v1.2` | `../genome-assembly-annotation/results_v1.2` -> `GenomeMaker/results_v1.2` | all |
| `GenomeMaker` | clone of the pipeline repo | `comparison/` (database paths in `params/Porites_sp_v1_2.yaml` only) |
| `interproscan_data` | InterProScan 5.73-104.0 `data/` folder | `comparison/` (`params/Porites_sp_v1_2.yaml`) |
| `ont_reads.fastq.gz` | ONT reads of the sequenced individual | `introgression/` |
| `introgression_refs/` | NCBI FASTAs of the reference panel (`introgression/references.tsv`) | `introgression/` |
| `busco_refs/` | `busco_<name>/` runs of the reference genomes | `species_tree/` |
| `busco_downloads/` | offline BUSCO lineage datasets (anthozoa_odb12) | `species_tree/` |

Files of `results_v1.2` that the analyses read (all in `release_curation/Porites_sp_v1.2/`):
`genome/Porites_sp._Indonesia_ID_GA_Por_15.scaffolds.fa` (all three), `genome/...scaffolds.agp` and
`mitogenome/mitogenome.fasta` (introgression), `annotation/Porites_sp._Indonesia_ID_GA_Por_15.gff3`
(comparison).

```bash
ln -s ../genome-assembly-annotation/results_v1.2 data/results_v1.2
ln -s <path>/GenomeMaker data/GenomeMaker
```
