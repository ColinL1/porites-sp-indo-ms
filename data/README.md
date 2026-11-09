# data/ (not tracked)

Links to inputs that live outside this repo. Everything here is ignored by git; recreate the links on a new machine.

| Link | Target | Used by |
|---|---|---|
| `results_v1.2` | `../genome-assembly-annotation/results_v1.2` -> `GenomeMaker/results_v1.2` | all |
| `GenomeMaker` | clone of the pipeline repo | `comparison/` |
| `ont_reads.fastq.gz` | ONT reads of the sequenced individual | `introgression/` |
| `introgression_refs/` | NCBI FASTAs of the reference panel (`introgression/references.tsv`) | `introgression/` |
| `busco_refs/` | `busco_<name>/` runs of the reference genomes | `species_tree/` |
| `busco_downloads/` | offline BUSCO lineage datasets (anthozoa_odb12) | `species_tree/` |

Files of `results_v1.2` that the analyses read: `final/genome/Porites_sp.fa` (+ `scaffold_name_map.tsv`),
`final/mitogenome/mitogenome.fasta`, `08_functional_annotation/.../Porites_sp.gbk`.

```bash
ln -s ../genome-assembly-annotation/results_v1.2 data/results_v1.2
ln -s <path>/GenomeMaker data/GenomeMaker
```
