# genome-assembly-annotation

The assembly and annotation pipeline is a separate repository, **GenomeMaker** (`github.com/ColinL1/GenomeMaker`).
This folder only holds a link, `results_v1.2`, to that pipeline's v1.2 results (untracked), so that the analyses here
read the final genome, annotation and mitogenome from one place:

```bash
ln -s <path>/GenomeMaker/results_v1.2 genome-assembly-annotation/results_v1.2
```
