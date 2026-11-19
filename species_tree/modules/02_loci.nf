// modules/02_loci.nf -- Loci complete and single-copy in every genome of a taxon set (protein sequences).
// Included by ../main.nf; results layout: the output block there.

process SHARED_LOCI {
    conda 'conda-forge::python=3.12.13'
    container 'community.wave.seqera.io/library/python:3.12.13--bf5c18fa662d47d8'
    tag "${taxa}"

    input:
    tuple val(taxa), val(names), path(busco_dirs, stageAs: 'b??/*')

    output:
    tuple val(taxa), path('loci'), emit: loci
    tuple val(taxa), path('shared_genes.txt'), emit: genes

    script:
    def pairs = [names, busco_dirs].transpose().collect { n, d -> "${n}=${d}" }.join(' ')
    """
    shared_loci.py ${params.lineage} loci ${pairs}
    mv loci/shared_genes.txt .
    """

    stub:
    """
    mkdir loci; for i in 1 2 3; do printf '>a\\nM\\n' > loci/l\$i.fasta; done; touch shared_genes.txt
    """
}
