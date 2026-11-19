// modules/06_astral.nf -- ASTRAL-IV species tree from the gene trees, rooted on the outgroup.
// Included by ../main.nf; results layout: the output block there.

process ASTRAL {
    conda 'bioconda::aster=1.25'
    container 'quay.io/biocontainers/aster:1.25--hc52dbad_1'
    tag "${taxa}"

    input:
    tuple val(taxa), path(gene_trees)

    output:
    tuple val(taxa), path("${taxa}_astral.tre"), emit: tree
    tuple val(taxa), path("${taxa}_astral.log"), emit: log

    script:
    """
    astral4 -i ${gene_trees} -o ${taxa}_astral.tre -u 2 --root ${params.outgroup} -t ${task.cpus} 2> ${taxa}_astral.log
    """

    stub:
    """
    echo '(${params.outgroup},(a,b));' > ${taxa}_astral.tre; touch ${taxa}_astral.log
    """
}
