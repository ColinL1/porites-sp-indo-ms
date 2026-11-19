// modules/05_gene_trees.nf -- One ML tree per locus (IQ-TREE -S).
// Included by ../main.nf; results layout: the output block there.

process GENE_TREES {
    conda 'bioconda::iqtree=3.1.3'
    container 'quay.io/biocontainers/iqtree:3.1.3--h8471819_0'
    tag "${taxa}"

    input:
    tuple val(taxa), path(trimmed)

    output:
    tuple val(taxa), path("${taxa}_loci.treefile"), emit: trees
    tuple val(taxa), path("${taxa}_loci.{iqtree,log,best_scheme.nex}"), emit: files

    script:
    """
    iqtree -S ${trimmed} -m MFP -T ${task.cpus} --prefix ${taxa}_loci
    """

    stub:
    """
    echo '(${params.outgroup},(a,b));' > ${taxa}_loci.treefile; touch ${taxa}_loci.iqtree ${taxa}_loci.log ${taxa}_loci.best_scheme.nex
    """
}
