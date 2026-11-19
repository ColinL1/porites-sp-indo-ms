// modules/07_plots.nf -- Analysis plots of a taxon set: quartet support, tanglegram, normalised RF.
// Included by ../main.nf; results layout: the output block there.

process PLOT {
    conda 'conda-forge::r-base=4.5.2 conda-forge::r-ape=5.8_1 conda-forge::r-phytools=2.5_2 conda-forge::r-phangorn=2.12.1'
    container 'community.wave.seqera.io/library/r-base_r-ape_r-phytools_r-phangorn:be3443325fe730fa'
    tag "${taxa}"

    input:
    tuple val(taxa), path(astral), path(concat)

    output:
    tuple val(taxa), path("${taxa}_*")

    script:
    """
    plot_trees.R ${astral} ${concat} ${params.outgroup} ${taxa}
    """

    stub:
    """
    touch ${taxa}_astral_quartets.pdf ${taxa}_tanglegram.pdf ${taxa}_rf.tsv
    """
}
