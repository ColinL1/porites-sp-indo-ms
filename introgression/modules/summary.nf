// modules/summary.nf -- One row per test: summary/results.tsv.
// Included by ../main.nf; results layout: the output block there.

process SUMMARY {
    container params.tools_container

    input:
    path het,     stageAs: 'r/01_heterozygosity/*'
    path phase,   stageAs: 'r/02_phasing/*'
    path paint,   stageAs: 'r/04_ancestry_painting/*'
    path dstats,  stageAs: 'r/05_d_statistics/*'
    path twisst,  stageAs: 'r/06_twisst/*'
    path robust,  stageAs: 'r/07_twisst_robustness/*'
    path hyde,    stageAs: 'r/08_hyde/*'
    path pca,     stageAs: 'r/09_pca_admixture/*'
    path mito,    stageAs: 'r/10_mitogenome/*'
    path references

    output:
    path 'summary/results.tsv'

    script:
    """
    summary.py --outdir r --references ${references} --target_name ${params.name}
    mkdir summary; cp r/summary/results.tsv summary/
    """

    stub:
    """
    mkdir summary; printf 'test\\tresult\\tfolder\\n' > summary/results.tsv
    """
}
