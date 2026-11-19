// modules/08_hyde.nf -- HyDe hybridization test, one run per non-parent species as outgroup.
// Included by ../main.nf; results layout: the output block there.

process HYDE {
    container params.tools_container

    input:
    tuple val(set), path(vcf), path(tbi)
    path references

    output:
    path 'hyde_summary.tsv', emit: summary
    path '{seqs.txt,map.txt,hyde_*-out*.txt}', emit: raw

    script:
    """
    hyde.py --vcf ${vcf} --references ${references} --target_name ${params.name} --outdir .
    """

    stub:
    """
    printf 'outgroup\\tP1\\tHybrid\\tP2\\tZscore\\tPvalue\\tGamma\\n' > hyde_summary.tsv; touch seqs.txt map.txt
    """
}
