// modules/07_twisst_robustness.nf -- Twisst over the grid of window sizes, tree methods and target representations.
// Included by ../main.nf; results layout: the output block there.

process ROBUSTNESS_TREES {
    container params.tools_container

    input:
    path geno
    tuple path(phased), path(phased_tbi)
    path blocks
    path references
    val og

    output:
    path '*_trees.nwk', emit: trees
    path '*_{tree_window.txt,windows.tsv,groups.tsv}', emit: meta

    script:
    """
    twisst_robustness.py --geno ${geno} --phased_vcf ${phased} --block_list ${blocks} --references ${references} \\
      --target_name ${params.name} --outgroup ${og} --outdir . --threads ${task.cpus}
    """

    stub:
    """
    for m in rand2 phased2 cons1; do touch \${m}_groups.tsv; for w in 50 200; do for t in NJ ML; do
      touch \${m}_W\${w}_\${t}_trees.nwk \${m}_W\${w}_\${t}_tree_window.txt \${m}_W\${w}_\${t}_windows.tsv; done; done; done
    """
}

process TWISST_CELL {
    tag "${cell}"
    container params.tools_container

    input:
    tuple val(cell), path(trees), val(p1), val(p2), val(og)
    path meta

    output:
    path "${cell}_{weights.csv.gz,windows.tsv}", includeInputs: true, emit: weights
    tuple val(cell), path("${cell}_{weights.csv.gz,weights_per_tree.csv.gz,topos.nwk}"), emit: by_cell

    script:
    def mode = cell.tokenize('_')[0]
    """
    python \$TWISST -t ${cell}_trees.nwk -w ${cell}_weights_per_tree.csv.gz --outputTopos ${cell}_topos.nwk \\
      -g ${p1} -g ${p2} -g ${params.name} -g ${og} --method complete --groupsFile ${mode}_groups.tsv --silent
    twisst_average.py ${cell}_weights_per_tree.csv.gz ${cell}_tree_window.txt ${cell}_weights.csv.gz
    """

    stub:
    """
    touch ${cell}_weights.csv.gz ${cell}_weights_per_tree.csv.gz ${cell}_topos.nwk
    """
}

process ROBUSTNESS_SUMMARY {
    container params.tools_container

    input:
    path weights, stageAs: 'grid/*'
    path depth
    path aln, stageAs: 'aln/*'
    path references
    path agp

    output:
    path 'robustness_summary.tsv', emit: summary
    path 'context_summary.txt', emit: context
    path 'context.tsv', emit: context_windows

    script:
    def agp_opt = agp.name != 'NO_AGP' ? "--agp ${agp}" : ''
    """
    robustness_summary.py --dir grid --references ${references} --target_name ${params.name}
    context.py --dir grid --depth ${depth} --aln_dir aln --references ${references} --target_name ${params.name} ${agp_opt}
    mv grid/robustness_summary.tsv grid/context.tsv grid/context_summary.txt .
    """

    stub:
    """
    printf 'representation\\twindow_snps\\tmethod\\tn_windows\\tT1\\tT2\\tT3\\tplurality\\n' > robustness_summary.tsv
    printf 'context\\tr_T1\\tr_T2\\tr_T3\\ndepth\\t+0.010\\t+0.010\\t+0.010\\n' > context_summary.txt; touch context.tsv
    """
}
