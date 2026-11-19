// modules/06_twisst.nf -- Genotypes, Twisst topology weights per outgroup, and their plots.
// Included by ../main.nf; results layout: the output block there.

process GENO {
    container params.tools_container

    input:
    tuple val(set), path(vcf), path(tbi)

    output:
    path 'joint_quartet.geno.gz'

    script:
    """
    python \$PARSEVCF -i ${vcf} -o joint_quartet.geno.gz
    """

    stub:
    """
    touch joint_quartet.geno.gz
    """
}

process TWISST {
    tag "${og}"
    container params.tools_container

    input:
    tuple val(og), val(p1), val(p2), path(geno), path(references)

    output:
    path "${og}_{weights.csv.gz,windows.tsv}", emit: weights
    tuple val(og), path("${og}_*"), emit: by_outgroup

    script:
    """
    twisst_trees.py --geno ${geno} --references ${references} --target_name ${params.name} --outgroup ${og} \\
      --window 50 --prefix ${og}
    python \$TWISST -t ${og}_trees.nwk -w ${og}_weights_per_tree.csv.gz --outputTopos ${og}_topos.nwk \\
      -g ${p1} -g ${p2} -g ${params.name} -g ${og} --method complete --groupsFile ${og}_groups.tsv --silent
    twisst_average.py ${og}_weights_per_tree.csv.gz ${og}_tree_window.txt ${og}_weights.csv.gz
    """

    stub:
    """
    touch ${og}_weights.csv.gz ${og}_windows.tsv ${og}_trees.nwk
    """
}

process TWISST_PLOT {
    container params.tools_container

    input:
    path weights, stageAs: 'tw/*'
    path references
    tuple path(target), path(fai)

    output:
    path 'twisst_summary.tsv', emit: summary
    path '*.png', emit: plots

    script:
    """
    twisst_plot.py --dir tw --references ${references} --target_name ${params.name} --fai ${fai}
    mv tw/twisst_summary.tsv tw/*.png .
    """

    stub:
    """
    printf 'outgroup\\tn_windows\\tT1_x\\tT2_x\\tT3_x\\tplurality\\n' > twisst_summary.tsv; touch twisst_painting.png
    """
}
