// modules/joint_vcf.nf -- Joint VCFs of the target and the references (quartet: steps 05-07; panel: 08-09).
// Included by ../main.nf; results layout: the output block there.

process JOINT_VCF {
    tag "${set}"
    container params.tools_container

    input:
    val set
    path references
    path aln, stageAs: 'aln/*'
    tuple path(pileup), path(pileup_tbi)
    tuple path(target), path(fai)

    output:
    tuple val(set), path("joint_${set}.vcf.gz"), path("joint_${set}.vcf.gz.tbi")

    script:
    """
    joint_vcf.py --set ${set} --references ${references} --target_name ${params.name} --aln_dir aln \\
      --clair3_vcf ${pileup} --target ${target} --out joint_${set}.vcf
    bgzip -f joint_${set}.vcf; tabix -f -p vcf joint_${set}.vcf.gz
    """

    stub:
    """
    touch joint_${set}.vcf.gz joint_${set}.vcf.gz.tbi
    """
}
