// modules/01_heterozygosity.nf -- Reads mapped to the target, SNVs called (Clair3) and filtered, read depth.
// Included by ../main.nf; results layout: the output block there.

process MAP_READS {
    tag "${params.name}"
    container params.tools_container

    input:
    tuple path(target), path(fai)
    path reads

    output:
    tuple path('reads.bam'), path('reads.bam.bai'), emit: bam
    path 'reads.meandepth.txt', emit: meandepth
    path 'reads.flagstat.txt', emit: flagstat

    script:
    """
    map_reads.sh ${target} ${reads} ${params.name} reads.bam ${task.cpus}
    """

    stub:
    """
    touch reads.bam reads.bam.bai reads.flagstat.txt; printf 'mean_depth_all_positions\\t18.7\\ttotal_pos\\t4\\n' > reads.meandepth.txt
    """
}

process CLAIR3 {
    tag "${params.name}"
    container 'quay.io/biocontainers/clair3:2.0.3--py311hbc58adc_0'

    input:
    tuple path(target), path(fai)
    tuple path(bam), path(bai)
    path bed

    output:
    tuple path('clair3/pileup.vcf.gz'), path('clair3/pileup.vcf.gz.tbi')

    script:
    """
    run_clair3.sh --bam_fn=${bam} --ref_fn=${target} --threads=${task.cpus} --platform=ont \\
      --model_path=${params.clair3_model} --output=clair3 --bed_fn=${bed} --min_coverage=4 \\
      --sample_name=${params.name}
    """

    stub:
    """
    mkdir clair3; touch clair3/pileup.vcf.gz clair3/pileup.vcf.gz.tbi
    """
}

process SNV_FILTER {
    tag "${params.name}"
    container params.tools_container

    input:
    tuple path(pileup), path(tbi)
    path bed
    path meandepth

    output:
    tuple path('snvs.vcf.gz'), path('snvs.vcf.gz.tbi'), emit: snvs
    path 'heterozygosity.tsv', emit: het

    script:
    """
    bcftools view -f PASS -v snps -m2 -M2 -e 'QUAL<10' ${pileup} -Oz -o snvs.vcf.gz
    tabix -f -p vcf snvs.vcf.gz
    het=\$(bcftools view -H -g het snvs.vcf.gz | wc -l)
    hom=\$(bcftools view -H -g hom snvs.vcf.gz | wc -l)
    bp=\$(awk '{s+=\$3-\$2} END{print s}' ${bed})
    printf "het_snvs\\t%s\\nhom_alt_snvs\\t%s\\nscaffold_bp\\t%s\\nhet_per_bp\\t%s\\nmean_depth\\t%s\\n" \\
      \$het \$hom \$bp \$(awk -v h=\$het -v b=\$bp 'BEGIN{printf "%.6f", h/b}') \$(cut -f2 ${meandepth}) > heterozygosity.tsv
    """

    stub:
    """
    touch snvs.vcf.gz snvs.vcf.gz.tbi
    printf 'het_snvs\\t1\\nhom_alt_snvs\\t1\\nscaffold_bp\\t4\\nhet_per_bp\\t0.25\\nmean_depth\\t18.7\\n' > heterozygosity.tsv
    """
}

process MOSDEPTH {
    tag "${params.name}"
    container 'quay.io/biocontainers/mosdepth:0.3.14--h87be163_2'

    input:
    tuple path(bam), path(bai)

    output:
    path 'depth.regions.bed.gz', emit: regions
    path 'depth.*', emit: files

    script:
    """
    mosdepth -t ${task.cpus} -n -b 5000 depth ${bam}
    """

    stub:
    """
    touch depth.regions.bed.gz depth.mosdepth.summary.txt
    """
}
