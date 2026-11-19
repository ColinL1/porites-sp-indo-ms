// modules/02_phasing.nf -- Phasing per scaffold (WhatsHap), phase stats and the haplotagged BAM.
// Included by ../main.nf; results layout: the output block there.

process SPLIT_SNVS {
    container params.tools_container

    input:
    tuple path(snvs), path(tbi)
    path bed

    output:
    path 'per/*.in.vcf.gz'

    script:
    """
    mkdir per
    cut -f1 ${bed} | xargs -P ${task.cpus} -I{} sh -c 'bcftools view -r {} ${snvs} -Oz -o per/{}.in.vcf.gz'
    """

    stub:
    """
    mkdir per; cut -f1 ${bed} | while read c; do touch per/\$c.in.vcf.gz; done
    """
}

process PHASE_SCAFFOLD {
    tag "${chrom}"
    container 'quay.io/biocontainers/whatshap:2.8--py311he264feb_0'

    input:
    tuple val(chrom), path(vcf)
    tuple path(target), path(fai)
    tuple path(bam), path(bai)

    output:
    path "${chrom}.phased.vcf.gz"

    script:
    """
    if [ "\$(zcat ${vcf} | grep -vc '^#')" -lt 2 ]; then
      cp ${vcf} ${chrom}.phased.vcf.gz
    else
      whatshap phase --ignore-read-groups --chromosome ${chrom} -o ${chrom}.phased.vcf.gz \\
        --reference ${target} ${vcf} ${bam}
    fi
    """

    stub:
    """
    touch ${chrom}.phased.vcf.gz
    """
}

process PHASE_CONCAT {
    container params.tools_container

    input:
    path vcfs
    path bed

    output:
    tuple path('phased.vcf.gz'), path('phased.vcf.gz.tbi')

    script:
    """
    for f in ${vcfs}; do tabix -f -p vcf \$f; done
    bcftools concat -Oz -o phased.vcf.gz \$(cut -f1 ${bed} | sed 's/\$/.phased.vcf.gz/')
    tabix -f -p vcf phased.vcf.gz
    """

    stub:
    """
    touch phased.vcf.gz phased.vcf.gz.tbi
    """
}

process PHASE_STATS_HAPLOTAG {
    tag "${params.name}"
    container 'quay.io/biocontainers/whatshap:2.8--py311he264feb_0'

    input:
    tuple path(phased), path(tbi)
    tuple path(target), path(fai)
    tuple path(bam), path(bai)

    output:
    tuple path('haplotagged.bam'), path('haplotagged.bam.bai'), emit: bam
    path 'phase_block_list.tsv', emit: blocks
    path 'phase_block_stats.tsv', emit: stats
    path 'phase_stats.txt', emit: summary

    script:
    """
    whatshap stats ${phased} --tsv=phase_block_stats.tsv --block-list=phase_block_list.tsv \\
      --chr-lengths <(cut -f1,2 ${fai}) > phase_stats.txt 2>&1
    whatshap haplotag --ignore-read-groups --skip-missing-contigs -o haplotagged.bam --reference ${target} \\
      ${phased} ${bam}
    python -c "import pysam; pysam.index('-@', '${task.cpus}', 'haplotagged.bam')"
    """

    stub:
    """
    touch haplotagged.bam haplotagged.bam.bai phase_block_list.tsv phase_stats.txt
    printf 'sample\\tchromosome\\tblocks\\tphased_fraction\\tblock_n50\\nT\\tALL\\t1\\t0.99\\t21000\\n' > phase_block_stats.tsv
    """
}
