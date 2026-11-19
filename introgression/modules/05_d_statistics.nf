// modules/05_d_statistics.nf -- Dsuite Dtrios, one run per outgroup.
// Included by ../main.nf; results layout: the output block there.

process DSUITE {
    tag "${og}"
    container params.tools_container

    input:
    tuple val(og), val(p1), val(p2), path(vcf), path(tbi)

    output:
    path "${og}"

    script:
    """
    mkdir ${og}
    bcftools view -s ${p1},${p2},${params.name},${og} ${vcf} -Oz -o quartet.vcf.gz
    printf "%s\\t%s\\n%s\\t%s\\n%s\\t%s\\n%s\\tOutgroup\\n" ${p1} ${p1} ${p2} ${p2} ${params.name} ${params.name} ${og} > ${og}/SETS.txt
    cd ${og} && Dsuite Dtrios ../quartet.vcf.gz SETS.txt > dtrios.log 2>&1
    """

    stub:
    """
    mkdir ${og}; touch ${og}/SETS.txt ${og}/dtrios.log
    printf 'P1\\tP2\\tP3\\tDstatistic\\tZ-score\\tp-value\\tf4-ratio\\tBBAA\\tABBA\\tBABA\\n${p1}\\t${p2}\\t${params.name}\\t0.01\\t1\\t0.3\\t0.001\\t3\\t1\\t1\\n' > ${og}/SETS_BBAA.txt
    """
}
