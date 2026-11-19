// modules/03_reference_alignment.nf -- Every reference assembly projected onto the target (minimap2, paftools).
// Included by ../main.nf; results layout: the output block there.

process ALIGN_REFERENCE {
    tag "${name}"
    container params.tools_container

    input:
    tuple val(name), path(fna), val(preset)
    tuple path(target), path(fai)

    output:
    path "${name}.{callable.bed,vcf.gz,vcf.gz.tbi}", emit: files
    tuple val(name), path("${name}.{callable.bed,vcf.gz,vcf.gz.tbi}"), emit: by_name
    path "${name}.stats.tsv", emit: stats

    script:
    """
    mkdir work
    align_reference.sh ${name} ${fna} ${preset} ${target} . work ${task.cpus}
    mv work/${name}.stats.tsv .
    """

    stub:
    """
    touch ${name}.callable.bed ${name}.vcf.gz ${name}.vcf.gz.tbi; printf '${name}\\t${preset}\\t1\\t1\\n' > ${name}.stats.tsv
    """
}
