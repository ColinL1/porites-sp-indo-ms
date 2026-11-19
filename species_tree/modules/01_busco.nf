// modules/01_busco.nf -- BUSCO (genome mode, offline) on the new assembly.
// Included by ../main.nf; results layout: the output block there.

process BUSCO_NEW {
    conda 'bioconda::busco=6.1.0'
    container 'quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1'
    tag "${name}"

    input:
    tuple val(name), path(fasta)

    output:
    tuple val(name), path("busco_${name}")

    script:
    """
    busco -i ${fasta} -m genome -l ${params.lineage} -c ${task.cpus} -o busco_${name} \\
      --offline --download_path ${params.busco_download_path}
    """

    stub:
    """
    mkdir -p busco_${name}/run_${params.lineage}/busco_sequences/single_copy_busco_sequences
    """
}
