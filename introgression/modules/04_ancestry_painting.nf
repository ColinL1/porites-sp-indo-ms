// modules/04_ancestry_painting.nf -- Diagnostic sites of the two parents and the painting of the target haplotypes.
// Included by ../main.nf; results layout: the output block there.

process ANCESTRY_PAINTING {
    tag "${params.name}"
    container params.tools_container

    input:
    path references
    path aln, stageAs: 'aln/*'
    tuple path(target), path(fai)
    tuple path(bam), path(bai)

    output:
    path 'ancestry_summary.txt', emit: summary
    path '*.{tsv,tsv.gz,png}', emit: files

    script:
    """
    diagnostic_sites.py --references ${references} --aln_dir aln --target ${target} --out diagnostic_sites.tsv.gz
    ancestry_painting.py --diag diagnostic_sites.tsv.gz --bam ${bam} --fai ${fai} --references ${references} \\
      --outdir . --threads ${task.cpus}
    """

    stub:
    """
    touch diagnostic_sites.tsv.gz ancestry_painting.png
    printf 'diagnostic sites: 1\\nsites with >= 1 haplotype-resolved call: 1\\n  HOM_P1 1 (50.0%%)\\n  HOM_P2 1 (50.0%%)\\nallele matches genome-wide: P1=1  P2=1  P1_fraction=0.5\\nphase blocks with >= 10 resolved sites: 1\\n' > ancestry_summary.txt
    """
}
