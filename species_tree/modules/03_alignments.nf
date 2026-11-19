// modules/03_alignments.nf -- Per-locus alignment (MAFFT) and trimming (trimAl).
// Included by ../main.nf; results layout: the output block there.

process ALIGN_TRIM {
    conda 'conda-forge::mafft=7.526 bioconda::trimal=1.5.1'
    container 'community.wave.seqera.io/library/mafft_trimal:e2650fe4ee2c93cc'
    tag "${taxa}"

    input:
    tuple val(taxa), path(loci)

    output:
    tuple val(taxa), path('trimmed')

    script:
    """
    mkdir aligned trimmed
    ls ${loci}/*.fasta | xargs -P ${task.cpus} -I{} sh -c '
      g=\$(basename {} .fasta)
      mafft --auto --quiet {} > aligned/\$g.aln
      trimal -in aligned/\$g.aln -out trimmed/\$g.trim.aln -automated1'
    """

    stub:
    """
    mkdir trimmed; touch trimmed/l1.trim.aln
    """
}
