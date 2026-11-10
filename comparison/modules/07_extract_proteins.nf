// ============================================================
// 07_extract_proteins.nf
// Protein FASTA from a gene-model GFF3 (input of the functional annotation).
// Copied from GenomeMaker/modules/07_structural_annotation.nf (EXTRACT_PROTEINS), with the
// genome name as a process input instead of the global params.species_name, so that
// several genomes can be annotated in one run.
// ============================================================

nextflow.enable.dsl = 2

process EXTRACT_PROTEINS {
    tag   "${sp}"

    conda 'bioconda::gffread=0.12.7'
    container 'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4'

    input:
    tuple val(sp), path(gff3), path(assembly)

    output:
    tuple val(sp), path("${sp}.proteins.faa"), emit: proteins_faa
    path "versions.yml", emit: versions

    script:
    """
    gffread ${gff3} \\
      -g ${assembly} \\
      -y ${sp}.proteins.faa

    if [ ! -s "${sp}.proteins.faa" ]; then
        echo "ERROR: ${sp}.proteins.faa is empty" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "\$(gffread --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${sp}.proteins.faa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "stub"
    END_VERSIONS
    """
}
