// ============================================================
// ncbi_download.nf
// Genome FASTA + GFF3 of one NCBI accession (only for genomes listed by accession).
// Genomes given as fasta + gff paths in the params file skip this step, which is also
// the way around compute nodes without internet access.
// ============================================================

nextflow.enable.dsl = 2

process NCBI_DOWNLOAD {
    tag   "${sp}"

    conda 'conda-forge::ncbi-datasets-cli=18.39.0'
    container 'docker.io/staphb/ncbi-datasets:18.39.0'

    input:
    tuple val(sp), val(accession)

    output:
    tuple val(sp), path("${accession}.genomic.fna"), path("${accession}.genomic.gff"), emit: genome
    path "versions.yml", emit: versions

    script:
    """
    datasets download genome accession ${accession} \\
      --include genome,gff3 --no-progressbar --filename ncbi_dataset.zip
    unzip -q -o ncbi_dataset.zip

    if [ ! -s ncbi_dataset/data/${accession}/genomic.gff ]; then
        echo "ERROR: ${accession} has no GFF3 annotation on NCBI" >&2
        exit 1
    fi
    mv ncbi_dataset/data/${accession}/*_genomic.fna ${accession}.genomic.fna
    mv ncbi_dataset/data/${accession}/genomic.gff   ${accession}.genomic.gff
    rm -rf ncbi_dataset ncbi_dataset.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets-cli: "\$(datasets --version 2>&1 | grep -oim1 '[0-9][0-9.]*' || echo unknown)"
    END_VERSIONS
    """

    stub:
    """
    touch ${accession}.genomic.fna ${accession}.genomic.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets-cli: "stub"
    END_VERSIONS
    """
}
