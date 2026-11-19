// modules/04_concat.nf -- Concatenated supermatrix (AMAS) and partitioned ML tree (IQ-TREE).
// Included by ../main.nf; results layout: the output block there.

process CONCAT_TREE {
    conda 'bioconda::amas=1.0 bioconda::iqtree=3.1.3'
    container 'community.wave.seqera.io/library/amas_iqtree:58bbd76a2b8ed70b'
    tag "${taxa}"

    input:
    tuple val(taxa), path(trimmed)

    output:
    tuple val(taxa), path("${taxa}_concat.treefile"), emit: tree
    tuple val(taxa), path("${taxa}_concat.*"), path('supermatrix.fasta'), path('partitions.txt'), emit: files

    script:
    """
    AMAS.py concat -i ${trimmed}/*.trim.aln -f fasta -d aa -t supermatrix.fasta -p partitions.txt
    sed -i 's/^/AUTO, /' partitions.txt
    iqtree -s supermatrix.fasta -p partitions.txt -m MFP -bb 1000 -alrt 1000 -T ${task.cpus} -pre ${taxa}_concat
    """

    stub:
    """
    touch supermatrix.fasta partitions.txt ${taxa}_concat.iqtree
    echo '(${params.outgroup},(a,b));' > ${taxa}_concat.treefile
    """
}
