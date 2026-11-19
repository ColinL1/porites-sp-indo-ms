// modules/10_mitogenome.nf -- Mitogenome of the target (given, or assembled from the reads) and the mito tree.
// Included by ../main.nf; results layout: the output block there.

process MITO_READS {
    container params.tools_container

    input:
    path bait
    path reads

    output:
    path 'mito_reads.fastq.gz'

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${bait} ${reads} 2>/dev/null | samtools view -b -F 0x904 - | samtools sort -o mito.bam
    samtools fastq mito.bam | gzip > mito_reads.fastq.gz
    """

    stub:
    """
    touch mito_reads.fastq.gz
    """
}

process MITO_SUBSAMPLE {
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    input:
    path fq

    output:
    path 'mito_sub.fastq.gz'

    script:
    """
    seqkit sample -s11 -p 0.15 ${fq} -o mito_sub.fastq.gz
    """

    stub:
    """
    touch mito_sub.fastq.gz
    """
}

process MITO_FLYE {
    container 'quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1'

    input:
    path fq

    output:
    path 'assembly.fasta', emit: fasta
    path 'flye_assembly_info.txt', emit: info

    script:
    """
    flye --nano-hq ${fq} --out-dir flye --threads ${task.cpus} --min-overlap 1000
    cp flye/assembly.fasta assembly.fasta; cp flye/assembly_info.txt flye_assembly_info.txt
    """

    stub:
    """
    touch assembly.fasta flye_assembly_info.txt
    """
}

process MITO_SEQUENCES {
    container params.tools_container

    input:
    path target_mito
    val source
    path bait
    path references

    output:
    path 'mito_seqs.fasta', emit: fasta
    path 'mito_sources.tsv', emit: sources

    script:
    """
    mito_sequences.py --references ${references} --target_name ${params.name} --target_mito ${target_mito} \\
      --target_source "${source}" --bait ${bait} --outdir . --threads ${task.cpus}
    """

    stub:
    """
    printf '>a\\nA\\n' > mito_seqs.fasta; printf 'sample\\tsource\\tlength\\n${params.name}\\tgiven\\t18000\\n' > mito_sources.tsv
    """
}

process MITO_ALIGN {
    container 'quay.io/biocontainers/mafft:7.525--h031d066_1'

    input:
    path seqs

    output:
    path 'mito_aln.fasta', emit: aln
    path 'mafft.log', emit: log

    script:
    """
    mafft --auto --adjustdirection --thread ${task.cpus} ${seqs} > mito_aln.fasta 2> mafft.log
    sed -i 's/^>_R_/>/' mito_aln.fasta
    """

    stub:
    """
    cp ${seqs} mito_aln.fasta; touch mafft.log
    """
}

process MITO_TREE {
    container 'quay.io/biocontainers/iqtree:3.1.3--h8471819_0'

    input:
    path aln

    output:
    path 'mito_tree.treefile', emit: tree
    path 'mito_tree.*', emit: files

    script:
    """
    iqtree -s ${aln} -m MFP -B 1000 -seed 1 -T ${task.cpus} --prefix mito_tree -redo > iqtree.log 2>&1
    """

    stub:
    """
    echo '(a,b,c);' > mito_tree.treefile; touch mito_tree.iqtree
    """
}

process MITO_SUMMARY {
    container params.tools_container

    input:
    path tree

    output:
    path 'mito_tree_summary.txt'

    script:
    """
    mito_tree_summary.py --tree ${tree} --root ${params.mito_root} --target_name ${params.name} > mito_tree_summary.txt
    """

    stub:
    """
    echo '${params.name} sister group: x (support 50)' > mito_tree_summary.txt
    """
}
