// BUSCO single-copy phylogeny of the new Porites assembly, two taxon sets (small, large) x
// two approaches (concatenated ML, ASTRAL-IV) = 4 trees, plus the quartet and tanglegram plots.
//
//   nextflow run main.nf --new_fasta <final assembly .fa>            (BUSCO is run on it)
//   nextflow run main.nf --new_busco_dir <existing BUSCO genome-mode output dir>
//
// Steps as in the earlier runs (earlier single_ref and astral_expanded runs, see archive/), except
// ASTRAL-IV is run unconstrained (the earlier runs passed the concatenated tree with -c, which in
// ASTRAL-IV is a topology constraint).

nextflow.enable.dsl = 2

process BUSCO_NEW {
    tag "${name}"
    publishDir "${params.outdir}/busco", mode: 'copy'

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

process SHARED_LOCI {
    tag "${taxa}"
    publishDir path: { "${params.outdir}/${taxa}" }, mode: 'copy', pattern: 'shared_genes.txt'

    input:
    tuple val(taxa), val(names), path(busco_dirs, stageAs: 'b??/*')

    output:
    tuple val(taxa), path('loci'), emit: loci
    path 'shared_genes.txt'

    script:
    def pairs = [names, busco_dirs].transpose().collect { n, d -> "${n}=${d}" }.join(' ')
    """
    shared_loci.py ${params.lineage} loci ${pairs}
    mv loci/shared_genes.txt .
    """

    stub:
    """
    mkdir loci; for i in 1 2 3; do printf '>a\\nM\\n' > loci/l\$i.fasta; done; touch shared_genes.txt
    """
}

process ALIGN_TRIM {
    tag "${taxa}"
    publishDir path: { "${params.outdir}/${taxa}" }, mode: 'copy'

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

process CONCAT_TREE {
    tag "${taxa}"
    publishDir path: { "${params.outdir}/${taxa}/concat" }, mode: 'copy'

    input:
    tuple val(taxa), path(trimmed)

    output:
    tuple val(taxa), path("${taxa}_concat.treefile"), emit: tree
    path "${taxa}_concat.*"
    path 'supermatrix.fasta'
    path 'partitions.txt'

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

process GENE_TREES {
    tag "${taxa}"
    publishDir path: { "${params.outdir}/${taxa}/gene_trees" }, mode: 'copy'

    input:
    tuple val(taxa), path(trimmed)

    output:
    tuple val(taxa), path("${taxa}_loci.treefile"), emit: trees
    path "${taxa}_loci.{iqtree,log,best_scheme.nex}"

    script:
    """
    iqtree -S ${trimmed} -m MFP -T ${task.cpus} --prefix ${taxa}_loci
    """

    stub:
    """
    echo '(${params.outgroup},(a,b));' > ${taxa}_loci.treefile; touch ${taxa}_loci.iqtree ${taxa}_loci.log ${taxa}_loci.best_scheme.nex
    """
}

process ASTRAL {
    tag "${taxa}"
    publishDir path: { "${params.outdir}/${taxa}/astral" }, mode: 'copy'

    input:
    tuple val(taxa), path(gene_trees)

    output:
    tuple val(taxa), path("${taxa}_astral.tre"), emit: tree
    path "${taxa}_astral.log"

    script:
    """
    astral4 -i ${gene_trees} -o ${taxa}_astral.tre -u 2 --root ${params.outgroup} -t ${task.cpus} 2> ${taxa}_astral.log
    """

    stub:
    """
    echo '(${params.outgroup},(a,b));' > ${taxa}_astral.tre; touch ${taxa}_astral.log
    """
}

process PLOT {
    tag "${taxa}"
    publishDir path: { "${params.outdir}/${taxa}/plots" }, mode: 'copy'

    input:
    tuple val(taxa), path(astral), path(concat)

    output:
    path "${taxa}_*"

    script:
    """
    plot_trees.R ${astral} ${concat} ${params.outgroup} ${taxa}
    """

    stub:
    """
    touch ${taxa}_astral_quartets.pdf ${taxa}_tanglegram.pdf ${taxa}_rf.tsv
    """
}

workflow {
    if (!(params.new_fasta || params.new_busco_dir)) {
        error "Give the new assembly with --new_fasta, or its BUSCO output with --new_busco_dir"
    }
    def new_busco = params.new_busco_dir
        ? channel.value([params.new_name, file(params.new_busco_dir, checkIfExists: true)])
        : BUSCO_NEW(channel.of([params.new_name, file(params.new_fasta, checkIfExists: true)])).first()

    def refs = channel.fromPath(params.references, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')

    def sets = channel.of('small', 'large')
        .combine(refs)
        .filter { taxa, r -> r[taxa] == 'yes' }
        .map { taxa, r -> [taxa, r.name, file(r.busco_dir.replace("@BUSCO_ROOT@", params.busco_root.toString()), checkIfExists: true)] }
        .mix(channel.of('small', 'large').combine(new_busco))
        .groupTuple()

    SHARED_LOCI(sets)
    ALIGN_TRIM(SHARED_LOCI.out.loci)
    CONCAT_TREE(ALIGN_TRIM.out)
    GENE_TREES(ALIGN_TRIM.out)
    ASTRAL(GENE_TREES.out.trees)
    PLOT(ASTRAL.out.tree.join(CONCAT_TREE.out.tree))
}
