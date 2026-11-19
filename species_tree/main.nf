// BUSCO single-copy phylogeny of the new Porites assembly, two taxon sets (small, large) x
// two approaches (concatenated ML, ASTRAL-IV) = 4 trees, plus the quartet and tanglegram plots.
//
//   nextflow run main.nf -params-file params/Porites_sp_v1_2.yaml
// The params file gives name and either fasta (BUSCO is run on it) or busco_dir
// (an existing BUSCO genome-mode output).
//
// Steps as in the earlier runs (earlier single_ref and astral_expanded runs, see archive/), except
// ASTRAL-IV is run unconstrained (the earlier runs passed the concatenated tree with -c, which in
// ASTRAL-IV is a topology constraint).

// Processes: modules/ (one file per results step); their scripts: bin/. Tools: one pinned container per process
// (Apptainer); -profile conda uses the same versions as conda packages.

nextflow.enable.dsl = 2

include { BUSCO_NEW } from './modules/01_busco.nf'
include { SHARED_LOCI } from './modules/02_loci.nf'
include { ALIGN_TRIM } from './modules/03_alignments.nf'
include { CONCAT_TREE } from './modules/04_concat.nf'
include { GENE_TREES } from './modules/05_gene_trees.nf'
include { ASTRAL } from './modules/06_astral.nf'
include { PLOT } from './modules/07_plots.nf'

workflow {

    main:
        if (!params.name) {
            error "Give the tip label of the new genome with name (params file)"
        }
        if (!(params.fasta || params.busco_dir)) {
            error "Give the new assembly with --fasta, or its BUSCO output with --busco_dir"
        }
        def new_busco = params.busco_dir
            ? channel.value([params.name, file(params.busco_dir, checkIfExists: true)])
            : BUSCO_NEW(channel.of([params.name, file(params.fasta, checkIfExists: true)])).first()

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

    // what is published, by output name; where each output goes: the output block below
    publish:
    inputs      = channel.fromPath(params.references).collectFile(name: 'references.tsv')   // a copy: files outside work/ are not published
                      .mix(channel.value(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params)))
                               .collectFile(name: 'params.json'))
    busco       = params.busco_dir ? channel.empty() : new_busco
    loci        = SHARED_LOCI.out.genes
    alignments  = ALIGN_TRIM.out
    concat      = CONCAT_TREE.out.files
    gene_trees  = GENE_TREES.out.trees.join(GENE_TREES.out.files)
    astral      = ASTRAL.out.tree.join(ASTRAL.out.log)
    plots       = PLOT.out
}

// ---------------------------------------------------------------- published results
// One entry per output of the publish: section above. Paths are relative to outputDir
// (params.outdir, nextflow.config); a closure gets one channel value and returns its folder.
// Change a folder here, nowhere else. Files are copied (workflow.output.mode).
output {
    inputs      { path '00_inputs' }
    busco       { path '01_busco' }
    loci        { path { taxa, _files -> "02_loci/${taxa}" } }
    alignments  { path { taxa, _files -> "03_alignments/${taxa}" } }
    concat      { path { taxa, _f1, _f2, _f3 -> "04_concat/${taxa}" } }
    gene_trees  { path { taxa, _f1, _f2 -> "05_gene_trees/${taxa}" } }
    astral      { path { taxa, _f1, _f2 -> "06_astral/${taxa}" } }
    plots       { path { taxa, _files -> "07_plots/${taxa}" } }
}
