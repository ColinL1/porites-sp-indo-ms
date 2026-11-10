// ============================================================
// comparison/main.nf -- the new Porites genome against other annotated genomes
// ============================================================
//
//   nextflow run main.nf -params-file params/Porites_sp_v1_2.yaml
//   nextflow run main.nf -params-file params/Porites_sp_v1_2.yaml -profile stub -stub    (wiring test)
//   nextflow run main.nf -params-file params/Porites_sp_v1_2_binac2.yaml -profile binac2 (BinAC 2)
//
// Every genome listed under `genomes:` in the params file (the new one by fasta + gff, the NCBI
// references by accession or by fasta + gff) goes through the same steps, in parallel:
//   NCBI download (accession only) -> gene-model prep -> drop untranslatable mRNAs -> proteins
//   -> eggNOG-mapper / InterProScan / Phobius -> funannotate annotate
// so that every genome carries the same functional evidence. `funannotate compare` then runs once
// over all GenBank files. Results go to <outdir>/<NN_step>/, laid out by the output block at the
// end of this file.
//
// Stage 07/08 modules are copies of GenomeMaker's (modules/), with the genome name as an input;
// nothing is imported from GenomeMaker at run time.

nextflow.enable.dsl = 2

include { NCBI_DOWNLOAD } from './modules/ncbi_download.nf'
include { ANNOTATE_REF } from './subworkflows/annotate_ref.nf'
include { COMPARE } from './subworkflows/compare.nf'

workflow {

    main:
        def genomes = (params.genomes ?: []) as List
        def extra_gbks = (params.compare_gbks ?: []) as List

        // ---- validate the genome list before anything starts
        def names = [] as List
        genomes.each { g ->
            if (!g.name || !(g.name.toString() ==~ /[A-Za-z0-9_]+/)) {
                error "each genome needs a `name` of letters, digits and '_' (funannotate -s name), got: ${g}"
            }
            if (!g.accession && !(g.fasta && g.gff)) {
                error "genome ${g.name}: give either `accession` or both `fasta` and `gff`"
            }
            if (g.name in names) {
                error "genome name used twice: ${g.name}"
            }
            names << g.name
        }
        if (genomes.size() + extra_gbks.size() < 2) {
            error "funannotate compare needs at least two genomes (params.genomes plus params.compare_gbks), got ${genomes.size() + extra_gbks.size()}"
        }

        def funodb_ch   = channel.value(file(params.funannotate_db, checkIfExists: true))
        def eggnogdb_ch = channel.value(file(params.eggnog_db, checkIfExists: true))

        // ---- inputs: paths as given, or downloaded by accession
        def local_ch = channel.fromList(genomes.findAll { g -> g.fasta && g.gff })
            .map { g -> tuple(g.name.toString(), file(g.fasta, checkIfExists: true), file(g.gff, checkIfExists: true)) }
        def download_ch = channel.fromList(genomes.findAll { g -> !(g.fasta && g.gff) })
            .map { g -> tuple(g.name.toString(), g.accession.toString()) }
        NCBI_DOWNLOAD(download_ch)
        def genomes_ch = local_ch.mix(NCBI_DOWNLOAD.out.genome)

        // ---- same annotation for every genome
        ANNOTATE_REF(genomes_ch, funodb_ch, eggnogdb_ch)

        // ---- compare: the new annotations plus any ready-made funannotate GenBank files
        def gbk_ch = ANNOTATE_REF.out.gbk
            .map { _sp, gbk -> gbk }
            .mix(channel.fromList(extra_gbks).map { p -> file(p, checkIfExists: true) })
        COMPARE(gbk_ch, funodb_ch)

    // what is published, by output name; where each output goes: the output block below
    publish:
    inputs       = channel.value(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params)))
                       .collectFile(name: 'params.json')
    ncbi         = NCBI_DOWNLOAD.out.genome
    prepared     = ANNOTATE_REF.out.prepared
    seqmap       = ANNOTATE_REF.out.seqmap.mix(ANNOTATE_REF.out.prep_summary)
    proteins     = ANNOTATE_REF.out.proteins
    eggnog       = ANNOTATE_REF.out.eggnog
    iprscan      = ANNOTATE_REF.out.iprscan
    phobius      = ANNOTATE_REF.out.phobius
    annotation   = ANNOTATE_REF.out.annotate_dir
    stats        = ANNOTATE_REF.out.stats_json
    compare      = COMPARE.out.compare_dir
}

// ---------------------------------------------------------------- published results
// One entry per output of the publish: section above. Paths are relative to outputDir
// (params.outdir, nextflow.config); a closure gets one channel value and returns its folder.
// Change a folder here, nowhere else. Files are copied (workflow.output.mode).
output {
    inputs       { path '00_inputs' }
    ncbi         { path '00_inputs/ncbi' }
    prepared     { path { sp, _files -> "01_prepared/${sp}" } }
    seqmap       { path '01_prepared' }
    proteins     { path { sp, _files -> "02_proteins/${sp}" } }
    eggnog       { path { sp, _files -> "03_functional_annotation/${sp}/eggnog" } }
    iprscan      { path { sp, _files -> "03_functional_annotation/${sp}/interproscan" } }
    phobius      { path { sp, _files -> "03_functional_annotation/${sp}/phobius" } }
    annotation   { path { sp, _files -> "03_functional_annotation/${sp}/funannotate" } }
    stats        { path { sp, _files -> "03_functional_annotation/${sp}/funannotate" } }
    compare      { path '04_compare' }
}
