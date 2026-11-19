// Hybridization / introgression tests of one target assembly against the reference panel in
// references.tsv: is the target a hybrid of parent1 x parent2, or introgressed?
//
//   nextflow run main.nf -params-file params/Porites_sp_v1_2.yaml
//
// Processes: modules/ (one file per results step); their scripts: bin/. Results go to <outdir>/<NN_step>/, laid
// out by the output block at the end of this file. Tools run in exact-version BioContainers images; the multi-tool
// steps and Dsuite / HyDe / Twisst / genomics_general run in the custom image containers/introgression_tools.def
// (params.tools_container).

nextflow.enable.dsl = 2

include { INPUTS } from './modules/00_inputs.nf'
include { MAP_READS; CLAIR3; SNV_FILTER; MOSDEPTH } from './modules/01_heterozygosity.nf'
include { SPLIT_SNVS; PHASE_SCAFFOLD; PHASE_CONCAT; PHASE_STATS_HAPLOTAG } from './modules/02_phasing.nf'
include { ALIGN_REFERENCE } from './modules/03_reference_alignment.nf'
include { ANCESTRY_PAINTING } from './modules/04_ancestry_painting.nf'
include { JOINT_VCF } from './modules/joint_vcf.nf'
include { DSUITE } from './modules/05_d_statistics.nf'
include { GENO; TWISST; TWISST_PLOT } from './modules/06_twisst.nf'
include { ROBUSTNESS_TREES; TWISST_CELL; ROBUSTNESS_SUMMARY } from './modules/07_twisst_robustness.nf'
include { HYDE } from './modules/08_hyde.nf'
include { PCA; PLINK_BED; ADMIXTURE; ADMIXTURE_CV } from './modules/09_pca_admixture.nf'
include { MITO_READS; MITO_SUBSAMPLE; MITO_FLYE; MITO_SEQUENCES; MITO_ALIGN; MITO_TREE; MITO_SUMMARY } from './modules/10_mitogenome.nf'
include { SUMMARY } from './modules/summary.nf'

// ---------------------------------------------------------------- workflow
workflow {
    main:
    if (!(params.name ==~ /[A-Za-z0-9_]+/)) error "name must match [A-Za-z0-9_]+"
    ['fasta', 'reads', 'references'].each { p -> if (!params[p]) error "--${p} is required" }

    def refs_tsv  = file(params.references, checkIfExists: true)
    def refs      = channel.fromPath(refs_tsv).splitCsv(header: true, sep: '\t')
    def parents   = refs.toList().map { rows ->
        [rows.find { r -> r.role == 'parent1' }.name, rows.find { r -> r.role == 'parent2' }.name] }
    def outgroups = refs.filter { r -> r.role == 'outgroup' }.map { r -> r.name }
    def scaffolds = params.scaffolds ? file(params.scaffolds, checkIfExists: true) : file("${projectDir}/assets/NO_SCAFFOLDS")
    def agp       = params.agp ? file(params.agp, checkIfExists: true) : file("${projectDir}/assets/NO_AGP")
    def bait      = file(params.mito_bait, checkIfExists: true)
    def reads     = file(params.reads, checkIfExists: true)

    INPUTS(file(params.fasta, checkIfExists: true), refs_tsv, scaffolds, "${projectDir}/mitogenomes")
    def target = INPUTS.out.target.first()
    def bed    = INPUTS.out.bed.first()
    def refs_v = INPUTS.out.references.first()

    // 01-02: reads
    MAP_READS(target, reads)
    def bam = MAP_READS.out.bam.first()
    CLAIR3(target, bam, bed)
    def pileup = CLAIR3.out.first()
    SNV_FILTER(pileup, bed, MAP_READS.out.meandepth)
    MOSDEPTH(bam)
    SPLIT_SNVS(SNV_FILTER.out.snvs, bed)
    PHASE_SCAFFOLD(SPLIT_SNVS.out.flatten().map { f -> [f.name - '.in.vcf.gz', f] }, target, bam)
    PHASE_CONCAT(PHASE_SCAFFOLD.out.collect(), bed)
    PHASE_STATS_HAPLOTAG(PHASE_CONCAT.out, target, bam)

    // 03: every reference projected onto the target
    ALIGN_REFERENCE(refs.map { r -> [r.name, file(r.fasta.replace("@REF_DIR@", file(params.ref_dir).toString()), checkIfExists: true), r.preset] }, target)
    def aln = ALIGN_REFERENCE.out.files.flatten().collect()
    def aln_stats = ALIGN_REFERENCE.out.stats.collectFile(name: 'alignment_stats.tsv',
        seed: "name\tpreset\tvariants\tcallable_bp\n", sort: { f -> f.name })

    // 04
    ANCESTRY_PAINTING(refs_v, aln, target, PHASE_STATS_HAPLOTAG.out.bam)

    // 05-07: quartet joint VCF
    JOINT_VCF(channel.of('quartet', 'panel'), refs_v, aln, pileup, target)
    def quartet = JOINT_VCF.out.filter { v -> v[0] == 'quartet' }.first()
    def panel   = JOINT_VCF.out.filter { v -> v[0] == 'panel' }.first()
    def og_p    = outgroups.combine(parents).map { og, p1, p2 -> [og, p1, p2] }

    DSUITE(og_p.combine(quartet).map { og, p1, p2, _s, v, t -> [og, p1, p2, v, t] })
    GENO(quartet)
    def geno = GENO.out.first()
    TWISST(og_p.combine(geno).combine(refs_v))
    TWISST_PLOT(TWISST.out.weights.flatten().collect(), refs_v, target)

    def rob_og = params.robustness_outgroup ? channel.value(params.robustness_outgroup) : outgroups.first()
    ROBUSTNESS_TREES(geno, PHASE_CONCAT.out.first(), PHASE_STATS_HAPLOTAG.out.blocks, refs_v, rob_og)
    def cells = ROBUSTNESS_TREES.out.trees.flatten()
        .map { f -> [f.name - '_trees.nwk', f] }
        .combine(parents).combine(rob_og)
    TWISST_CELL(cells, ROBUSTNESS_TREES.out.meta.first())
    ROBUSTNESS_SUMMARY(TWISST_CELL.out.weights.flatten().collect(), MOSDEPTH.out.regions, aln, refs_v, agp)

    // 08-09: panel joint VCF
    HYDE(panel, refs_v)
    PCA(panel, refs_v)
    PLINK_BED(PCA.out.vcf)
    ADMIXTURE(PLINK_BED.out.combine(channel.of(2, 3, 4, 5)))
    ADMIXTURE_CV(ADMIXTURE.out.log.collect(), PLINK_BED.out)

    // 10: mitogenome, given or assembled from the reads
    def target_mito
    def mito_source
    def flye_info = channel.empty()
    if (params.mito) {
        target_mito = channel.value(file(params.mito, checkIfExists: true))
        mito_source = 'given'
    } else {
        MITO_READS(bait, reads)
        MITO_SUBSAMPLE(MITO_READS.out)
        MITO_FLYE(MITO_SUBSAMPLE.out)
        target_mito = MITO_FLYE.out.fasta
        flye_info   = MITO_FLYE.out.info
        mito_source = 'assembled from reads (Flye)'
    }
    MITO_SEQUENCES(target_mito, mito_source, bait, refs_v)
    MITO_ALIGN(MITO_SEQUENCES.out.fasta)
    MITO_TREE(MITO_ALIGN.out.aln)
    MITO_SUMMARY(MITO_TREE.out.tree)

    // summary/results.tsv: one row per test
    SUMMARY(
        SNV_FILTER.out.het,
        PHASE_STATS_HAPLOTAG.out.stats,
        ANCESTRY_PAINTING.out.summary,
        DSUITE.out.collect(),
        TWISST_PLOT.out.summary,
        ROBUSTNESS_SUMMARY.out.summary.mix(ROBUSTNESS_SUMMARY.out.context).collect(),
        HYDE.out.summary,
        PCA.out.distances.mix(ADMIXTURE_CV.out.cv).collect(),
        MITO_SEQUENCES.out.sources.mix(MITO_SUMMARY.out).collect(),
        refs_v)

    // what is published, by output name; where each output goes: the output block below
    publish:
    inputs              = INPUTS.out.bed.mix(INPUTS.out.references,
                              channel.value(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params)))
                                  .collectFile(name: 'params.json'))
    heterozygosity      = SNV_FILTER.out.het.mix(SNV_FILTER.out.snvs, MAP_READS.out.flagstat, MAP_READS.out.meandepth)
    depth               = MOSDEPTH.out.files
    phasing             = PHASE_CONCAT.out.mix(PHASE_STATS_HAPLOTAG.out.summary, PHASE_STATS_HAPLOTAG.out.stats,
                              PHASE_STATS_HAPLOTAG.out.blocks)
    alignment_stats     = aln_stats
    reference_alignment = ALIGN_REFERENCE.out.by_name
    ancestry_painting   = ANCESTRY_PAINTING.out.summary.mix(ANCESTRY_PAINTING.out.files)
    d_statistics        = DSUITE.out
    twisst              = TWISST_PLOT.out.summary.mix(TWISST_PLOT.out.plots)
    twisst_outgroup     = TWISST.out.by_outgroup
    robustness          = ROBUSTNESS_SUMMARY.out.summary.mix(ROBUSTNESS_SUMMARY.out.context,
                              ROBUSTNESS_SUMMARY.out.context_windows)
    robustness_grid     = ROBUSTNESS_TREES.out.trees.mix(ROBUSTNESS_TREES.out.meta).flatten()
    robustness_cell     = TWISST_CELL.out.by_cell
    hyde                = HYDE.out.summary
    hyde_raw            = HYDE.out.raw
    pca_admixture       = PCA.out.distances.mix(PCA.out.files, ADMIXTURE_CV.out.cv, ADMIXTURE_CV.out.samples)
    admixture           = ADMIXTURE.out.q.mix(ADMIXTURE.out.log)
    mitogenome          = MITO_SUMMARY.out.mix(MITO_SEQUENCES.out.fasta, MITO_SEQUENCES.out.sources, MITO_ALIGN.out.aln,
                              flye_info)
    mito_tree           = MITO_TREE.out.files.mix(MITO_ALIGN.out.log)
    summary             = SUMMARY.out
    joint_vcfs          = JOINT_VCF.out.map { _set, vcf, tbi -> [vcf, tbi] }
}

// ---------------------------------------------------------------- published results
// One entry per output of the publish: section above. Paths are relative to outputDir
// (params.outdir, nextflow.config); a closure gets one channel value and returns its folder.
// Change a folder here, nowhere else. Files are copied (workflow.output.mode).
output {
    inputs              { path '00_inputs' }
    heterozygosity      { path '01_heterozygosity' }
    depth               { path '01_heterozygosity/depth' }
    phasing             { path '02_phasing' }
    alignment_stats     { path '03_reference_alignment' }
    reference_alignment { path { name, _files -> "03_reference_alignment/${name}" } }
    ancestry_painting   { path '04_ancestry_painting' }
    d_statistics        { path '05_d_statistics' }
    twisst              { path '06_twisst' }
    twisst_outgroup     { path { og, _files -> "06_twisst/${og}" } }
    robustness          { path '07_twisst_robustness' }
    // <rep>_W<snps>_<NJ|ML>_* files go to their cell's folder, <rep>_groups.tsv to grid/
    robustness_grid     { path { f -> (f.name =~ /^([A-Za-z0-9]+_W[0-9]+_(NJ|ML))_/).with { m -> m ? "07_twisst_robustness/grid/${m[0][1]}" : '07_twisst_robustness/grid' } } }
    robustness_cell     { path { cell, _files -> "07_twisst_robustness/grid/${cell}" } }
    hyde                { path '08_hyde' }
    hyde_raw            { path '08_hyde/raw' }
    pca_admixture       { path '09_pca_admixture' }
    admixture           { path '09_pca_admixture/admixture' }
    mitogenome          { path '10_mitogenome' }
    mito_tree           { path '10_mitogenome/iqtree' }
    summary             { path '.' }
    // the joint VCFs of 05-09 (BAMs and other intermediates stay in the work directory)
    joint_vcfs          { path '_intermediate' }
}
