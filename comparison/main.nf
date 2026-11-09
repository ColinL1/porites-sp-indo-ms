// ============================================================
// comparison/main.nf — Porites comparative analysis (funannotate compare)
// ============================================================
// Separate from the assembly / annotation pipeline (GenomeMaker), which it only imports from.
//
//   --step annotate_ref   Re-annotate the gene models of ONE annotated NCBI genome
//                         with the pipeline's own stage 08 (Phobius, InterProScan,
//                         eggNOG-mapper, funannotate annotate), so it carries the
//                         same evidence as the pipeline's genome.
//                         Needs --ref_fasta, --ref_gff (NCBI genomic.fna / .gff)
//                         and --species_name.
//
//   --step compare        funannotate compare over --compare_gbks (comma-separated
//                         funannotate GenBank files). --compare_proteinortho <tsv>
//                         hands in a precomputed ProteinOrtho table.
//
// Stage 08 is single-genome (fixed output names, params.species_name), so the
// reference genomes are annotated one run each rather than as one channel.
// Launch from a directory of its own, never the repo root, so that .nextflow/
// and work/ stay apart from the main pipeline's cache -- see
// run_compare.sh.
//
// projectDir is this folder: nextflow.config points the stage-08 databases at GenomeMaker and
// assets/ holds the three placeholder files the stage-08 subworkflow looks up in projectDir.
// ============================================================

nextflow.enable.dsl = 2

// Stage 08 and the protein extraction come from the GenomeMaker pipeline, reached through the
// link data/GenomeMaker (see ../data/README.md); nothing of the pipeline is copied or edited here.
include { PREP_REF_ANNOTATION; DROP_UNTRANSLATABLE; FUNANNOTATE_GENE_IDS; FUNANNOTATE_COMPARE } from './modules/10_comparison.nf'
include { EXTRACT_PROTEINS }      from '../data/GenomeMaker/modules/07_structural_annotation.nf'
include { FUNCTIONAL_ANNOTATION } from '../data/GenomeMaker/subworkflows/08_functional_annotation.nf'

params.step               = null
params.ref_fasta          = null
params.ref_gff            = null
params.compare_gbks       = null
// false only to resume a run made before DROP_UNTRANSLATABLE existed (P. harrisoni)
params.drop_untranslatable = true
params.compare_num_orthos = 500
params.compare_cpus       = 240
params.compare_proteinortho = null

workflow ANNOTATE_REF {

    main:
        if (!params.ref_fasta || !params.ref_gff) {
            error "--step annotate_ref needs --ref_fasta, --ref_gff and --species_name"
        }

        def funodb_ch   = channel.value(file(params.funannotate_db, checkIfExists: true))
        def eggnogdb_ch = channel.value(file(params.eggnog_db, checkIfExists: true))

        PREP_REF_ANNOTATION(
            channel.value(file(params.ref_fasta, checkIfExists: true)),
            channel.value(file(params.ref_gff, checkIfExists: true))
        )
        // `--drop_untranslatable false` on the command line arrives as the string 'false'
        def models_gff = PREP_REF_ANNOTATION.out.gff3
        if (params.drop_untranslatable && params.drop_untranslatable.toString() != 'false') {
            DROP_UNTRANSLATABLE(PREP_REF_ANNOTATION.out.gff3, PREP_REF_ANNOTATION.out.assembly)
            models_gff = DROP_UNTRANSLATABLE.out.gff3
        }
        EXTRACT_PROTEINS(models_gff, PREP_REF_ANNOTATION.out.assembly)
        // gene IDs only; transcript IDs (and so the proteins above) are unchanged
        FUNANNOTATE_GENE_IDS(models_gff)

        FUNCTIONAL_ANNOTATION(
            PREP_REF_ANNOTATION.out.assembly,
            FUNANNOTATE_GENE_IDS.out.gff3,
            EXTRACT_PROTEINS.out.proteins_faa,
            funodb_ch,
            eggnogdb_ch
        )
}

workflow COMPARE {

    main:
        if (!params.compare_gbks) {
            error "--step compare needs --compare_gbks (comma-separated funannotate .gbk files)"
        }
        def gbks = params.compare_gbks.toString().split(',')
            .collect { p -> p.trim() }.findAll { p -> p }
            .collect { p -> file(p, checkIfExists: true) }
        if (gbks.size() < 2) {
            error "funannotate compare needs at least two genomes, got ${gbks.size()}"
        }

        // optional: a ProteinOrtho table computed without -synteny (see the module)
        def po_tsv = params.compare_proteinortho
            ? file(params.compare_proteinortho, checkIfExists: true)
            : file("${projectDir}/assets/NO_FILE", checkIfExists: true)

        FUNANNOTATE_COMPARE(
            channel.value(gbks),
            channel.value(file(params.funannotate_db, checkIfExists: true)),
            channel.value(po_tsv)
        )
}

workflow {
    if (params.step == 'annotate_ref') {
        ANNOTATE_REF()
    }
    else if (params.step == 'compare') {
        COMPARE()
    }
    else {
        error "Set --step annotate_ref or --step compare (got: ${params.step})"
    }
}
