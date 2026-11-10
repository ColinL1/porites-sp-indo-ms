// ============================================================
// annotate_ref.nf
// Every genome (the new one and the NCBI references alike) is reduced to bare gene models and
// re-annotated with the same functional evidence, so that funannotate compare sees comparable
// annotations:
//   prep (names <= 16 chars, gene models only) -> drop untranslatable mRNAs -> proteins
//   -> funannotate-safe gene IDs -> eggNOG-mapper (+ InterProScan, Phobius) -> funannotate annotate
// ============================================================

nextflow.enable.dsl = 2

include { PREP_REF_ANNOTATION; DROP_UNTRANSLATABLE; FUNANNOTATE_GENE_IDS } from '../modules/10_comparison.nf'
include { EXTRACT_PROTEINS } from '../modules/07_extract_proteins.nf'
include { PHOBIUS_RUN; INTERPROSCAN; EGGNOG_MAPPER; FUNANNOTATE_ANNOTATE } from '../modules/08_functional_annotation.nf'

workflow ANNOTATE_REF {

    take:
        genomes         // [sp, fasta, gff]
        funannotate_db  // value channel
        eggnog_db       // value channel

    main:
        // distinct placeholder names: one shared NO_FILE causes an input file name collision
        def no_ipr = file("${projectDir}/assets/NO_IPRSCAN.xml", checkIfExists: true)
        def no_pho = file("${projectDir}/assets/NO_PHOBIUS.txt", checkIfExists: true)

        PREP_REF_ANNOTATION(genomes)

        // drop_untranslatable: false only to resume an annotation made before this step existed
        def models = PREP_REF_ANNOTATION.out.gff3
        if (params.drop_untranslatable) {
            DROP_UNTRANSLATABLE(PREP_REF_ANNOTATION.out.gff3.join(PREP_REF_ANNOTATION.out.assembly))
            models = DROP_UNTRANSLATABLE.out.gff3
        }
        EXTRACT_PROTEINS(models.join(PREP_REF_ANNOTATION.out.assembly))
        // gene IDs only; transcript IDs (and so the proteins above) are unchanged
        FUNANNOTATE_GENE_IDS(models)

        EGGNOG_MAPPER(EXTRACT_PROTEINS.out.proteins_faa, eggnog_db)

        // --- optional: InterProScan (local install, or container + matching data) ---
        def ipr_ch = EXTRACT_PROTEINS.out.proteins_faa.map { sp, _p -> tuple(sp, no_ipr) }
        def ipr_local     = params.interproscan_dir && file("${params.interproscan_dir}/interproscan.sh").exists()
        def ipr_container = !params.interproscan_dir && params.interproscan_data && file(params.interproscan_data).exists()
        if (ipr_local || ipr_container) {
            INTERPROSCAN(
                EXTRACT_PROTEINS.out.proteins_faa,
                channel.value(file(ipr_local ? params.interproscan_dir : params.interproscan_data))
            )
            ipr_ch = INTERPROSCAN.out.iprscan_xml
        } else {
            log.warn "InterProScan SKIPPED (interproscan_dir=${params.interproscan_dir}, interproscan_data=${params.interproscan_data}) -- " +
                     "the annotations carry NO InterPro domain assignments and none of the GO terms derived from them"
        }

        // --- optional: Phobius ---
        def pho_ch = EXTRACT_PROTEINS.out.proteins_faa.map { sp, _p -> tuple(sp, no_pho) }
        if (params.run_phobius) {
            PHOBIUS_RUN(
                EXTRACT_PROTEINS.out.proteins_faa,
                channel.value(file(params.phobius_dir, checkIfExists: true))
            )
            pho_ch = PHOBIUS_RUN.out.phobius_results
        } else {
            log.warn "Phobius SKIPPED (run_phobius=false) -- no signal-peptide or transmembrane predictions"
        }

        FUNANNOTATE_ANNOTATE(
            FUNANNOTATE_GENE_IDS.out.gff3
                .join(PREP_REF_ANNOTATION.out.assembly)
                .join(EGGNOG_MAPPER.out.eggnog_anno)
                .join(ipr_ch)
                .join(pho_ch),
            funannotate_db
        )

    emit:
        gbk          = FUNANNOTATE_ANNOTATE.out.gbk          // [sp, <sp>.gbk]
        prepared     = PREP_REF_ANNOTATION.out.gff3
        proteins     = EXTRACT_PROTEINS.out.proteins_faa
        eggnog       = EGGNOG_MAPPER.out.eggnog_anno
        iprscan      = ipr_ch
        phobius      = pho_ch
        annotate_dir = FUNANNOTATE_ANNOTATE.out.anno_dir
        stats_json   = FUNANNOTATE_ANNOTATE.out.stats_json
        seqmap       = PREP_REF_ANNOTATION.out.seqmap
        prep_summary = PREP_REF_ANNOTATION.out.summary
}
