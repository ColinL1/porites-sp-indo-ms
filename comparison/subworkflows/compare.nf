// ============================================================
// compare.nf
// funannotate compare over the funannotate GenBank files of all genomes.
// ProteinOrtho runs with synteny by default (params.compare_synteny). With compare_synteny = false
// a ProteinOrtho table made without -synteny must be handed in (params.compare_proteinortho).
// ============================================================

nextflow.enable.dsl = 2

include { FUNANNOTATE_COMPARE } from '../modules/10_comparison.nf'

workflow COMPARE {

    take:
        gbks            // channel of GenBank files (funannotate annotate output)
        funannotate_db  // value channel

    main:
        if (!params.compare_synteny && !params.compare_proteinortho) {
            error "compare_synteny = false needs compare_proteinortho (a ProteinOrtho table made without -synteny)"
        }
        // optional: a ProteinOrtho table computed without -synteny
        def po_tsv = params.compare_proteinortho
            ? file(params.compare_proteinortho, checkIfExists: true)
            : file("${projectDir}/assets/NO_FILE", checkIfExists: true)

        FUNANNOTATE_COMPARE(gbks.toSortedList { a, b -> a.name <=> b.name }, funannotate_db, channel.value(po_tsv))

    emit:
        compare_dir = FUNANNOTATE_COMPARE.out.compare_dir
        versions    = FUNANNOTATE_COMPARE.out.versions
}
