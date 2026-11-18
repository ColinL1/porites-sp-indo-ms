// ============================================================
// compare.nf
// funannotate compare over the funannotate GenBank files of all genomes.
// ProteinOrtho runs with synteny by default (params.compare_synteny), inside funannotate compare. With
// compare_synteny = false the pipeline makes the table without -synteny (PROTEINORTHO_NOSYNTENY) and hands
// it to compare; params.compare_proteinortho supplies a ready-made table instead.
// ============================================================

nextflow.enable.dsl = 2

include { FUNANNOTATE_COMPARE; PROTEINORTHO_NOSYNTENY } from '../modules/10_comparison.nf'

workflow COMPARE {

    take:
        gbks            // channel of GenBank files (funannotate annotate output)
        funannotate_db  // value channel

    main:
        def sorted_gbks = gbks.toSortedList { a, b -> a.name <=> b.name }
        def no_file = file("${projectDir}/assets/NO_FILE", checkIfExists: true)

        if (params.compare_proteinortho) {
            po_tsv = channel.value(file(params.compare_proteinortho, checkIfExists: true))
        }
        else if (params.compare_synteny.toString().toBoolean()) {
            po_tsv = channel.value(no_file)      // compare runs proteinortho -synteny itself
        }
        else {
            PROTEINORTHO_NOSYNTENY(sorted_gbks)
            po_tsv = PROTEINORTHO_NOSYNTENY.out.tsv
        }

        FUNANNOTATE_COMPARE(sorted_gbks, funannotate_db, po_tsv)

    emit:
        compare_dir = FUNANNOTATE_COMPARE.out.compare_dir
        versions    = FUNANNOTATE_COMPARE.out.versions
}
