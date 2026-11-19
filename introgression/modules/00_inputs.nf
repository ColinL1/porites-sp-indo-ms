// modules/00_inputs.nf -- Inputs: the target assembly (indexed), the scaffold set and the resolved reference panel.
// Included by ../main.nf; results layout: the output block there.

process INPUTS {
    container params.tools_container

    input:
    path fasta, stageAs: 'input/*'
    path references, stageAs: 'references.in.tsv'
    path scaffolds
    val mito_dir                  // folder of the reference mitogenomes (@MITO_DIR@ in references.tsv)

    output:
    tuple path('target.fa'), path('target.fa.fai'), emit: target
    path 'scaffolds.bed', emit: bed
    path 'references.tsv', emit: references

    script:
    def bed = scaffolds.name != 'NO_SCAFFOLDS' ?
        "cp ${scaffolds} scaffolds.bed" :
        "awk -v m=${params.min_scaffold_len} '\$2>=m {print \$1\"\\t0\\t\"\$2}' target.fa.fai > scaffolds.bed"
    // @REF_DIR@ / @MITO_DIR@ in references.tsv -> this machine's folders
    def refdir = file(params.ref_dir).toString()
    def resolveCmd = "sed -e 's#@REF_DIR@#${refdir}#g' -e 's#@MITO_DIR@#${mito_dir}#g' ${references} > references.tsv"
    """
    ln -s \$(readlink -f ${fasta}) target.fa
    samtools faidx target.fa
    ${bed}
    [ -s scaffolds.bed ] || { echo "empty scaffold set" >&2; exit 1; }
    ${resolveCmd}
    """

    stub:
    def resolveCmd = "sed -e 's#@REF_DIR@#${file(params.ref_dir)}#g' -e 's#@MITO_DIR@#${mito_dir}#g' ${references} > references.tsv"
    """
    printf '>s1\\nACGT\\n' > target.fa; printf 's1\\t4\\t4\\t4\\t5\\n' > target.fa.fai
    printf 's1\\t0\\t4\\n' > scaffolds.bed; ${resolveCmd}
    """
}
