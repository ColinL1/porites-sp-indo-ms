#!/usr/bin/env bash
# Porites comparative analysis: re-annotate the NCBI references listed in
# genomes.tsv with the pipeline's stage 08, then run funannotate compare against
# the pipeline's own annotation.
#
#   ./run_compare.sh annotate [species_name ...]   all of genomes.tsv, or the named ones
#   ./run_compare.sh compare
#   ./run_compare.sh all
#
# GENOMES_TSV picks the genome list and COMPARE_NAME the compare run, so that more
# than one comparison can be kept side by side (results in compare_<name>/):
#   ./run_compare.sh compare                                   # Porites only
#   GENOMES_TSV=genomes_outgroup.tsv COMPARE_NAME=outgroup ./run_compare.sh compare
#
# The NCBI download uses the `ncbi_datasets` conda env (DATASETS_ENV); every
# other step runs in the pipeline's containers.
#
# Every Nextflow run is launched from its own directory under comparison/, so
# its .nextflow/ cache and work/ never mix with the main pipeline's in the repo
# root. Re-running a step resumes it.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROFILE="${PROFILE:-local}"
DATASETS_ENV="${DATASETS_ENV:-ncbi_datasets}"
# extra nextflow arguments for the annotate step, e.g. "--drop_untranslatable false"
EXTRA_ARGS="${EXTRA_ARGS:-}"
GENOMES_TSV="${GENOMES_TSV:-genomes.tsv}"
COMPARE_NAME="${COMPARE_NAME:-porites}"
# the new genome's annotation (GenomeMaker results, linked in ../data/)
OWN_GBK="${OWN_GBK:-$HERE/../data/results_v1.2/08_functional_annotation/funannotate/final_anno/annotate_results/Porites_sp.gbk}"

ref_gbk() {
    echo "$HERE/refs/$1/08_functional_annotation/funannotate/final_anno/annotate_results/$1.gbk"
}

# Genome FASTA + GFF3 for one accession -> refs/<species>/00_ncbi/
fetch() {
    local acc="$1" dir="$HERE/refs/$2/00_ncbi"
    if [ -s "$dir/$acc.genomic.fna" ] && [ -s "$dir/$acc.genomic.gff" ]; then
        return
    fi
    mkdir -p "$dir/download"
    (cd "$dir/download" \
        && conda run -n "$DATASETS_ENV" datasets download genome accession "$acc" \
            --include genome,gff3 --no-progressbar --filename ncbi_dataset.zip \
        && unzip -q -o ncbi_dataset.zip)
    if [ ! -s "$dir/download/ncbi_dataset/data/$acc/genomic.gff" ]; then
        echo "ERROR: $acc has no GFF3 annotation on NCBI" >&2
        exit 1
    fi
    mv "$dir/download/ncbi_dataset/data/$acc/"*_genomic.fna "$dir/$acc.genomic.fna"
    mv "$dir/download/ncbi_dataset/data/$acc/genomic.gff" "$dir/$acc.genomic.gff"
}

genomes() {
    grep -v '^#' "$HERE/$GENOMES_TSV" | awk -F'\t' 'NF >= 2'
}

annotate() {
    local wanted=("$@")
    while IFS=$'\t' read -r acc sp; do
        if [ "${#wanted[@]}" -gt 0 ] && [[ ! " ${wanted[*]} " == *" $sp "* ]]; then
            continue
        fi
        echo "=== annotate_ref: $sp ($acc) ==="
        mkdir -p "$HERE/refs/$sp"
        fetch "$acc" "$sp"
        # </dev/null: keep nextflow from consuming the genome list on stdin
        (cd "$HERE/refs/$sp" && nextflow run "$HERE/main.nf" \
            -profile "$PROFILE" \
            --step annotate_ref \
            --ref_fasta "$HERE/refs/$sp/00_ncbi/$acc.genomic.fna" \
            --ref_gff "$HERE/refs/$sp/00_ncbi/$acc.genomic.gff" \
            --species_name "$sp" \
            --outdir "$HERE/refs/$sp" \
            $EXTRA_ARGS \
            -resume </dev/null)
    done < <(genomes)
}

compare() {
    local gbks="$OWN_GBK"
    while IFS=$'\t' read -r acc sp; do
        local gbk
        gbk="$(ref_gbk "$sp")"
        if [ ! -s "$gbk" ]; then
            echo "ERROR: $sp ($acc) is not annotated yet: $gbk" >&2
            exit 1
        fi
        gbks="$gbks,$gbk"
    done < <(genomes)
    echo "=== compare ($COMPARE_NAME) ==="
    mkdir -p "$HERE/compare_$COMPARE_NAME"
    (cd "$HERE/compare_$COMPARE_NAME" && nextflow run "$HERE/main.nf" \
        -profile "$PROFILE" \
        --step compare \
        --compare_gbks "$gbks" \
        --outdir "$HERE/compare_$COMPARE_NAME" \
        $EXTRA_ARGS \
        -resume </dev/null)
}

case "${1:-}" in
    annotate) shift; annotate "$@" ;;
    compare)  compare ;;
    all)      annotate; compare ;;
    *) echo "usage: $0 annotate [species_name ...] | compare | all" >&2; exit 2 ;;
esac
