// ============================================================
// 10_comparison.nf
// Processes of the Porites comparative analysis (main.nf): reduce an annotated genome to bare
// gene models for re-annotation with the functional annotation (08_functional_annotation.nf),
// and run `funannotate compare`. Every process takes the genome name (sp) as an input, so all
// genomes run in one Nextflow run. Resources: process_resource.config.
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: PREP_REF_ANNOTATION
// Genome FASTA/GFF3 -> the inputs funannotate annotate expects: sequence names of at most
// 16 characters and a GFF3 holding only gene models (see bin/ncbi_ref_prep.pl).
// Isoforms are kept as delivered. Used for every genome, the new one included.
// ============================================================
process PREP_REF_ANNOTATION {
    tag   "${sp}"

    conda 'bioconda::agat=1.7.0'
    container 'quay.io/biocontainers/agat:1.7.0--pl5321hdfd78af_0'

    input:
    tuple val(sp), path(fasta), path(gff)

    output:
    tuple val(sp), path("${sp}.genome.fa"),   emit: assembly
    tuple val(sp), path("${sp}.models.gff3"), emit: gff3
    path "${sp}.seqmap.tsv", emit: seqmap
    path "${sp}.prep_summary.txt", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    ncbi_ref_prep.pl fasta ${fasta} ${sp}.genome.fa ${sp}.seqmap.tsv

    # Normalise first (complete gene/mRNA/exon/CDS hierarchy, one Parent per
    # line), then cut down to the models to re-annotate.
    agat_convert_sp_gxf2gxf.pl -g ${gff} -o normalised.gff3
    ncbi_ref_prep.pl gff normalised.gff3 ${sp}.seqmap.tsv \\
      > ${sp}.models.gff3 2> ${sp}.prep_summary.txt
    cat ${sp}.prep_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: "\$(agat_convert_sp_gxf2gxf.pl --version 2>&1 | grep -oim1 'v[0-9.]*' | sed 's/v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${sp}.genome.fa ${sp}.models.gff3 ${sp}.seqmap.tsv ${sp}.prep_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: DROP_UNTRANSLATABLE
// Remove mRNAs whose CDS does not translate cleanly from the genome (in-frame
// stop). RefSeq has a few hundred: models whose protein NCBI corrected relative
// to the assembly, and selenoproteins. gffread writes those stops as '.', which
// Phobius rejects outright, and funannotate would translate them the same way.
// Genes left without any transcript are removed with them.
// ============================================================
process DROP_UNTRANSLATABLE {
    tag   "${sp}"

    conda 'bioconda::gffread=0.12.7'
    container 'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4'

    input:
    tuple val(sp), path(gff3), path(assembly)

    output:
    tuple val(sp), path("${sp}.gff3"), emit: gff3
    path "${sp}.untranslatable.txt", emit: dropped
    path "versions.yml", emit: versions

    script:
    """
    gffread ${gff3} -g ${assembly} -y all_models.faa
    awk '/^>/{id=substr(\$1,2); next} /\\./{print id}' all_models.faa | sort -u > ${sp}.untranslatable.txt

    # ${gff3} carries only "ID=..;Parent=..[;product=..]" (bin/ncbi_ref_prep.pl).
    # pass 1: genes that keep a transcript; pass 2: print what survives.
    awk -F'\\t' '
      function attr(key,   n, a, i) {
          n = split(\$9, a, ";")
          for (i = 1; i <= n; i++) if (index(a[i], key "=") == 1) return substr(a[i], length(key) + 2)
          return ""
      }
      pass == 0 { bad[\$1] = 1; next }
      /^#/      { if (pass == 2) print; next }
      pass == 1 { if (\$3 != "gene" && \$3 != "exon" && \$3 != "CDS" && !(attr("ID") in bad)) keep[attr("Parent")] = 1; next }
      \$3 == "gene"                 { if (attr("ID") in keep) print; next }
      \$3 == "exon" || \$3 == "CDS" { if (!(attr("Parent") in bad)) print; next }
                                    { if (!(attr("ID") in bad)) print }
    ' pass=0 ${sp}.untranslatable.txt pass=1 ${gff3} pass=2 ${gff3} > ${sp}.gff3

    echo "mRNAs dropped (in-frame stop): \$(wc -l < ${sp}.untranslatable.txt)"
    echo "genes kept: \$(awk -F'\\t' '\$3=="gene"' ${sp}.gff3 | wc -l)"
    echo "mRNAs kept: \$(awk -F'\\t' '\$3=="mRNA"' ${sp}.gff3 | wc -l)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "\$(gffread --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${sp}.gff3 ${sp}.untranslatable.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: FUNANNOTATE_GENE_IDS
// Reshape gene IDs for `funannotate annotate`, which silently drops every Pfam/InterPro/eggNOG/GO
// annotation ("adding 0 valid annotations") when the gene ID does not fit the base it derives from
// each protein header "<transcript> <gene>". Rules and the two ID styles handled (NCBI "rna-...",
// funannotate "g1.t2"): bin/funannotate_gene_ids.sh. Transcript IDs are left alone.
// ============================================================
process FUNANNOTATE_GENE_IDS {
    tag   "${sp}"

    conda 'bioconda::gffread=0.12.7'
    container 'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4'

    input:
    tuple val(sp), path(gff3)

    output:
    tuple val(sp), path("${sp}.funannotate.gff3"), emit: gff3

    script:
    """
    funannotate_gene_ids.sh ${gff3} > ${sp}.funannotate.gff3
    """

    stub:
    """
    touch ${sp}.funannotate.gff3
    """
}

// ============================================================
// Process: PROTEINORTHO_NOSYNTENY
// The ProteinOrtho table for `funannotate compare` without -synteny (compare itself always adds it,
// and that step did not finish on fragmented assemblies, see FUNANNOTATE_COMPARE). Protein files
// are written by compare's own code (bin/gbk_to_proteinortho.py), then proteinortho runs with the
// call compare makes minus -synteny. The table goes back to compare with --proteinortho.
// ============================================================
process PROTEINORTHO_NOSYNTENY {
    tag   "proteinortho"

    conda 'bioconda::funannotate=1.8.17'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'

    input:
    path gbks   // sorted GenBank files, the order compare gets them in

    output:
    path "protortho_input/funannotate.proteinortho.tsv", emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    check_unique_locus_tags.sh ${gbks}

    gbk_to_proteinortho.py protortho_input ${gbks} > genome_names.txt

    cd protortho_input
    proteinortho -project=funannotate -cpus=${task.cpus} -singles -selfblast \
      \$(sed 's/\$/.faa/' ../genome_names.txt | tr '\\n' ' ')
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proteinortho: "\$(proteinortho --version 2>&1 | grep -oim1 '[0-9][0-9.]*' | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p protortho_input
    touch protortho_input/funannotate.proteinortho.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proteinortho: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: FUNANNOTATE_COMPARE
// Comparative genomics across funannotate-annotated genomes (GenBank files
// written by `funannotate annotate`).
// ============================================================
process FUNANNOTATE_COMPARE {
    tag   "funannotate_compare"

    conda 'bioconda::funannotate=1.8.17'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'

    input:
    path gbks
    path funannotate_db   // configured funannotate database (funannotate setup)
    path proteinortho_tsv // precomputed ProteinOrtho table, or a NO_ placeholder

    output:
    path "funannotate_compare", emit: compare_dir
    path "versions.yml",        emit: versions

    script:
    // Without a table, compare runs `proteinortho -synteny` itself (params.compare_synteny). On
    // fragmented assemblies that synteny (FFAdj-MCS) step did not finish: after 2.5 h only the
    // 3 pairs involving the chromosome-level P. lutea were done, one job had found 190 million
    // candidate runs and the 12 jobs held 340 GB. A table made without -synteny can be handed in
    // instead (params.compare_proteinortho).
    def po_arg = proteinortho_tsv.name.startsWith('NO_') ? '' : "--proteinortho ${proteinortho_tsv}"
    // --ml_method raxml: the image ships RAxML but not IQ-TREE (compare's default).
    // No --run_dnds: codeml is not in the image either.
    // No --outgroup: FUNANNOTATE_DB/outgroups only holds fungal BUSCO sets.
    """
    export FUNANNOTATE_DB=\$(readlink -f ${funannotate_db})

    funannotate compare \\
      -i ${gbks} \\
      -o funannotate_compare \\
      --ml_method raxml \\
      ${po_arg} \\
      --num_orthos ${params.compare_num_orthos} \\
      --cpus ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "\$(funannotate --version 2>&1 | grep -o '[0-9][0-9.]*\$' | tail -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p funannotate_compare
    echo "FUNANNOTATE_COMPARE stub completed" > funannotate_compare/index.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "stub"
    END_VERSIONS
    """
}
