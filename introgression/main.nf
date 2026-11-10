// Hybridization / introgression tests of one target assembly against the reference panel in
// references.tsv: is the target a hybrid of parent1 x parent2, or introgressed?
//
//   nextflow run main.nf -params-file params/Porites_sp_v1_2.yaml
//
// Same steps and scripts (bin/) as run_all.sh; what each step tests: ../docs/. Results go to
// <outdir>/<NN_step>/, laid out by the output block at the end of this file. Tools run in exact-version BioContainers images; the
// multi-tool steps and Dsuite / HyDe / Twisst / genomics_general run in the custom image
// containers/introgression_tools.def (params.tools_container).

nextflow.enable.dsl = 2


// ---------------------------------------------------------------- 00 inputs
process INPUTS {
    container params.tools_container

    input:
    path fasta, stageAs: 'input/*'
    path references, stageAs: 'references.in.tsv'
    path scaffolds

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
    def resolveCmd = "sed -e 's#@REF_DIR@#${refdir}#g' -e 's#@MITO_DIR@#${projectDir}/mitogenomes#g' ${references} > references.tsv"
    """
    ln -s \$(readlink -f ${fasta}) target.fa
    samtools faidx target.fa
    ${bed}
    [ -s scaffolds.bed ] || { echo "empty scaffold set" >&2; exit 1; }
    ${resolveCmd}
    """

    stub:
    def resolveCmd = "sed -e 's#@REF_DIR@#${file(params.ref_dir)}#g' -e 's#@MITO_DIR@#${projectDir}/mitogenomes#g' ${references} > references.tsv"
    """
    printf '>s1\\nACGT\\n' > target.fa; printf 's1\\t4\\t4\\t4\\t5\\n' > target.fa.fai
    printf 's1\\t0\\t4\\n' > scaffolds.bed; ${resolveCmd}
    """
}

// ---------------------------------------------------------------- 01 heterozygosity
process MAP_READS {
    tag "${params.new_name}"
    container params.tools_container
    cpus params.max_cpus

    input:
    tuple path(target), path(fai)
    path reads

    output:
    tuple path('reads.bam'), path('reads.bam.bai'), emit: bam
    path 'reads.meandepth.txt', emit: meandepth
    path 'reads.flagstat.txt', emit: flagstat

    script:
    """
    map_reads.sh ${target} ${reads} ${params.new_name} reads.bam ${task.cpus}
    """

    stub:
    """
    touch reads.bam reads.bam.bai reads.flagstat.txt; printf 'mean_depth_all_positions\\t18.7\\ttotal_pos\\t4\\n' > reads.meandepth.txt
    """
}

process CLAIR3 {
    tag "${params.new_name}"
    container 'quay.io/biocontainers/clair3:2.0.3--py311hbc58adc_0'
    cpus params.max_cpus

    input:
    tuple path(target), path(fai)
    tuple path(bam), path(bai)
    path bed

    output:
    tuple path('clair3/pileup.vcf.gz'), path('clair3/pileup.vcf.gz.tbi')

    script:
    """
    run_clair3.sh --bam_fn=${bam} --ref_fn=${target} --threads=${task.cpus} --platform=ont \\
      --model_path=${params.clair3_model} --output=clair3 --bed_fn=${bed} --min_coverage=4 \\
      --sample_name=${params.new_name}
    """

    stub:
    """
    mkdir clair3; touch clair3/pileup.vcf.gz clair3/pileup.vcf.gz.tbi
    """
}

process SNV_FILTER {
    tag "${params.new_name}"
    container params.tools_container

    input:
    tuple path(pileup), path(tbi)
    path bed
    path meandepth

    output:
    tuple path('snvs.vcf.gz'), path('snvs.vcf.gz.tbi'), emit: snvs
    path 'heterozygosity.tsv', emit: het

    script:
    """
    bcftools view -f PASS -v snps -m2 -M2 -e 'QUAL<10' ${pileup} -Oz -o snvs.vcf.gz
    tabix -f -p vcf snvs.vcf.gz
    het=\$(bcftools view -H -g het snvs.vcf.gz | wc -l)
    hom=\$(bcftools view -H -g hom snvs.vcf.gz | wc -l)
    bp=\$(awk '{s+=\$3-\$2} END{print s}' ${bed})
    printf "het_snvs\\t%s\\nhom_alt_snvs\\t%s\\nscaffold_bp\\t%s\\nhet_per_bp\\t%s\\nmean_depth\\t%s\\n" \\
      \$het \$hom \$bp \$(awk -v h=\$het -v b=\$bp 'BEGIN{printf "%.6f", h/b}') \$(cut -f2 ${meandepth}) > heterozygosity.tsv
    """

    stub:
    """
    touch snvs.vcf.gz snvs.vcf.gz.tbi
    printf 'het_snvs\\t1\\nhom_alt_snvs\\t1\\nscaffold_bp\\t4\\nhet_per_bp\\t0.25\\nmean_depth\\t18.7\\n' > heterozygosity.tsv
    """
}

process MOSDEPTH {
    tag "${params.new_name}"
    container 'quay.io/biocontainers/mosdepth:0.3.14--h87be163_2'
    cpus 4

    input:
    tuple path(bam), path(bai)

    output:
    path 'depth.regions.bed.gz', emit: regions
    path 'depth.*', emit: files

    script:
    """
    mosdepth -t ${task.cpus} -n -b 5000 depth ${bam}
    """

    stub:
    """
    touch depth.regions.bed.gz depth.mosdepth.summary.txt
    """
}

// ---------------------------------------------------------------- 02 phasing
process SPLIT_SNVS {
    container params.tools_container

    input:
    tuple path(snvs), path(tbi)
    path bed

    output:
    path 'per/*.in.vcf.gz'

    script:
    """
    mkdir per
    cut -f1 ${bed} | xargs -P ${task.cpus} -I{} sh -c 'bcftools view -r {} ${snvs} -Oz -o per/{}.in.vcf.gz'
    """

    stub:
    """
    mkdir per; cut -f1 ${bed} | while read c; do touch per/\$c.in.vcf.gz; done
    """
}

process PHASE_SCAFFOLD {
    tag "${chrom}"
    container 'quay.io/biocontainers/whatshap:2.8--py311he264feb_0'
    cpus 1

    input:
    tuple val(chrom), path(vcf)
    tuple path(target), path(fai)
    tuple path(bam), path(bai)

    output:
    path "${chrom}.phased.vcf.gz"

    script:
    """
    if [ "\$(zcat ${vcf} | grep -vc '^#')" -lt 2 ]; then
      cp ${vcf} ${chrom}.phased.vcf.gz
    else
      whatshap phase --ignore-read-groups --chromosome ${chrom} -o ${chrom}.phased.vcf.gz \\
        --reference ${target} ${vcf} ${bam}
    fi
    """

    stub:
    """
    touch ${chrom}.phased.vcf.gz
    """
}

process PHASE_CONCAT {
    container params.tools_container

    input:
    path vcfs
    path bed

    output:
    tuple path('phased.vcf.gz'), path('phased.vcf.gz.tbi')

    script:
    """
    for f in ${vcfs}; do tabix -f -p vcf \$f; done
    bcftools concat -Oz -o phased.vcf.gz \$(cut -f1 ${bed} | sed 's/\$/.phased.vcf.gz/')
    tabix -f -p vcf phased.vcf.gz
    """

    stub:
    """
    touch phased.vcf.gz phased.vcf.gz.tbi
    """
}

process PHASE_STATS_HAPLOTAG {
    tag "${params.new_name}"
    container 'quay.io/biocontainers/whatshap:2.8--py311he264feb_0'
    cpus 8

    input:
    tuple path(phased), path(tbi)
    tuple path(target), path(fai)
    tuple path(bam), path(bai)

    output:
    tuple path('haplotagged.bam'), path('haplotagged.bam.bai'), emit: bam
    path 'phase_block_list.tsv', emit: blocks
    path 'phase_block_stats.tsv', emit: stats
    path 'phase_stats.txt', emit: summary

    script:
    """
    whatshap stats ${phased} --tsv=phase_block_stats.tsv --block-list=phase_block_list.tsv \\
      --chr-lengths <(cut -f1,2 ${fai}) > phase_stats.txt 2>&1
    whatshap haplotag --ignore-read-groups --skip-missing-contigs -o haplotagged.bam --reference ${target} \\
      ${phased} ${bam}
    python -c "import pysam; pysam.index('-@', '${task.cpus}', 'haplotagged.bam')"
    """

    stub:
    """
    touch haplotagged.bam haplotagged.bam.bai phase_block_list.tsv phase_stats.txt
    printf 'sample\\tchromosome\\tblocks\\tphased_fraction\\tblock_n50\\nT\\tALL\\t1\\t0.99\\t21000\\n' > phase_block_stats.tsv
    """
}

// ---------------------------------------------------------------- 03 reference alignment
process ALIGN_REFERENCE {
    tag "${name}"
    container params.tools_container
    cpus 48

    input:
    tuple val(name), path(fna), val(preset)
    tuple path(target), path(fai)

    output:
    path "${name}.{callable.bed,vcf.gz,vcf.gz.tbi}", emit: files
    tuple val(name), path("${name}.{callable.bed,vcf.gz,vcf.gz.tbi}"), emit: by_name
    path "${name}.stats.tsv", emit: stats

    script:
    """
    mkdir work
    align_reference.sh ${name} ${fna} ${preset} ${target} . work ${task.cpus}
    mv work/${name}.stats.tsv .
    """

    stub:
    """
    touch ${name}.callable.bed ${name}.vcf.gz ${name}.vcf.gz.tbi; printf '${name}\\t${preset}\\t1\\t1\\n' > ${name}.stats.tsv
    """
}

// ---------------------------------------------------------------- 04 ancestry painting
process ANCESTRY_PAINTING {
    tag "${params.new_name}"
    container params.tools_container
    cpus params.max_cpus

    input:
    path references
    path aln, stageAs: 'aln/*'
    tuple path(target), path(fai)
    tuple path(bam), path(bai)

    output:
    path 'ancestry_summary.txt', emit: summary
    path '*.{tsv,tsv.gz,png}', emit: files

    script:
    """
    diagnostic_sites.py --references ${references} --aln_dir aln --target ${target} --out diagnostic_sites.tsv.gz
    ancestry_painting.py --diag diagnostic_sites.tsv.gz --bam ${bam} --fai ${fai} --references ${references} \\
      --outdir . --threads ${task.cpus}
    """

    stub:
    """
    touch diagnostic_sites.tsv.gz ancestry_painting.png
    printf 'diagnostic sites: 1\\nsites with >= 1 haplotype-resolved call: 1\\n  HOM_P1 1 (50.0%%)\\n  HOM_P2 1 (50.0%%)\\nallele matches genome-wide: P1=1  P2=1  P1_fraction=0.5\\nphase blocks with >= 10 resolved sites: 1\\n' > ancestry_summary.txt
    """
}

// ---------------------------------------------------------------- joint VCFs
process JOINT_VCF {
    tag "${set}"
    container params.tools_container

    input:
    val set
    path references
    path aln, stageAs: 'aln/*'
    tuple path(pileup), path(pileup_tbi)
    tuple path(target), path(fai)

    output:
    tuple val(set), path("joint_${set}.vcf.gz"), path("joint_${set}.vcf.gz.tbi")

    script:
    """
    joint_vcf.py --set ${set} --references ${references} --target_name ${params.new_name} --aln_dir aln \\
      --clair3_vcf ${pileup} --target ${target} --out joint_${set}.vcf
    bgzip -f joint_${set}.vcf; tabix -f -p vcf joint_${set}.vcf.gz
    """

    stub:
    """
    touch joint_${set}.vcf.gz joint_${set}.vcf.gz.tbi
    """
}

// ---------------------------------------------------------------- 05 D-statistics
process DSUITE {
    tag "${og}"
    container params.tools_container

    input:
    tuple val(og), val(p1), val(p2), path(vcf), path(tbi)

    output:
    path "${og}"

    script:
    """
    mkdir ${og}
    bcftools view -s ${p1},${p2},${params.new_name},${og} ${vcf} -Oz -o quartet.vcf.gz
    printf "%s\\t%s\\n%s\\t%s\\n%s\\t%s\\n%s\\tOutgroup\\n" ${p1} ${p1} ${p2} ${p2} ${params.new_name} ${params.new_name} ${og} > ${og}/SETS.txt
    cd ${og} && Dsuite Dtrios ../quartet.vcf.gz SETS.txt > dtrios.log 2>&1
    """

    stub:
    """
    mkdir ${og}; touch ${og}/SETS.txt ${og}/dtrios.log
    printf 'P1\\tP2\\tP3\\tDstatistic\\tZ-score\\tp-value\\tf4-ratio\\tBBAA\\tABBA\\tBABA\\n${p1}\\t${p2}\\t${params.new_name}\\t0.01\\t1\\t0.3\\t0.001\\t3\\t1\\t1\\n' > ${og}/SETS_BBAA.txt
    """
}

// ---------------------------------------------------------------- 06 Twisst
process GENO {
    container params.tools_container

    input:
    tuple val(set), path(vcf), path(tbi)

    output:
    path 'joint_quartet.geno.gz'

    script:
    """
    python \$PARSEVCF -i ${vcf} -o joint_quartet.geno.gz
    """

    stub:
    """
    touch joint_quartet.geno.gz
    """
}

process TWISST {
    tag "${og}"
    container params.tools_container
    cpus 1

    input:
    tuple val(og), val(p1), val(p2), path(geno), path(references)

    output:
    path "${og}_{weights.csv.gz,windows.tsv}", emit: weights
    tuple val(og), path("${og}_*"), emit: by_outgroup

    script:
    """
    twisst_trees.py --geno ${geno} --references ${references} --target_name ${params.new_name} --outgroup ${og} \\
      --window 50 --prefix ${og}
    python \$TWISST -t ${og}_trees.nwk -w ${og}_weights_per_tree.csv.gz --outputTopos ${og}_topos.nwk \\
      -g ${p1} -g ${p2} -g ${params.new_name} -g ${og} --method complete --groupsFile ${og}_groups.tsv --silent
    twisst_average.py ${og}_weights_per_tree.csv.gz ${og}_tree_window.txt ${og}_weights.csv.gz
    """

    stub:
    """
    touch ${og}_weights.csv.gz ${og}_windows.tsv ${og}_trees.nwk
    """
}

process TWISST_PLOT {
    container params.tools_container

    input:
    path weights, stageAs: 'tw/*'
    path references
    tuple path(target), path(fai)

    output:
    path 'twisst_summary.tsv', emit: summary
    path '*.png', emit: plots

    script:
    """
    twisst_plot.py --dir tw --references ${references} --target_name ${params.new_name} --fai ${fai}
    mv tw/twisst_summary.tsv tw/*.png .
    """

    stub:
    """
    printf 'outgroup\\tn_windows\\tT1_x\\tT2_x\\tT3_x\\tplurality\\n' > twisst_summary.tsv; touch twisst_painting.png
    """
}

// ---------------------------------------------------------------- 07 Twisst robustness
process ROBUSTNESS_TREES {
    container params.tools_container
    cpus params.max_cpus

    input:
    path geno
    tuple path(phased), path(phased_tbi)
    path blocks
    path references
    val og

    output:
    path '*_trees.nwk', emit: trees
    path '*_{tree_window.txt,windows.tsv,groups.tsv}', emit: meta

    script:
    """
    twisst_robustness.py --geno ${geno} --phased_vcf ${phased} --block_list ${blocks} --references ${references} \\
      --target_name ${params.new_name} --outgroup ${og} --outdir . --threads ${task.cpus}
    """

    stub:
    """
    for m in rand2 phased2 cons1; do touch \${m}_groups.tsv; for w in 50 200; do for t in NJ ML; do
      touch \${m}_W\${w}_\${t}_trees.nwk \${m}_W\${w}_\${t}_tree_window.txt \${m}_W\${w}_\${t}_windows.tsv; done; done; done
    """
}

process TWISST_CELL {
    tag "${cell}"
    container params.tools_container
    cpus 1

    input:
    tuple val(cell), path(trees), val(p1), val(p2), val(og)
    path meta

    output:
    path "${cell}_{weights.csv.gz,windows.tsv}", includeInputs: true, emit: weights
    tuple val(cell), path("${cell}_{weights.csv.gz,weights_per_tree.csv.gz,topos.nwk}"), emit: by_cell

    script:
    def mode = cell.tokenize('_')[0]
    """
    python \$TWISST -t ${cell}_trees.nwk -w ${cell}_weights_per_tree.csv.gz --outputTopos ${cell}_topos.nwk \\
      -g ${p1} -g ${p2} -g ${params.new_name} -g ${og} --method complete --groupsFile ${mode}_groups.tsv --silent
    twisst_average.py ${cell}_weights_per_tree.csv.gz ${cell}_tree_window.txt ${cell}_weights.csv.gz
    """

    stub:
    """
    touch ${cell}_weights.csv.gz ${cell}_weights_per_tree.csv.gz ${cell}_topos.nwk
    """
}

process ROBUSTNESS_SUMMARY {
    container params.tools_container

    input:
    path weights, stageAs: 'grid/*'
    path depth
    path aln, stageAs: 'aln/*'
    path references
    path agp

    output:
    path 'robustness_summary.tsv', emit: summary
    path 'context_summary.txt', emit: context
    path 'context.tsv', emit: context_windows

    script:
    def agp_opt = agp.name != 'NO_AGP' ? "--agp ${agp}" : ''
    """
    robustness_summary.py --dir grid --references ${references} --target_name ${params.new_name}
    context.py --dir grid --depth ${depth} --aln_dir aln --references ${references} --target_name ${params.new_name} ${agp_opt}
    mv grid/robustness_summary.tsv grid/context.tsv grid/context_summary.txt .
    """

    stub:
    """
    printf 'representation\\twindow_snps\\tmethod\\tn_windows\\tT1\\tT2\\tT3\\tplurality\\n' > robustness_summary.tsv
    printf 'context\\tr_T1\\tr_T2\\tr_T3\\ndepth\\t+0.010\\t+0.010\\t+0.010\\n' > context_summary.txt; touch context.tsv
    """
}

// ---------------------------------------------------------------- 08 HyDe
process HYDE {
    container params.tools_container
    cpus 4

    input:
    tuple val(set), path(vcf), path(tbi)
    path references

    output:
    path 'hyde_summary.tsv', emit: summary
    path '{seqs.txt,map.txt,hyde_*-out*.txt}', emit: raw

    script:
    """
    hyde.py --vcf ${vcf} --references ${references} --target_name ${params.new_name} --outdir .
    """

    stub:
    """
    printf 'outgroup\\tP1\\tHybrid\\tP2\\tZscore\\tPvalue\\tGamma\\n' > hyde_summary.tsv; touch seqs.txt map.txt
    """
}

// ---------------------------------------------------------------- 09 PCA / ADMIXTURE
process PCA {
    container params.tools_container

    input:
    tuple val(set), path(vcf), path(tbi)
    path references

    output:
    tuple path('panel_filtered.vcf.gz'), path('panel_filtered.vcf.gz.tbi'), emit: vcf
    path 'pca_distances.txt', emit: distances
    path 'pca.{tsv,png}', emit: files

    script:
    """
    bcftools view -m2 -M2 -v snps -e 'F_MISSING>0.10' -q 0.05:minor ${vcf} -Oz -o panel_filtered.vcf.gz
    tabix -f -p vcf panel_filtered.vcf.gz
    pca_admixture.py --vcf panel_filtered.vcf.gz --references ${references} --target_name ${params.new_name} \\
      --outdir . --pca_only
    """

    stub:
    """
    touch panel_filtered.vcf.gz panel_filtered.vcf.gz.tbi pca.tsv pca.png
    printf 'x\\nT distance\\n  lobata 1.0\\n  harrisoni 1.0\\nlobata cluster radius: 0.4; T is outside it\\n' > pca_distances.txt
    """
}

process PLINK_BED {
    container 'quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0'
    cpus 4

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path('ad.bed'), path('ad.bim'), path('ad.fam')

    script:
    """
    plink2 --vcf ${vcf} --set-all-var-ids @:# --allow-extra-chr --double-id --max-alleles 2 \\
      --threads ${task.cpus} --make-bed --out ad
    awk 'BEGIN{OFS="\\t"} {\$1="1"; print}' ad.bim > ad.bim.tmp && mv ad.bim.tmp ad.bim
    """

    stub:
    """
    touch ad.bed ad.bim ad.fam
    """
}

process ADMIXTURE {
    tag "K=${k}"
    container 'quay.io/biocontainers/admixture:1.3.0--0'
    cpus 1

    input:
    tuple path(bed), path(bim), path(fam), val(k)

    output:
    path "admixture_K${k}.log", emit: log
    path "ad.${k}.Q", emit: q
    path "ad.${k}.P"

    script:
    """
    admixture --cv ${bed} ${k} > admixture_K${k}.log 2>&1
    """

    stub:
    """
    echo 'CV error (K=${k}): 0.${k}' > admixture_K${k}.log; touch ad.${k}.Q ad.${k}.P
    """
}

process ADMIXTURE_CV {
    container params.tools_container

    input:
    path logs
    tuple path(bed), path(bim), path(fam)

    output:
    path 'admixture_cv.tsv', emit: cv
    path 'admixture_samples.txt', emit: samples

    script:
    """
    { printf "K\\tcv_error\\n"; grep -h "CV error" ${logs} | sed -E 's/.*K=([0-9]+)\\): */\\1\\t/' | sort -n; } > admixture_cv.tsv
    awk '{print \$2}' ${fam} > admixture_samples.txt
    """
}

// ---------------------------------------------------------------- 10 mitogenome
process MITO_READS {
    container params.tools_container
    cpus params.max_cpus

    input:
    path bait
    path reads

    output:
    path 'mito_reads.fastq.gz'

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${bait} ${reads} 2>/dev/null | samtools view -b -F 0x904 - | samtools sort -o mito.bam
    samtools fastq mito.bam | gzip > mito_reads.fastq.gz
    """

    stub:
    """
    touch mito_reads.fastq.gz
    """
}

process MITO_SUBSAMPLE {
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    input:
    path fq

    output:
    path 'mito_sub.fastq.gz'

    script:
    """
    seqkit sample -s11 -p 0.15 ${fq} -o mito_sub.fastq.gz
    """

    stub:
    """
    touch mito_sub.fastq.gz
    """
}

process MITO_FLYE {
    container 'quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1'
    cpus 16

    input:
    path fq

    output:
    path 'assembly.fasta', emit: fasta
    path 'flye_assembly_info.txt', emit: info

    script:
    """
    flye --nano-hq ${fq} --out-dir flye --threads ${task.cpus} --min-overlap 1000
    cp flye/assembly.fasta assembly.fasta; cp flye/assembly_info.txt flye_assembly_info.txt
    """

    stub:
    """
    touch assembly.fasta flye_assembly_info.txt
    """
}

process MITO_SEQUENCES {
    container params.tools_container
    cpus 16

    input:
    path target_mito
    val source
    path bait
    path references

    output:
    path 'mito_seqs.fasta', emit: fasta
    path 'mito_sources.tsv', emit: sources

    script:
    """
    mito_sequences.py --references ${references} --target_name ${params.new_name} --target_mito ${target_mito} \\
      --target_source "${source}" --bait ${bait} --outdir . --threads ${task.cpus}
    """

    stub:
    """
    printf '>a\\nA\\n' > mito_seqs.fasta; printf 'sample\\tsource\\tlength\\n${params.new_name}\\tgiven\\t18000\\n' > mito_sources.tsv
    """
}

process MITO_ALIGN {
    container 'quay.io/biocontainers/mafft:7.525--h031d066_1'
    cpus 8

    input:
    path seqs

    output:
    path 'mito_aln.fasta', emit: aln
    path 'mafft.log', emit: log

    script:
    """
    mafft --auto --adjustdirection --thread ${task.cpus} ${seqs} > mito_aln.fasta 2> mafft.log
    sed -i 's/^>_R_/>/' mito_aln.fasta
    """

    stub:
    """
    cp ${seqs} mito_aln.fasta; touch mafft.log
    """
}

process MITO_TREE {
    container 'quay.io/biocontainers/iqtree:3.1.3--h8471819_0'
    cpus 16

    input:
    path aln

    output:
    path 'mito_tree.treefile', emit: tree
    path 'mito_tree.*', emit: files

    script:
    """
    iqtree -s ${aln} -m MFP -B 1000 -seed 1 -T ${task.cpus} --prefix mito_tree -redo > iqtree.log 2>&1
    """

    stub:
    """
    echo '(a,b,c);' > mito_tree.treefile; touch mito_tree.iqtree
    """
}

process MITO_SUMMARY {
    container params.tools_container

    input:
    path tree

    output:
    path 'mito_tree_summary.txt'

    script:
    """
    mito_tree_summary.py --tree ${tree} --root ${params.mito_root} --target_name ${params.new_name} > mito_tree_summary.txt
    """

    stub:
    """
    echo '${params.new_name} sister group: x (support 50)' > mito_tree_summary.txt
    """
}

// ---------------------------------------------------------------- summary
process SUMMARY {
    container params.tools_container

    input:
    path het,     stageAs: 'r/01_heterozygosity/*'
    path phase,   stageAs: 'r/02_phasing/*'
    path paint,   stageAs: 'r/04_ancestry_painting/*'
    path dstats,  stageAs: 'r/05_d_statistics/*'
    path twisst,  stageAs: 'r/06_twisst/*'
    path robust,  stageAs: 'r/07_twisst_robustness/*'
    path hyde,    stageAs: 'r/08_hyde/*'
    path pca,     stageAs: 'r/09_pca_admixture/*'
    path mito,    stageAs: 'r/10_mitogenome/*'
    path references

    output:
    path 'summary/results.tsv'

    script:
    """
    summary.py --outdir r --references ${references} --target_name ${params.new_name}
    mkdir summary; cp r/summary/results.tsv summary/
    """

    stub:
    """
    mkdir summary; printf 'test\\tresult\\tfolder\\n' > summary/results.tsv
    """
}

// ---------------------------------------------------------------- workflow
workflow {
    main:
    if (!(params.new_name ==~ /[A-Za-z0-9_]+/)) error "new_name must match [A-Za-z0-9_]+"
    ['new_fasta', 'reads', 'references'].each { if (!params[it]) error "--${it} is required" }

    def refs_tsv  = file(params.references, checkIfExists: true)
    def refs      = channel.fromPath(refs_tsv).splitCsv(header: true, sep: '\t')
    def parents   = refs.toList().map { rows ->
        [rows.find { it.role == 'parent1' }.name, rows.find { it.role == 'parent2' }.name] }
    def outgroups = refs.filter { it.role == 'outgroup' }.map { it.name }
    def scaffolds = params.scaffolds ? file(params.scaffolds, checkIfExists: true) : file("${projectDir}/assets/NO_SCAFFOLDS")
    def agp       = params.agp ? file(params.agp, checkIfExists: true) : file("${projectDir}/assets/NO_AGP")
    def bait      = file(params.mito_bait, checkIfExists: true)
    def reads     = file(params.reads, checkIfExists: true)

    INPUTS(file(params.new_fasta, checkIfExists: true), refs_tsv, scaffolds)
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
    PHASE_SCAFFOLD(SPLIT_SNVS.out.flatten().map { [it.name - '.in.vcf.gz', it] }, target, bam)
    PHASE_CONCAT(PHASE_SCAFFOLD.out.collect(), bed)
    PHASE_STATS_HAPLOTAG(PHASE_CONCAT.out, target, bam)

    // 03: every reference projected onto the target
    ALIGN_REFERENCE(refs.map { [it.name, file(it.fasta.replace("@REF_DIR@", file(params.ref_dir).toString()), checkIfExists: true), it.preset] }, target)
    def aln = ALIGN_REFERENCE.out.files.flatten().collect()
    def aln_stats = ALIGN_REFERENCE.out.stats.collectFile(name: 'alignment_stats.tsv',
        seed: "name\tpreset\tvariants\tcallable_bp\n", sort: { it.name })

    // 04
    ANCESTRY_PAINTING(refs_v, aln, target, PHASE_STATS_HAPLOTAG.out.bam)

    // 05-07: quartet joint VCF
    JOINT_VCF(channel.of('quartet', 'panel'), refs_v, aln, pileup, target)
    def quartet = JOINT_VCF.out.filter { it[0] == 'quartet' }.first()
    def panel   = JOINT_VCF.out.filter { it[0] == 'panel' }.first()
    def og_p    = outgroups.combine(parents).map { og, p1, p2 -> [og, p1, p2] }

    DSUITE(og_p.combine(quartet).map { og, p1, p2, s, v, t -> [og, p1, p2, v, t] })
    GENO(quartet)
    def geno = GENO.out.first()
    TWISST(og_p.combine(geno).combine(refs_v))
    TWISST_PLOT(TWISST.out.weights.flatten().collect(), refs_v, target)

    def rob_og = params.robustness_outgroup ? channel.value(params.robustness_outgroup) : outgroups.first()
    ROBUSTNESS_TREES(geno, PHASE_CONCAT.out.first(), PHASE_STATS_HAPLOTAG.out.blocks, refs_v, rob_og)
    def cells = ROBUSTNESS_TREES.out.trees.flatten()
        .map { [it.name - '_trees.nwk', it] }
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
    if (params.new_mito) {
        target_mito = channel.value(file(params.new_mito, checkIfExists: true))
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
