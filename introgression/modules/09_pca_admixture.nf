// modules/09_pca_admixture.nf -- PCA and ADMIXTURE on the panel joint VCF.
// Included by ../main.nf; results layout: the output block there.

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
    pca_admixture.py --vcf panel_filtered.vcf.gz --references ${references} --target_name ${params.name} \\
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
