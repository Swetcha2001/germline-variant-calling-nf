/*
 * Functional consequence annotation with bcftools csq (Ensembl GFF3)
 * and clinical significance lookup against ClinVar.
 */
process ANNOTATE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/variants/annotated", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path ref
    path gff
    path clinvar
    path clinvar_tbi

    output:
    tuple val(meta), path("${meta.id}.annotated.vcf.gz"), path("${meta.id}.annotated.vcf.gz.tbi"), emit: vcf

    script:
    def run_csq = gff.name != 'NO_FILE'
    def run_clinvar = clinvar.name != 'NO_FILE'
    """
    bcftools view -f PASS ${vcf} -Oz -o pass.vcf.gz && bcftools index -t pass.vcf.gz
    cur=pass.vcf.gz

    if ${run_csq}; then
        bcftools csq -f genome.fa -g ${gff} --phase a -l \$cur -Oz -o csq.vcf.gz
        bcftools index -t csq.vcf.gz; cur=csq.vcf.gz
    fi

    if ${run_clinvar}; then
        bcftools annotate -a ${clinvar} \\
            -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/GENEINFO,INFO/CLNVID:=INFO/ALLELEID \\
            \$cur -Oz -o clin.vcf.gz
        bcftools index -t clin.vcf.gz; cur=clin.vcf.gz
    fi

    cp \$cur ${meta.id}.annotated.vcf.gz
    bcftools index -t -f ${meta.id}.annotated.vcf.gz
    """
}
