/*
 * GATK-recommended hard filters (used when too few samples for VQSR/CNN scoring).
 * https://gatk.broadinstitute.org/hc/en-us/articles/360035890471
 */
process HARD_FILTER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/variants/filtered", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path ref

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${meta.id}.bcftools_stats.txt"),                                     emit: stats

    script:
    """
    gatk SelectVariants -R genome.fa -V ${vcf} --select-type-to-include SNP   -O snps.vcf.gz
    gatk SelectVariants -R genome.fa -V ${vcf} --select-type-to-include INDEL -O indels.vcf.gz

    gatk VariantFiltration -R genome.fa -V snps.vcf.gz -O snps.filt.vcf.gz \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "QUAL < 30.0" --filter-name "QUAL30" \\
        -filter "SOR > 3.0" --filter-name "SOR3" \\
        -filter "FS > 60.0" --filter-name "FS60" \\
        -filter "MQ < 40.0" --filter-name "MQ40" \\
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

    gatk VariantFiltration -R genome.fa -V indels.vcf.gz -O indels.filt.vcf.gz \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "QUAL < 30.0" --filter-name "QUAL30" \\
        -filter "FS > 200.0" --filter-name "FS200" \\
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

    gatk MergeVcfs -I snps.filt.vcf.gz -I indels.filt.vcf.gz -O ${meta.id}.filtered.vcf.gz
    bcftools stats -f PASS ${meta.id}.filtered.vcf.gz > ${meta.id}.bcftools_stats.txt
    """
}
