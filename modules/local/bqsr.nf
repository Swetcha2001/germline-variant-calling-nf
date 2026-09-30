process BQSR {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref
    path known_sites
    path known_sites_tbi

    output:
    tuple val(meta), path("${meta.id}.recal.bam"), path("${meta.id}.recal.bai"), emit: bam
    tuple val(meta), path("${meta.id}.recal.table"),                           emit: table

    script:
    """
    gatk BaseRecalibrator -R genome.fa -I ${bam} --known-sites ${known_sites} -O ${meta.id}.recal.table
    gatk ApplyBQSR -R genome.fa -I ${bam} --bqsr-recal-file ${meta.id}.recal.table -O ${meta.id}.recal.bam
    """
}
