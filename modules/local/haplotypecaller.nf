process HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/variants/raw", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref
    val intervals

    output:
    tuple val(meta), path("${meta.id}.raw.vcf.gz"), path("${meta.id}.raw.vcf.gz.tbi"), emit: vcf

    script:
    def iv = intervals ? "-L ${intervals}" : ''
    """
    gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" HaplotypeCaller \\
        -R genome.fa -I ${bam} -O ${meta.id}.raw.vcf.gz ${iv} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --standard-min-confidence-threshold-for-calling 20
    """
}
