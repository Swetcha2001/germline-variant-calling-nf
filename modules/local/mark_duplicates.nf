process MARK_DUPLICATES {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.md.bam"), path("${meta.id}.md.bai"), emit: bam
    tuple val(meta), path("${meta.id}.md.metrics.txt"),                       emit: metrics

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" MarkDuplicates \\
        -I ${bam} -O ${meta.id}.md.bam -M ${meta.id}.md.metrics.txt \\
        --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT
    """
}
