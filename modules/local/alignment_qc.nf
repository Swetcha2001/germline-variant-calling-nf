process ALIGNMENT_QC {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/qc/alignment", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref

    output:
    tuple val(meta), path("${meta.id}.{stats,flagstat,idxstats}.txt"), emit: stats

    script:
    """
    samtools stats --reference genome.fa ${bam} > ${meta.id}.stats.txt
    samtools flagstat ${bam} > ${meta.id}.flagstat.txt
    samtools idxstats ${bam} > ${meta.id}.idxstats.txt
    """
}
