process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path ref

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam

    script:
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}_lib1\\tPU:${meta.id}.1"
    """
    bwa mem -t ${task.cpus} -K 100000000 -Y -R '${rg}' genome.fa ${reads[0]} ${reads[1]} \\
        | samtools sort -@ ${task.cpus} -m 1G -o ${meta.id}.sorted.bam -
    samtools index ${meta.id}.sorted.bam
    """
}
