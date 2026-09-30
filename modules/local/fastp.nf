process FASTP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.trim.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"),           emit: json
    tuple val(meta), path("${meta.id}.fastp.html"),           emit: html

    script:
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${meta.id}_R1.trim.fastq.gz -O ${meta.id}_R2.trim.fastq.gz \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred 15 --length_required 36 \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
}
