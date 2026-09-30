process VARIANT_SUMMARY {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.variant_table.tsv"), path("${meta.id}.variant_summary.md"), emit: report

    script:
    """
    summarize_variants.py --vcf ${vcf} --sample ${meta.id} \\
        --table ${meta.id}.variant_table.tsv --report ${meta.id}.variant_summary.md
    """
}
