/*
 * Compare calls to a GIAB truth set inside high-confidence regions.
 * Variants are split into biallelic records, left-normalised against the reference
 * and matched on CHROM/POS/REF/ALT (bcftools isec), reported separately for SNVs and indels.
 */
process BENCHMARK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/benchmark", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path ref
    path truth_vcf
    path truth_tbi
    path truth_bed
    val intervals

    output:
    tuple val(meta), path("${meta.id}.benchmark.tsv"), emit: metrics

    script:
    """
    benchmark_vcf.py --query ${vcf} --truth ${truth_vcf} --bed ${truth_bed} \\
        --fasta genome.fa --region '${intervals}' --out ${meta.id}.benchmark.tsv
    """
}
