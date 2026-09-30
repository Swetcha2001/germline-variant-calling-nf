#!/usr/bin/env nextflow
/*
 * germline-variant-calling-nf
 * Short-read germline SNV/indel calling following GATK best practices:
 *   FastQC -> fastp -> BWA-MEM -> MarkDuplicates -> (BQSR) -> HaplotypeCaller
 *   -> hard filtering -> bcftools csq consequence annotation -> ClinVar lookup
 *   -> optional benchmarking against a GIAB truth set -> MultiQC
 */
nextflow.enable.dsl = 2

include { FASTQC                  } from './modules/local/fastqc'
include { FASTP                   } from './modules/local/fastp'
include { PREPARE_REFERENCE       } from './modules/local/prepare_reference'
include { BWA_MEM                 } from './modules/local/bwa_mem'
include { MARK_DUPLICATES         } from './modules/local/mark_duplicates'
include { BQSR                    } from './modules/local/bqsr'
include { ALIGNMENT_QC            } from './modules/local/alignment_qc'
include { HAPLOTYPECALLER         } from './modules/local/haplotypecaller'
include { HARD_FILTER             } from './modules/local/hard_filter'
include { ANNOTATE                } from './modules/local/annotate'
include { VARIANT_SUMMARY         } from './modules/local/variant_summary'
include { BENCHMARK               } from './modules/local/benchmark'
include { MULTIQC                 } from './modules/local/multiqc'

def helpMessage() {
    log.info """
    germline-variant-calling-nf
    ===========================
    Usage:
      nextflow run main.nf -profile test,conda
      nextflow run main.nf --input samples.csv --fasta ref.fa [options] -profile docker

    Required:
      --input           CSV with columns: sample,fastq_1,fastq_2
      --fasta           Reference genome FASTA (plain or .gz)

    Optional:
      --intervals       BED/interval string to restrict calling (e.g. 20:10000000-12000000)
      --known_sites     VCF of known variants (dbSNP/Mills); enables BQSR
      --gff             Ensembl GFF3 for bcftools csq consequence annotation
      --clinvar         bgzipped + indexed ClinVar VCF for clinical significance lookup
      --truth_vcf       GIAB truth VCF (bgzipped + indexed) for benchmarking
      --truth_bed       GIAB high-confidence BED for benchmarking
      --outdir          Output directory [${params.outdir}]
    """.stripIndent()
}

// Relative paths in the samplesheet are resolved against the launch dir, then the project dir
def resolvePath(String p) {
    def f = file(p)
    if (f.exists()) return f
    def g = file("${projectDir}/${p}")
    if (g.exists()) return g
    error "File not found: ${p}"
}

workflow {
    if (params.help) { helpMessage(); exit 0 }
    if (!params.input) { error "Please provide --input samplesheet (see --help)" }
    if (!params.fasta) { error "Please provide --fasta reference (see --help)" }

    ch_reads = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            [meta, [resolvePath(row.fastq_1), resolvePath(row.fastq_2)]]
        }

    // Reference: decompress, index (bwa, samtools faidx, GATK dict) once
    PREPARE_REFERENCE(file(params.fasta, checkIfExists: true))
    ch_ref = PREPARE_REFERENCE.out.ref.collect()   // [fasta, fai, dict, bwa index files...]

    FASTQC(ch_reads)
    FASTP(ch_reads)
    BWA_MEM(FASTP.out.reads, ch_ref)
    MARK_DUPLICATES(BWA_MEM.out.bam)

    if (params.known_sites) {
        ks     = file(params.known_sites, checkIfExists: true)
        ks_tbi = file("${params.known_sites}.tbi", checkIfExists: true)
        BQSR(MARK_DUPLICATES.out.bam, ch_ref, ks, ks_tbi)
        ch_final_bam = BQSR.out.bam
    } else {
        ch_final_bam = MARK_DUPLICATES.out.bam
    }

    ALIGNMENT_QC(ch_final_bam, ch_ref)
    HAPLOTYPECALLER(ch_final_bam, ch_ref, params.intervals ?: '')
    HARD_FILTER(HAPLOTYPECALLER.out.vcf, ch_ref)

    gff     = params.gff     ? file(params.gff, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")
    clinvar = params.clinvar ? file(params.clinvar, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")
    clinvar_tbi = params.clinvar ? file("${params.clinvar}.tbi", checkIfExists: true) : file("${projectDir}/assets/NO_FILE_TBI")
    ANNOTATE(HARD_FILTER.out.vcf, ch_ref, gff, clinvar, clinvar_tbi)
    VARIANT_SUMMARY(ANNOTATE.out.vcf)

    ch_bench = Channel.empty()
    if (params.truth_vcf && params.truth_bed) {
        BENCHMARK(
            HARD_FILTER.out.vcf,
            ch_ref,
            file(params.truth_vcf, checkIfExists: true),
            file("${params.truth_vcf}.tbi", checkIfExists: true),
            file(params.truth_bed, checkIfExists: true),
            params.intervals ?: ''
        )
        ch_bench = BENCHMARK.out.metrics
    }

    ch_mqc = FASTQC.out.zip.map { it[1] }
        .mix(FASTP.out.json.map { it[1] })
        .mix(MARK_DUPLICATES.out.metrics.map { it[1] })
        .mix(ALIGNMENT_QC.out.stats.map { it[1] })
        .mix(HARD_FILTER.out.stats.map { it[1] })
        .collect()
    MULTIQC(ch_mqc)
}

workflow.onComplete {
    log.info(workflow.success ? "\nDone. Results in: ${params.outdir}\n" : "\nPipeline failed - see .nextflow.log\n")
}
