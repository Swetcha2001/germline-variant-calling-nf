process PREPARE_REFERENCE {
    tag "${fasta.name}"
    label 'process_medium'
    storeDir "${params.ref_cache}"

    input:
    path fasta

    output:
    path "genome.*", emit: ref

    script:
    def decompress = fasta.name.endsWith('.gz') ? "gzip -cd ${fasta} > genome.fa" : "cp -L ${fasta} genome.fa"
    """
    ${decompress}
    samtools faidx genome.fa
    gatk CreateSequenceDictionary -R genome.fa -O genome.dict
    bwa index genome.fa
    """
}
