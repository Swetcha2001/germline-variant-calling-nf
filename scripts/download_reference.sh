#!/usr/bin/env bash
# Download GRCh37 chr20 reference, annotation, ClinVar and GIAB HG001 truth set into ./ref
# Usage: bash scripts/download_reference.sh
set -euo pipefail
mkdir -p ref && cd ref

ENS=https://ftp.ensembl.org/pub
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh37
CLINVAR=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz

# 1. Reference: Ensembl GRCh37 chromosome 20 (contig name "20", same as b37/hs37d5)
[ -s Homo_sapiens.GRCh37.75.dna.chromosome.20.fa.gz ] || \
  curl -fsSL -o Homo_sapiens.GRCh37.75.dna.chromosome.20.fa.gz \
  $ENS/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.20.fa.gz

# 2. Gene models for bcftools csq (Ensembl GRCh37 release 87 GFF3, chr20)
[ -s Homo_sapiens.GRCh37.87.chr20.gff3.gz ] || \
  curl -fsSL -o Homo_sapiens.GRCh37.87.chr20.gff3.gz \
  $ENS/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chromosome.20.gff3.gz

# 3. ClinVar (GRCh37) restricted to chr20 via remote tabix query
if [ ! -s clinvar_GRCh37_chr20.vcf.gz ]; then
  bcftools view -r 20 $CLINVAR -Oz -o clinvar_GRCh37_chr20.vcf.gz
  bcftools index -t clinvar_GRCh37_chr20.vcf.gz
fi

# 4. GIAB HG001 v4.2.1 benchmark (truth VCF + high-confidence BED), chr20 only
if [ ! -s HG001_GRCh37_chr20_v4.2.1_benchmark.vcf.gz ]; then
  bcftools view -r 20 $GIAB/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \
    -Oz -o HG001_GRCh37_chr20_v4.2.1_benchmark.vcf.gz
  bcftools index -t HG001_GRCh37_chr20_v4.2.1_benchmark.vcf.gz
fi
[ -s HG001_GRCh37_chr20_v4.2.1_benchmark.bed ] || \
  curl -fsSL $GIAB/HG001_GRCh37_1_22_v4.2.1_benchmark.bed | awk '$1=="20"' > HG001_GRCh37_chr20_v4.2.1_benchmark.bed

rm -f clinvar.vcf.gz.tbi HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi
ls -lh
