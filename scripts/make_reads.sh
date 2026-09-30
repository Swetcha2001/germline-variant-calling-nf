#!/usr/bin/env bash
# Extract paired-end reads for a region from the public GIAB HG001 (NA12878) 30x Illumina BAM
# and write them as FASTQ, so the pipeline starts from raw reads.
# Usage: bash scripts/make_reads.sh 20:10000000-12000000 benchmark_data
set -euo pipefail
REGION=${1:?region e.g. 20:10000000-12000000}
OUT=${2:?output dir}
BAM=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam
TAG=$(echo "$REGION" | sed 's/[:,-]/_/g')
mkdir -p "$OUT"

# properly paired, primary, non-supplementary reads only; collate by name, then split mates
samtools view -u -f 2 -F 0x900 "$BAM" "$REGION" \
  | samtools collate -Ou - "$OUT/tmp_$TAG" \
  | samtools fastq -N -1 "$OUT/NA12878_${TAG}_R1.fastq.gz" -2 "$OUT/NA12878_${TAG}_R2.fastq.gz" \
                   -0 /dev/null -s /dev/null -

printf "sample,fastq_1,fastq_2\nNA12878,%s,%s\n" \
  "$OUT/NA12878_${TAG}_R1.fastq.gz" "$OUT/NA12878_${TAG}_R2.fastq.gz" > "$OUT/samplesheet.csv"
echo "Wrote $(( $(zcat "$OUT/NA12878_${TAG}_R1.fastq.gz" | wc -l) / 4 )) read pairs to $OUT"
