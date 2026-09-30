#!/usr/bin/env python3
"""
Benchmark a query VCF against a truth VCF (e.g. GIAB) within high-confidence regions.

Both VCFs are restricted to PASS calls, split to biallelic records and left-normalised
against the reference (bcftools norm), then matched on CHROM/POS/REF/ALT with bcftools isec.
Metrics (TP/FP/FN, precision, recall, F1) are reported separately for SNVs and indels.

This is an allele-level exact-match comparison; representation-aware tools such as
hap.py/vcfeval can match a few additional complex indels, so indel scores here are conservative.
"""
import argparse, os, subprocess, tempfile, sys


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def count(vcf, vtype):
    out = subprocess.run(f"bcftools view -H -v {vtype} {vcf} | wc -l", shell=True,
                         check=True, capture_output=True, text=True, executable="/bin/bash")
    return int(out.stdout.strip())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--query", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bed", required=True, help="high-confidence regions BED")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--region", default="", help="optional region, e.g. 20:10000000-12000000")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    tmp = tempfile.mkdtemp(prefix="bench_")
    bed = os.path.join(tmp, "eval.bed")
    if a.region:
        chrom, rng = a.region.split(":")
        start, end = [int(x.replace(",", "")) for x in rng.split("-")]
        with open(a.bed) as fin, open(bed, "w") as fout:
            for line in fin:
                c, s, e = line.split("\t")[:3]
                s, e = int(s), int(e)
                if c == chrom and e > start - 1 and s < end:
                    fout.write(f"{c}\t{max(s, start - 1)}\t{min(e, end)}\n")
    else:
        run(f"cp {a.bed} {bed}")

    for name, src in (("query", a.query), ("truth", a.truth)):
        run(f"bcftools view -f PASS,. -T {bed} {src} "
            f"| bcftools norm -m -any -f {a.fasta} -c s "
            f"| bcftools view -e 'ALT=\"*\"' "
            f"| bcftools sort -Oz -o {tmp}/{name}.vcf.gz 2>/dev/null")
        run(f"bcftools index -t {tmp}/{name}.vcf.gz")

    # 0000 = private to query (FP), 0001 = private to truth (FN), 0002 = shared (TP, query records)
    run(f"bcftools isec -c none -p {tmp}/isec {tmp}/query.vcf.gz {tmp}/truth.vcf.gz")

    rows = []
    for label, vtype in (("SNV", "snps"), ("INDEL", "indels")):
        fp = count(f"{tmp}/isec/0000.vcf", vtype)
        fn = count(f"{tmp}/isec/0001.vcf", vtype)
        tp = count(f"{tmp}/isec/0002.vcf", vtype)
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        rows.append((label, tp, fp, fn, prec, rec, f1))

    bp = sum(int(l.split()[2]) - int(l.split()[1]) for l in open(bed))
    with open(a.out, "w") as f:
        f.write(f"# high-confidence bases evaluated: {bp}\n")
        f.write("type\tTP\tFP\tFN\tprecision\trecall\tF1\n")
        for r in rows:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.4f}\t{r[5]:.4f}\t{r[6]:.4f}\n")
    sys.stdout.write(open(a.out).read())


if __name__ == "__main__":
    main()
