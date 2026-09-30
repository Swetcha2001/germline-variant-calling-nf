#!/usr/bin/env python3
"""
Turn an annotated VCF (bcftools csq + ClinVar) into a flat variant table and a short
Markdown summary: variant counts, Ti/Tv, consequence classes, genes hit by protein-altering
variants, and any ClinVar pathogenic / likely pathogenic matches.
"""
import argparse, collections, subprocess


def query(vcf, fmt):
    cmd = ["bcftools", "query", "-f", fmt, vcf]
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.splitlines()


def has_tag(vcf, tag):
    hdr = subprocess.run(["bcftools", "view", "-h", vcf], check=True, capture_output=True, text=True).stdout
    return f"ID={tag}," in hdr


TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
PROTEIN_ALTERING = {"missense", "stop_gained", "stop_lost", "start_lost", "frameshift",
                    "inframe_insertion", "inframe_deletion", "splice_acceptor", "splice_donor",
                    "splice_region"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()

    csq = "%INFO/BCSQ" if has_tag(a.vcf, "BCSQ") else "."
    clnsig = "%INFO/CLNSIG" if has_tag(a.vcf, "CLNSIG") else "."
    clndn = "%INFO/CLNDN" if has_tag(a.vcf, "CLNDN") else "."
    fmt = f"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT]\t[%DP]\t{csq}\t{clnsig}\t{clndn}\n"

    n_snv = n_indel = ti = tv = het = hom = 0
    consequences = collections.Counter()
    genes = collections.Counter()
    clinical = []
    with open(a.table, "w") as out:
        out.write("chrom\tpos\tref\talt\tqual\tgenotype\tdepth\tgene\tconsequence\tprotein_change\tclinvar_significance\tclinvar_condition\n")
        for line in query(a.vcf, fmt):
            chrom, pos, ref, alt, qual, gt, dp, bcsq, sig, dn = line.split("\t")
            alt1 = alt.split(",")[0]
            if len(ref) == 1 and len(alt1) == 1:
                n_snv += 1
                ti += (ref, alt1) in TRANSITIONS
                tv += (ref, alt1) not in TRANSITIONS
            else:
                n_indel += 1
            if gt.replace("|", "/") in ("1/1",):
                hom += 1
            elif "1" in gt:
                het += 1

            gene = cons = aa = "."
            if bcsq not in (".", ""):
                first = bcsq.split(",")[0].split("|")
                cons = first[0].lstrip("@*")
                gene = first[1] if len(first) > 1 else "."
                aa = first[5] if len(first) > 5 else "."
                for c in cons.split("&"):
                    consequences[c] += 1
                if any(c in PROTEIN_ALTERING for c in cons.split("&")):
                    genes[gene] += 1
            else:
                consequences["intergenic/non-coding (no BCSQ)"] += 1

            if sig not in (".", "") and ("athogenic" in sig) and "onflicting" not in sig:
                clinical.append((chrom, pos, ref, alt, gene, aa, sig, dn))
            out.write("\t".join([chrom, pos, ref, alt, qual, gt, dp, gene, cons, aa, sig, dn]) + "\n")

    titv = ti / tv if tv else float("nan")
    with open(a.report, "w") as r:
        r.write(f"# Variant summary: {a.sample}\n\n")
        r.write("| Metric | Value |\n|---|---|\n")
        r.write(f"| PASS SNVs | {n_snv:,} |\n| PASS indels | {n_indel:,} |\n")
        r.write(f"| Ti/Tv (SNVs) | {titv:.2f} |\n| Heterozygous | {het:,} |\n| Homozygous alt | {hom:,} |\n\n")
        r.write("## Consequence classes (bcftools csq)\n\n| Consequence | Count |\n|---|---|\n")
        for c, n in consequences.most_common(15):
            r.write(f"| {c} | {n:,} |\n")
        r.write("\n## Genes with protein-altering variants\n\n")
        r.write(", ".join(f"{g} ({n})" for g, n in genes.most_common(30)) or "None")
        r.write("\n\n## ClinVar pathogenic / likely pathogenic matches\n\n")
        if clinical:
            r.write("| Variant | Gene | Protein | ClinVar | Condition |\n|---|---|---|---|---|\n")
            for c in clinical:
                r.write(f"| {c[0]}:{c[1]} {c[2]}>{c[3]} | {c[4]} | {c[5]} | {c[6]} | {c[7].replace('_', ' ')} |\n")
        else:
            r.write("None found in the analysed region.\n")


if __name__ == "__main__":
    main()
