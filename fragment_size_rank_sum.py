

import pysam
import argparse
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

def compute_z_score(ref_sizes, alt_sizes):
    if not ref_sizes or not alt_sizes:
        return None

    try:
        u_stat, _ = mannwhitneyu(ref_sizes, alt_sizes, alternative='two-sided')
        n1, n2 = len(ref_sizes), len(alt_sizes)
        mean_u = n1 * n2 / 2
        std_u = ((n1 * n2) * (n1 + n2 + 1) / 12) ** 0.5
        z = (u_stat - mean_u) / std_u
        return round(z, 4)
    except Exception as e:
        print("⚠️ Error computing Z-score:", e)
        return None

def get_base_from_read(read, pos):
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos == pos and qpos is not None:
            try:
                return read.query_sequence[qpos]
            except:
                return None
    return None

def analyze_insert_sizes(vcf_path, bam_path, output_path):
    vcf = pysam.VariantFile(vcf_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    results = []

    for rec in vcf:
        chrom, pos, ref, alts = rec.chrom, rec.pos, rec.ref, rec.alts

        for alt in alts:
            ref_sizes = []
            alt_sizes = []

            fetch_end = pos + max(len(ref), len(alt))
            for read in bam.fetch(chrom, pos - 1, fetch_end):
                if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
                    continue

                base = get_base_from_read(read, pos - 1)
                if base is None:
                    continue

                insert_size = abs(read.template_length)
                if insert_size == 0:
                    continue  # avoid garbage values

                if base == ref:
                    ref_sizes.append(insert_size)
                elif base == alt or len(alt) != len(ref):  # allow indel
                    alt_sizes.append(insert_size)

            z = compute_z_score(ref_sizes, alt_sizes)
            results.append({
                "Chr": chrom,
                "Pos": pos,
                "Ref": ref,
                "Alt": alt,
                "Z_score": z if z is not None else "NA",
                "N_ref": len(ref_sizes),
                "N_alt": len(alt_sizes)
            })
            
      
    df = pd.DataFrame(results)
    df.to_csv(output_path, sep="\t", index=False)
    print(f"✅ Output written to: {output_path}")
    print("this command for cleaning null values on Bash: awk 'NR==1 || !/NA(\t|$)/' insert_size_ranksum.tsv > insert_size_ranksum.cleaned.tsv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert size Z-score analysis on consensus reads (VCF + BAM).")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--bam", required=True, help="Input BAM file (consensus reads)")
    parser.add_argument("--out", default="insert_size_rank_sum.tsv", help="Output TSV file")
    args = parser.parse_args()

    analyze_insert_sizes(args.vcf, args.bam, args.out)


