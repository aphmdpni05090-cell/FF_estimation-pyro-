import pysam
import pandas as pd
from collections import defaultdict

# Get input from user
input_vcf_path = input(str('Enter the input VCF file:'))
output_vcf_file=input(str('Enter the onput VCF file:'))
# Open input VCF
vcf_in = pysam.VariantFile(input_vcf_path,"r")

# Prepare output VCF writer

vcf_out = pysam.VariantFile(output_vcf_file, "wz", header=vcf_in.header)

# Step 1: Build dictionary of positions and associated (REF, ALT)
variants_by_pos = defaultdict(list)
for record in vcf_in.fetch():
    for alt in record.alts:
        variants_by_pos[(record.chrom, record.pos)].append((record.ref, alt))

# Step 2: Identify positions to filter
positions_to_filter = set()
for (chrom, pos), alleles in variants_by_pos.items():
    if len(alleles) > 1:
        # Check if at least one is an indel
        if any(len(ref) != len(alt) for ref, alt in alleles):
            positions_to_filter.add((chrom, pos))

# Save removed positions to file
df = pd.DataFrame(list(positions_to_filter), columns=["CHROM", "POS"])
df.to_csv("removed_pos_indel.txt", sep='\t', index=False)

# Step 3: Rewind and write filtered records
vcf_in = pysam.VariantFile(input_vcf_path, "r")  # Reopen for fresh iterator

for record in vcf_in.fetch():
    if (record.chrom, record.pos) not in positions_to_filter:
        vcf_out.write(record)

# Cleanup
vcf_in.close()
vcf_out.close()
