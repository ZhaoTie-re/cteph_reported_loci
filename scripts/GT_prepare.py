import argparse

parser = argparse.ArgumentParser(description='VCF file to GT matrix')
parser.add_argument('--vcfPath', type=str, help='Path to the VCF file')
parser.add_argument('--GTMatrix', type=str, help='Name of the output GT matrix file')
parser.add_argument('--nocallFigure', type=str, help='Name of the output no-call proportion figure')

args = parser.parse_args()

vcfPath = args.vcfPath
GTMatrix = args.GTMatrix
nocallFigure = args.nocallFigure

# %%
import pandas as pd
import pysam

# Open the VCF file
vcf = pysam.VariantFile(vcfPath)

# Get the sample names
samples = list(vcf.header.samples)

# Create a dictionary to store the genotype data
data = {}

# Iterate over the records in the VCF file
for record in vcf:
    # Get the ID of the record
    id = record.id

    # Create a dictionary to store the genotype data for this record
    record_data = {}

    # Iterate over the samples
    for sample in samples:
        # Get the genotype data for this sample
        genotype = record.samples[sample]['GT']

        # Convert the genotype data to a string
        if genotype == (0, 0):
            genotype_str = 'hom-ref'
        elif genotype == (0, 1) or genotype == (1, 0):
            genotype_str = 'het'
        elif genotype == (1, 1):
            genotype_str = 'hom-alt'
        else:
            genotype_str = 'no-call'

        # Add the genotype data to the dictionary
        record_data[sample] = genotype_str

    # Add the record data to the main data dictionary
    data[id] = record_data

# Convert the data dictionary to a DataFrame
df = pd.DataFrame(data).transpose()
df.to_csv(GTMatrix)

reported_loci = df.index.to_list()
with open('reported_loci_ls.txt', 'w') as f:
    for item in reported_loci:
        f.write("%s\n" % item)

# %%
import matplotlib.pyplot as plt

# Calculate the proportion of 'hom', 'het' and 'no-call' for each position
df['hom_ref_proportion'] = df.apply(lambda row: (row == 'hom-ref').mean(), axis=1)
df['hom_alt_proportion'] = df.apply(lambda row: (row == 'hom-alt').mean(), axis=1)
df['het_proportion'] = df.apply(lambda row: (row == 'het').mean(), axis=1)
df['no_call_proportion'] = df.apply(lambda row: (row == 'no-call').mean(), axis=1)

# Sort the DataFrame 
df = df.sort_values('hom_ref_proportion', ascending=False)

# Plot the proportion of 'hom-ref', 'hom-alt', 'het' and 'no-call'
plt.style.use('fivethirtyeight')
plt.figure(figsize=(10, 6))
plt.bar(df.index, df['hom_ref_proportion'], label='hom-ref')
plt.bar(df.index, df['hom_alt_proportion'], bottom=df['hom_ref_proportion'], label='hom-alt')
plt.bar(df.index, df['het_proportion'], bottom=df['hom_ref_proportion']+df['hom_alt_proportion'], label='het')
plt.bar(df.index, df['no_call_proportion'], bottom=df['hom_ref_proportion']+df['hom_alt_proportion']+df['het_proportion'], label='no-call')
plt.ylabel('Proportion')
plt.xlabel('')
plt.title('Proportion of hom-ref, hom-alt, het and no-call for each SNP')
plt.xticks(rotation=60, ha='right')
plt.legend(loc='upper left', bbox_to_anchor=(1,1))

# Save the figure as a PDF file
plt.tight_layout()
plt.savefig(nocallFigure, format='pdf')