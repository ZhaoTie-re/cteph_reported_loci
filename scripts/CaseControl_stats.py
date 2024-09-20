import argparse

parser = argparse.ArgumentParser(description='Statistics pertain to the case-control design VCF file')
parser.add_argument('--case_vcfPath', type=str, help='Path to the case VCF files')
parser.add_argument('--control_vcfPath', type=str, help='Path to the control VCF files')
parser.add_argument('--lociKnown', type=str, help='Path to the reported loci file (excel)')

args = parser.parse_args()
# %%
lociKnown_path = args.lociKnown
cteph_vcf_file = args.case_vcfPath
naga_vcf_file = args.control_vcfPath

# %%
import pandas as pd

lociKnown = pd.read_excel(lociKnown_path)

# %%
import pysam

cteph_vcf = pysam.VariantFile(cteph_vcf_file)

# 创建一个新的DataFrame来存储结果
cteph_result = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'cteph total ref reads', 'cteph total alt reads', 'cteph total reads', 'cteph allele frequency'])

# 遍历lociKnown的每一行
for index, row in lociKnown.iterrows():
    chrom = row['CHROM']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']

    total_ref_reads = 0
    total_alt_reads = 0
    total_reads = 0

    # 在vcf文件中查找匹配的记录
    for rec in cteph_vcf.fetch(chrom, pos-1, pos):
        if rec.ref == ref and alt in rec.alts:
            # 遍历每个样本
            for sample in rec.samples:
                # 累加ref reads和alt reads
                total_ref_reads += rec.samples[sample]['AD'][0]
                total_alt_reads += rec.samples[sample]['AD'][1]
                total_reads += rec.samples[sample]['DP']

    # 计算allele frequency
    allele_frequency = total_alt_reads / total_reads if total_reads > 0 else 0

    # 将结果添加到result DataFrame
    cteph_result.loc[len(cteph_result)] = [chrom, pos, ref, alt, total_ref_reads, total_alt_reads, total_reads, allele_frequency]

# 输出结果
print(cteph_result)

# %%
import pysam

cteph_vcf = pysam.VariantFile(cteph_vcf_file)

# 创建一个新的DataFrame来存储结果
tommo_result = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'TOMMO_AF'])

# 遍历lociKnown的每一行
for index, row in lociKnown.iterrows():
    chrom = row['CHROM']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']

    # 在vcf文件中查找匹配的记录
    for rec in cteph_vcf.fetch(chrom, pos-1, pos):
        if rec.ref == ref and alt in rec.alts:
            # 获取TOMMO_AF值
            tommo_af = rec.info.get('TOMMO_AF')[0]

            # 将结果添加到result DataFrame
            tommo_result.loc[len(tommo_result)] = [chrom, pos, ref, alt, tommo_af]

# 输出结果
print(tommo_result)



# %%
import pysam

naga_vcf = pysam.VariantFile(naga_vcf_file)

# 创建一个新的DataFrame来存储结果
naga_result = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'naga total ref reads', 'naga total alt reads', 'naga total reads', 'naga allele frequency'])

# 遍历lociKnown的每一行
for index, row in lociKnown.iterrows():
    chrom = row['CHROM']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']

    total_ref_reads = 0
    total_alt_reads = 0
    total_reads = 0

    # 在vcf文件中查找匹配的记录
    for rec in naga_vcf.fetch(chrom, pos-1, pos):
        if rec.ref == ref and alt in rec.alts:
            # 遍历每个样本
            for sample in rec.samples:
                # 累加ref reads和alt reads
                total_ref_reads += rec.samples[sample]['AD'][0]
                total_alt_reads += rec.samples[sample]['AD'][1]
                total_reads += rec.samples[sample]['DP']

    # 计算allele frequency
    allele_frequency = total_alt_reads / total_reads if total_reads > 0 else 0

    # 将结果添加到result DataFrame
    naga_result.loc[len(naga_result)] = [chrom, pos, ref, alt, total_ref_reads, total_alt_reads, total_reads, allele_frequency]

# 输出结果
print(naga_result)

# %%
merged_result = pd.merge(naga_result, cteph_result, on=['CHROM', 'POS', 'REF', 'ALT'], how='outer')
merged_result = pd.merge(merged_result, tommo_result, on=['CHROM', 'POS', 'REF', 'ALT'], how='outer')

# %%
from scipy.stats import chi2_contingency, fisher_exact

# Initialize new columns
merged_result['Odds Ratio'] = 0
merged_result['Test Used'] = ""
merged_result['p-value'] = 0

# Define function to decide which test to use and calculate p-value
def calculate_stats(row):
    contingency_table = [
        [row['cteph total alt reads'], row['cteph total ref reads']],
        [row['naga total alt reads'], row['naga total ref reads']]
    ]
    
    # Calculate the Odds Ratio
    if row['cteph total ref reads'] == 0 or row['naga total ref reads'] == 0:
        or_value = float('inf')  # avoid division by zero
    else:
        or_value = (row['cteph total alt reads'] * row['naga total ref reads']) / \
                   (row['cteph total ref reads'] * row['naga total alt reads'])
    
    # Calculate expected counts and decide on which test to use
    chi2, p_value, _, expected_counts = chi2_contingency(contingency_table)
    if (expected_counts < 5).any():
        # Use Fisher's exact test if any expected count is less than 5
        test_used = "Fisher"
        _, p_value = fisher_exact(contingency_table)
    else:
        # Use Chi-square test
        test_used = "Chi-square"
    
    return pd.Series([or_value, test_used, p_value])

# Apply the function
merged_result[['Odds Ratio', 'Test Used', 'p-value']] = merged_result.apply(calculate_stats, axis=1)
merged_result.head()


# %%
merged_result = pd.merge(merged_result, lociKnown, on=['CHROM', 'POS', 'REF', 'ALT'], how='outer')

# %%
merged_result.columns

# %%
result = merged_result[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'GENE', 'RS', 'COUNTRY', 'cteph allele frequency', 'naga allele frequency', 'TOMMO_AF', 'Odds Ratio', 'p-value', 'Test Used']]
result.to_csv('case_control_stats.csv', index=False)
