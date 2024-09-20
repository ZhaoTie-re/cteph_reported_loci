import argparse

parser = argparse.ArgumentParser(description='Extract records from VCF files')
parser.add_argument('--vcfPath', type=str, help='Path to the VCF files')
parser.add_argument('--lociKnown', type=str, help='Path to the reported loci file (excel)')
parser.add_argument('--outputFile', type=str, help='Name of the output file')

args = parser.parse_args()
# %%
cteph_vcf_path = args.vcfPath
reported_loci_path = args.lociKnown
output_file = args.outputFile
# %%
import pandas as pd
import pysam
import os

reported_loci = pd.read_excel(reported_loci_path)
chrom_pos_ref_alt = reported_loci[['CHROM', 'POS', 'REF', 'ALT']]


# 初始化一个空列表来存储找到的记录
found_records = []

# 处理 DataFrame 中的每一行
for idx, row in reported_loci.iterrows():
    # 构造 VCF 文件的文件名
    vcf_filename = os.path.join(cteph_vcf_path, f"{row['CHROM'].replace('chr', '')}.tommo.vcf.gz")
    
    # 打开 VCF 文件
    vcf_reader = pysam.VariantFile(vcf_filename)
    
    # 搜索与匹配位置、参考等位基因和替代等位基的记录
    for record in vcf_reader.fetch(row['CHROM'], row['POS']-1, row['POS']):
        if record.ref == row['REF'] and record.alts[0] == row['ALT']:
            # 如果找到，将记录添加到列表中
            found_records.append(record)
            break

found_records = sorted(found_records, key=lambda record: (int(record.chrom[3:]), record.pos))

# 将找到的记录保存到新的 VCF 文件中
with pysam.VariantFile(output_file, 'w', header=vcf_reader.header) as new_vcf:
    for record in found_records:
        new_vcf.write(record)

pysam.tabix_index(output_file, preset='vcf')
