import argparse

parser = argparse.ArgumentParser(description='Protein variantion associations')
parser.add_argument('--GTMatrix', type=str, help='Path to the GT matrix file')
parser.add_argument('--loci', type=str, help='Reported loci')   # e.g. chr9:133261703:A:G
parser.add_argument('--protPath', type=str, help='Path to the protein data')

# %%
args = parser.parse_args()

GT_matrix_file = args.GTMatrix
protPath = args.protPath
loci = args.loci

# GT_matrix_file = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_reported_loci/result_loci_PRO_assoc/06.GT_prepare_naga/GTmatrix_df_naga.csv'
# protPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Proteome/case_control_prepare'
# loci = 'chr9:133261703:A:G'

# %%
import pandas as pd
import os

df_PRO = pd.read_csv(os.path.join(protPath, 'CTEPH_day0_vs_NAGAHAMA_Pro_log2.csv'), index_col=0)
df_PRO.index = df_PRO.index.str.replace('_day0', '')

df_SP = pd.read_csv(os.path.join(protPath, 'sample.csv'), index_col=0)
df_SP = df_SP.drop('SampleCode', axis=1)
df_SP.index = df_SP.index.str.replace('_day0', '')

df_GT = pd.read_csv(GT_matrix_file, index_col=0)
df_GT = df_GT.T

# %%
common_index = df_GT.index.intersection(df_PRO.index)

df_GT = df_GT.loc[common_index]
df_PRO = df_PRO.loc[common_index]
df_SP = df_SP.loc[common_index]

df_GT = df_GT.sort_index()
df_PRO = df_PRO.sort_index()
df_SP = df_SP.sort_index()

df_GT_loci = df_GT[loci]

# %%
genotype_mapping = {'hom-ref': 0, 'het': 1, 'hom-alt': 2, 'no-call': pd.NA}
df_GT_loci = df_GT_loci.replace(genotype_mapping)

df_merged = df_SP.join(df_GT_loci).join(df_PRO)

# %%
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import numpy as np

# 初始化结果列表
results = []

# 遍历每个蛋白质表达量
for protein in df_PRO.columns:
    # 提取自变量和因变量
    X = df_merged[['sex_male1_female2', 'age', loci]]
    y = df_merged[protein]
    
    # 删除X中loci为NA的行以及y中对应的行
    X.dropna(subset=[loci], inplace=True)
    y = y[X.index]

    # 添加截距项
    X = sm.add_constant(X)
    
    # 将X和y中的数据转换为数值
    X = X.apply(pd.to_numeric, errors='coerce')
    y = pd.to_numeric(y, errors='coerce')

    # 使用sklearn进行拟合
    model_sklearn = LinearRegression().fit(X, y)
    beta = model_sklearn.coef_[2]  # 基因型的系数

    model_sm = sm.OLS(y, X).fit()
    p_value = model_sm.pvalues[loci]

    # 计算OR
    odds_ratio = np.exp(beta)
    
    results.append([protein, beta, odds_ratio, p_value])

df_results = pd.DataFrame(results, columns=['UniProt', 'Beta', 'OR', 'p-value'])
pvals_fdr_bh = multipletests(df_results['p-value'], method='fdr_bh')[1]
df_results['FDR'] = pvals_fdr_bh

loci_name = loci.replace(':', '_')

df_results.to_csv(f'{loci_name}_PRO_result.csv', index=False) #Save the results to a CSV file


# %%
import numpy as np
import matplotlib.pyplot as plt

# 计算-log(p-value)
df_results['neg_log_pvalue'] = -np.log10(df_results['p-value'])
neg_log_pvalues = df_results['neg_log_pvalue'].values

uniform_pvals = np.linspace(1/len(neg_log_pvalues), 1, len(neg_log_pvalues))
theoretical_neg_log_pvals = -np.log10(uniform_pvals)

theoretical_neg_log_pvals_sorted = np.sort(theoretical_neg_log_pvals)
observed_neg_log_pvals_sorted = np.sort(neg_log_pvalues)

plt.style.use('default')
plt.figure(figsize=(6, 6))
plt.scatter(theoretical_neg_log_pvals_sorted, observed_neg_log_pvals_sorted, edgecolor='k', facecolor='none')

# 添加y=x线
plt.plot([np.min((theoretical_neg_log_pvals_sorted.min(), observed_neg_log_pvals_sorted.min())), 
          np.max((theoretical_neg_log_pvals_sorted.max(), observed_neg_log_pvals_sorted.max()))], 
         [np.min((theoretical_neg_log_pvals_sorted.min(), observed_neg_log_pvals_sorted.min())), 
          np.max((theoretical_neg_log_pvals_sorted.max(), observed_neg_log_pvals_sorted.max()))], 'r-')

plt.title(f'QQ Plot of {loci}')
plt.xlabel('Theoretical -log10(p-values)')
plt.ylabel('Observed -log10(p-values)')
plt.grid(True)

# Save the figure as a PDF file
plt.tight_layout()
plt.savefig(f'{loci_name}_QQ_plot.pdf', format='pdf') # Save the figure as a PDF file