import argparse

parser = argparse.ArgumentParser(description='Protein variantion associations with annotation')
parser.add_argument('--anno_file', type=str, help='Path to the annotation file (annot_protein.csv)')
parser.add_argument('--reported_loci', type=str, help='Path to the reported loci file (cteph_reported_loci.xlsx)')
parser.add_argument('--Pro_ass_file', type=str, help='Path to the protein association result file (chr9_127824409_G_A_PRO_result.csv)')
parser.add_argument('--loci', type=str, help='Reported loci (e.g. chr9:127824409:G:A)')

args = parser.parse_args()

anno_file_path = args.anno_file
reported_loci_path = args.reported_loci
Pro_file_path = args.Pro_ass_file
loci = args.loci

print(args)

# %%
import pandas as pd

anno_file = pd.read_csv(anno_file_path, index_col=0)
reported_loci = pd.read_excel(reported_loci_path)
Pro_file = pd.read_csv(Pro_file_path)

loci_save = loci.replace(':', '_')
gene = reported_loci[reported_loci['ID'] == loci]['GENE'].values[0]
rs_id = reported_loci[reported_loci['ID'] == loci]['RS'].values[0]
country = reported_loci[reported_loci['ID'] == loci]['COUNTRY'].values[0]

# %%
merged_file = Pro_file.merge(anno_file, on='UniProt', how='left')
merged_file.to_csv(f"{loci_save}_PRO_result_anno.csv", index=False) #保存文件

# %%
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# 计算-log(p-value)
merged_file['neg_log_pvalue'] = -np.log10(merged_file['p-value'])
neg_log_pvalues = merged_file['neg_log_pvalue'].values

uniform_pvals = np.linspace(1/len(neg_log_pvalues), 1, len(neg_log_pvalues))
theoretical_neg_log_pvals = -np.log10(uniform_pvals)

theoretical_neg_log_pvals_sorted = np.sort(theoretical_neg_log_pvals)
observed_neg_log_pvals_sorted = np.sort(neg_log_pvalues)

entrezgenesymbol_values = merged_file['EntrezGeneSymbol'].values
beta_values = merged_file['Beta'].values
fdr_values = merged_file['FDR'].values

# 获取排序的索引
sort_index = np.argsort(neg_log_pvalues)

# 使用索引对其他数组进行排序
entrezgenesymbol_values_sorted = entrezgenesymbol_values[sort_index]
beta_values_sorted = beta_values[sort_index]
fdr_values_sorted = fdr_values[sort_index]

lower_quartile = np.percentile(beta_values_sorted, 25)
upper_quartile = np.percentile(beta_values_sorted, 75)

cmap = plt.get_cmap('RdBu_r')
norm = plt.Normalize(vmin=lower_quartile, vmax=-lower_quartile)
colors = cmap(norm(beta_values_sorted))

hex_colors = [mcolors.rgb2hex(color) for color in colors]

# 创建散点图
trace0 = go.Scatter(
    x = theoretical_neg_log_pvals_sorted,
    y = observed_neg_log_pvals_sorted,
    mode = 'markers',
    marker = dict(color = hex_colors, size = 7),
    name = 'Observed vs Theoretical',
    text = entrezgenesymbol_values_sorted,  # 添加注释
    customdata = np.stack((beta_values_sorted, fdr_values_sorted), axis=-1),  # 添加OR和FDR值
    hovertemplate = 'EntrezGeneSymbol: %{text}<br>Beta: %{customdata[0]}<br>FDR: %{customdata[1]}<extra></extra>'  # 自定义悬停信息，显示注释、OR和FDR值
)

# 创建y=x线
trace1 = go.Scatter(
    x = [np.min((theoretical_neg_log_pvals_sorted.min(), observed_neg_log_pvals_sorted.min())), 
         np.max((theoretical_neg_log_pvals_sorted.max(), observed_neg_log_pvals_sorted.max()))],
    y = [np.min((theoretical_neg_log_pvals_sorted.min(), observed_neg_log_pvals_sorted.min())), 
         np.max((theoretical_neg_log_pvals_sorted.max(), observed_neg_log_pvals_sorted.max()))],
    mode = 'lines',
    line = dict(color = 'red'),
    name = 'y=x'
)

trace2 = go.Scatter(
    x = [0, -np.log10(0.05), -np.log10(0.05), 0, 0],
    y = [0, 0, -np.log10(0.05), -np.log10(0.05), 0],
    mode = 'lines',
    fill = 'toself',
    fillcolor = 'rgba(128, 128, 128, 0.2)',  # 半透明灰色
    line = dict(width = 0),  # 不显示边线
    showlegend = False  # 不显示图例
)

data = [trace2, trace0, trace1]

# 查找entrezgenesymbol_values_sorted中等于gene的位置
index = np.where(entrezgenesymbol_values_sorted == gene)[0]
if len(index) > 0:
    index = index[0]
    # 创建一个只包含一个点的散点图
    trace3 = go.Scatter(
        x = [theoretical_neg_log_pvals_sorted[index]],
        y = [observed_neg_log_pvals_sorted[index]],
        mode = 'markers',
        marker = dict(color = 'red', size = 15, symbol = 'star'),
        name = gene,
        text = [entrezgenesymbol_values_sorted[index]],  # 添加注释
        customdata = np.array([[beta_values_sorted[index], fdr_values_sorted[index]]]),  # 添加OR和FDR值
        hovertemplate = 'EntrezGeneSymbol: %{text}<br>Beta: %{customdata[0]}<br>FDR: %{customdata[1]}<extra></extra>'  # 自定义悬停信息，显示注释、OR和FDR值
    )
    data.append(trace3)

    # 添加箭头
    layout_annotations = [
        go.layout.Annotation(
            x = theoretical_neg_log_pvals_sorted[index],
            y = observed_neg_log_pvals_sorted[index],
            text = gene,
            showarrow = True,
            arrowhead = 1,
            ax = -50,
            ay = -50
        )
    ]

layout = go.Layout(
    title = go.layout.Title(
        text = f'QQ Plot of {loci} ({rs_id}) in <b><i>{gene}</i></b>',
        x = 0.5,  # 将标题的x位置设置为0.5
        xanchor = 'center'  # 将标题的x锚点设置为'center'
    ),
    xaxis = dict(title = 'Theoretical -log10(p-values)'),
    yaxis = dict(title = 'Observed -log10(p-values)', scaleanchor="x", scaleratio=1),
    showlegend = False,
    autosize=False,
    width=600,
    height=600
)
fig = go.Figure(data=data, layout=layout)

pio.write_html(fig, f"{loci_save}_PRO_result_qq_plot.html") # 保存为html文件