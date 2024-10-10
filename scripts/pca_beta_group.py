import argparse
import json

parser = argparse.ArgumentParser(description='PCA of beta values')  
parser.add_argument('--group', type=str, help='Group name')
parser.add_argument('--assocSummary', type=str, help='List of mutation and association summary csv path')
parser.add_argument('--info_path', type=str, help='Path to the loci info xlsx file')

args = parser.parse_args()

# %%
import pandas as pd
import os

group = args.group
assocSummary = json.loads(args.assocSummary)
loci_info_path = args.info_path

# %%
loci_info = pd.read_excel(loci_info_path, sheet_name='Sheet1')

# %%
result_df = pd.DataFrame()

for mutation, csv_path in assocSummary:
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        if 'UniProt' in df.columns and 'Beta' in df.columns:
            beta_values = df.set_index('UniProt')['Beta']
            result_df = pd.concat([result_df, pd.DataFrame([beta_values], index=[mutation])])

result_df.index.name = 'ID'
result_df.to_csv(f'{group}_beta_values.csv')
# result_df.fillna(value=pd.NA, inplace=True)


# %%
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import plotly.io as pio


pca = PCA()
pca_result = pca.fit_transform(result_df)
pca_df = pd.DataFrame(pca_result)
pca_df.index = result_df.index

assert pca_df.shape[1] >= 3, "PCA result has less than three principal components."


pca_df_3d = pca_df.iloc[:, :3]
pca_df_3d.columns = ['PC1', 'PC2', 'PC3']


explained_variance_ratio = pca.explained_variance_ratio_[:3]
total_explained_variance = explained_variance_ratio.sum() * 100


hover_text = loci_info.set_index('ID').loc[pca_df_3d.index, ['GENE', 'RS', 'COUNTRY']].apply(lambda x: 'ID:' + x.name + '<br>' + '<br>'.join([f'{col}:{val}' for col, val in zip(x.index, x.astype(str))]), axis=1)
point_text = loci_info.set_index('ID').loc[pca_df_3d.index, 'GENE']


fig = go.Figure(data=[go.Scatter3d(
    x=pca_df_3d['PC1'],
    y=pca_df_3d['PC2'],
    z=pca_df_3d['PC3'],
    mode='markers+text',
    text=point_text,
    hovertext=hover_text,
    hoverinfo='text',
    marker=dict(
        size=3,
        colorscale='Viridis',
        opacity=0.8
    )
)])


fig.update_layout(
    title=f"PCA of {group} - Total explained variance: {total_explained_variance:.2f}%",
    scene=dict(
        xaxis_title=f"PC1 ({explained_variance_ratio[0]*100:.2f}%)",
        yaxis_title=f"PC2 ({explained_variance_ratio[1]*100:.2f}%)",
        zaxis_title=f"PC3 ({explained_variance_ratio[2]*100:.2f}%)"
    ),
    width=600,
    height=600
)

pio.write_html(fig, f"{group}_pca_beta.html")
