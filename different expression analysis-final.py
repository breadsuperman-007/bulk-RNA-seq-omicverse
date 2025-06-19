import omicverse as ov
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns

#————————————————————Step1:读取文件——————————————————————
data=ov.utils.read('/Users/zhanghuitian/Desktop/GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv',index_col=0)
print(data.head())
#————————————————————Step2:下载基因ID与基因名的对应关系——————————————————————
#通过 omicverse下载基因ID与基因名的对应关系（第一次用的时候下载即可，这后都不用再下载）
# ov.utils.download_geneid_annotation_pair()
#————————————————————Step3:将基因ID转化为基因名——————————————————————
#调取小鼠基因组版本：GRCm39，与IDmapping
data=ov.bulk.Matrix_ID_mapping(data,'genesets/pair_GRCh38.tsv')
print(data.head())
#————————————————————Step4:将data转换成DESeq对象——————————————————————
dds=ov.bulk.pyDEG(data)
#去除重复的index
dds.drop_duplicates_index()
print('... drop_duplicates_index success')
#————————————————————Step5:进行差异表达分析——————————————————————
#下载GSE213001_series_matrix.txt.gz
#整理成一个excel，包含样品信息和分组信息
file_path = '/Users/zhanghuitian/Desktop/groupname.xlsx'
groupname = pd.read_excel(file_path)
print(groupname.head())

sample_titles = groupname["Sample_title"]
source_names = groupname["Sample_source_name_ch1"]

is_ipf = source_names.str.contains("IPF", case=False, na=False)
is_ndc = source_names.str.contains("NDC", case=False, na=False)

ipf_samples = sample_titles[is_ipf].reset_index(drop=True)
ndc_samples = sample_titles[is_ndc].reset_index(drop=True)

treatment_groups = ipf_samples.tolist()
control_groups = ndc_samples.tolist()

print("实验组（IPF）Sample_title：")
print(treatment_groups)

print("\n对照组（NDC）Sample_title：")
print(control_groups)

#使用Deseq2进行差异分析，得到dataframe包含FC，p值，基因名等信息
result=dds.deg_analysis(treatment_groups,control_groups,method='DEseq2')
print(result.head())
print(result.shape)

#只保留表达量不低的基因（BaseMean > 2），去掉可能的低表达噪声或测不准的基因
result=result.loc[result['log2(BaseMean)']>1]
print(result.shape)
print(result.columns)
# ————————————————————Step6:设置FC和p值后绘图——————————————————————
# ----------------------- 数据处理 -----------------------
# 计算 -log10(pvalue)
result['-log10(pvalue)'] = -np.log10(result['pvalue'])

# 设定阈值
fc_threshold = 1.5
pval_threshold = 0.05

# 根据阈值标记显著性
result['significant'] = 'no'
result.loc[(result['log2FoldChange'] >= fc_threshold) & (result['pvalue'] <= pval_threshold), 'significant'] = 'up'
result.loc[(result['log2FoldChange'] <= -fc_threshold) & (result['pvalue'] <= pval_threshold), 'significant'] = 'down'
#显示上调和下调的基因数量
n_up = (result['significant'] == 'up').sum()
n_down = (result['significant'] == 'down').sum()
print(f"显著上调基因数：{n_up} 个")
print(f"显著下调基因数：{n_down} 个")

# 选择fold change最大和最小的5个基因
top5_up = result[result['significant'] == 'up'].sort_values('log2FoldChange', ascending=False).head(5).index.tolist()
top5_down = result[result['significant'] == 'down'].sort_values('log2FoldChange').head(5).index.tolist()
highlight_genes = top5_up + top5_down

# ----------------------- 绘图 -----------------------
plt.figure(figsize=(7, 6))

# 散点图
sns.scatterplot(
    data=result,
    x='log2FoldChange',
    y='-log10(pvalue)',
    hue='significant',
    palette={'up': 'red', 'down': 'blue', 'no': 'gray'},
    alpha=0.6,
    edgecolor=None,
    legend=False
)

# 添加横纵坐标阈值线
plt.axhline(-np.log10(pval_threshold), color='black', linestyle='--', linewidth=0.8)
plt.axvline(fc_threshold, color='black', linestyle='--', linewidth=0.8)
plt.axvline(-fc_threshold, color='black', linestyle='--', linewidth=0.8)

# 标注基因，添加连线，避免重叠
texts = []
for gene in highlight_genes:
    x = result.loc[gene, 'log2FoldChange']
    y = result.loc[gene, '-log10(pvalue)']
    plt.scatter(x, y, color='orange', s=50)  # 重点标记基因点
    texts.append(plt.text(x, y, gene, fontsize=10))
    plt.plot([x, x], [y, y + 1.5], color='gray', linestyle='--', linewidth=0.8)

adjust_text(texts, only_move={'points':'y', 'texts':'y'}, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

# 设置纵坐标从0开始
plt.ylim(0, None)

plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot')

plt.show()


sns.set_palette(['blue', 'red'])
dds.plot_boxplot(
    genes=['ATP1A1', 'COL1A2'],
    treatment_groups=treatment_groups,
    control_groups=control_groups,
    figsize=(4, 4),
    fontsize=12,
    legend_bbox=(1, 1)  #图例位置
)
plt.show()

#————————————————————Step7:富集分析——————————————————————

deg_genes = dds.result.loc[dds.result['sig'] != 'normal'].index.tolist()

pathway_dict_wiki = ov.utils.geneset_prepare('/Users/zhanghuitian/Downloads/WikiPathway_2021_Human.txt', organism='human')

enr_wiki = ov.bulk.geneset_enrichment(
    gene_list=deg_genes,
    pathways_dict=pathway_dict_wiki,
    pvalue_type='auto',     # 自动选择合适的p值方法
    organism='human'
)
print(enr_wiki.head())

top_n = 20
df_plot = enr_wiki.sort_values('logp', ascending=False).head(top_n)

plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=df_plot,
    x='logp',              # 横轴：富集显著性（-log10 P 值）
    y='Term',              # 纵轴：通路名
    size='num',            # 点大小：显著基因数
    hue='logp',            # 点颜色：富集显著性
    palette='Reds',
    sizes=(40, 400)
)

plt.xlabel('-log10(P-value)')
plt.ylabel('WikiPathway Term')
plt.title('WikiPathway Enrichment Analysis')
plt.tight_layout()
plt.legend(title='Gene Count', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()




