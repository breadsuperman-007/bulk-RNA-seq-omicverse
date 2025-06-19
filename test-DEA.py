import omicverse as ov
import matplotlib.pyplot as plt
import seaborn as sns

ov.utils.ov_plot_set()
#————————————————————Step1:读取文件——————————————————————
data=ov.utils.read('/Users/zhanghuitian/Desktop/data_counts.txt',index_col=0,header=1)

#去掉 .bam 后缀
data.columns=[i.split('/')[-1].replace('.bam','') for i in data.columns]
print(data.head())
#————————————————————Step2:下载基因ID与基因名的对应关系——————————————————————
#通过 omicverse下载基因ID与基因名的对应关系
# ov.utils.download_geneid_annotation_pair()
#————————————————————Step3:将基因ID转化为基因名——————————————————————
#调取小鼠基因组版本：GRCm39，与IDmapping
data=ov.bulk.Matrix_ID_mapping(data,'genesets/pair_GRCm39.tsv')
print(data.head())
#————————————————————Step4:将data转换成DESeq对象——————————————————————
dds=ov.bulk.pyDEG(data)
#去除重复的index
dds.drop_duplicates_index()
print('... drop_duplicates_index success')
#————————————————————Step5:进行差异表达分析——————————————————————
treatment_groups=['4-3','4-4']
control_groups=['1--1','1--2']
result=dds.deg_analysis(treatment_groups,control_groups,method='DEseq2')
print(result.head())

print(result.shape)
result=result.loc[result['log2(BaseMean)']>1]
print(result.shape)
#————————————————————Step6:设置FC和p值后绘图——————————————————————
#设置 Fold Change 阈值 -1 means automatically calculates
dds.foldchange_set(fc_threshold=-1,
                   pval_threshold=0.05,
                   logp_max=10)

dds.plot_volcano(title='DEG Analysis',figsize=(4,4),
                 plot_genes_num=8,plot_genes_fontsize=12,)
# plt.show()

 # dds.plot_boxplot(genes=['Ckap2','Lef1'],treatment_groups=treatment_groups,
 #                control_groups=control_groups,figsize=(2,3),fontsize=12,
 #                 legend_bbox=(2,0.55))
# plt.show()


# Step 1: 提取差异表达基因
deg_genes = dds.result.loc[dds.result['sig'] != 'normal'].index.tolist()

pathway_dict_wiki = ov.utils.geneset_prepare('/Users/zhanghuitian/Downloads/WikiPathways_2019_Mouse.txt', organism='mouse')

enr_wiki = ov.bulk.geneset_enrichment(
    gene_list=deg_genes,
    pathways_dict=pathway_dict_wiki,
    pvalue_type='auto',     # 自动选择合适的p值方法
    organism='mouse'
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

# # Step 2: 通路富集
# # KEGG 通路富集
# pathway_dict_kegg = ov.utils.geneset_prepare('/Users/zhanghuitian/Downloads/KEGG_2019_Mouse.txt', organism='Mouse')
#
# # GO 富集
# pathway_dict_go = ov.utils.geneset_prepare('/Users/zhanghuitian/Downloads/GO_Biological_Process_2021.txt', organism='Mouse')
# pathway_dict_go_mf = ov.utils.geneset_prepare('/Users/zhanghuitian/Downloads/GO_Molecular_Function_2021.txt', organism='Mouse')
# pathway_dict_go_cc = ov.utils.geneset_prepare('/Users/zhanghuitian/Downloads/GO_Cellular_Component_2021.txt', organism='Mouse')
#
# # Step 3: KEGG 富集分析
# enr_kegg = ov.bulk.geneset_enrichment(gene_list=deg_genes,
#                                       pathways_dict=pathway_dict_kegg,
#                                       pvalue_type='auto',
#                                       organism='mouse')
# print(enr_kegg.head())
# print(enr_kegg.columns)
# # # Step 4: GO 富集分析
# # enr_go = ov.bulk.geneset_enrichment(gene_list=deg_genes,
# #                                     pathways_dict=pathway_dict_go,
# #                                     pvalue_type='auto',
# #                                     organism='mouse')
# #
# # enr_go_mf = ov.bulk.geneset_enrichment(gene_list=deg_genes,
# #                                        pathways_dict=pathway_dict_go_mf,
# #                                        pvalue_type='auto',
# #                                        organism='mouse')
# #
# # enr_go_cc = ov.bulk.geneset_enrichment(gene_list=deg_genes,
# #                                        pathways_dict=pathway_dict_go_cc,
# #                                        pvalue_type='auto',
# #                                        organism='mouse')
#
# top_n = 20
# df_plot = enr_kegg.sort_values('logp', ascending=False).head(top_n)
#
# # 开始绘图
# plt.figure(figsize=(8, 6))
# bubble = sns.scatterplot(
#     data=df_plot,
#     x='logp',              # 横坐标：-log10(p值)
#     y='Term',              # 纵坐标：通路名
#     size='num',            # 点的大小：显著基因数
#     hue='logp',            # 点的颜色：富集显著性
#     palette='Reds',
#     sizes=(40, 400),       # 点的大小范围
#     legend='brief'
# )
#
# # 美化
# plt.xlabel('-log10(P-value)')
# plt.ylabel('KEGG Pathway')
# plt.title('KEGG Pathway Enrichment')
# plt.tight_layout()
# plt.legend(title='Gene count', bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.show()
#
#
#
# # # Step 5: 绘制富集分析结果
# # # KEGG 通路富集绘图
# # ov.bulk.geneset_plot(enr_kegg, figsize=(5, 10), fig_title='KEGG Pathway Enrichment', cmap='Reds')
# #
# # # GO 通路富集绘图 - Biological Process
# # ov.bulk.geneset_plot(enr_go, figsize=(5, 10), fig_title='GO Biological Process Enrichment', cmap='Blues')
# #
# # # GO 通路富集绘图 - Molecular Function
# # ov.bulk.geneset_plot(enr_go_mf, figsize=(5, 10), fig_title='GO Molecular Function Enrichment', cmap='Blues')
# #
# # # GO 通路富集绘图 - Cellular Component
# # ov.bulk.geneset_plot(enr_go_cc, figsize=(5, 10), fig_title='GO Cellular Component Enrichment', cmap='Blues')
# #
# # plt.show()
