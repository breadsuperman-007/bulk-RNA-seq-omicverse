# bulk-RNA-seq-omicverse
🫁 Bulk RNA-seq Analysis of Pulmonary Fibrosis Using OmicVerse
基于 OmicVerse 的肺纤维化组织转录组分析项目


📌 项目简介 | Project Overview

本项目旨在使用 [OmicVerse](https://omicverse.readthedocs.io) 工具，对公开数据集 **GSE213001** 中的肺纤维化（IPF）与正常对照组织（NDC）进行差异表达分析（DEG）与通路富集分析，探索与肺纤维化相关的关键调控因子。

This repository conducts a bulk RNA-seq analysis of pulmonary fibrosis patient samples using the OmicVerse framework. The pipeline includes DEG identification, visualization (volcano plots, boxplots), and pathway enrichment.

📌 项目流程 | Project Flow

差异表达分析主要流程（见 scripts/deg_analysis.py）
读取计数矩阵（raw count）
转换 Entrez ID 为基因名
构建 pyDEG 对象并执行 DESeq2 分析
设置 FC 和 p-value 阈值绘制火山图
可视化特定基因表达（boxplot）
富集分析（基于显著差异基因）
