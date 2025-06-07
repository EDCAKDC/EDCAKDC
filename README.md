# 👋 Hi there, I'm EDCAKDC

🎓 I'm currently learning bioinformatics, and I'd like to share the RNA-seq analysis code I've been studying.

---

## 🔬 Featured Project: RNA-seq DEG Analysis

A step-by-step R pipeline for analyzing differential gene expression using the **airway dataset**, visualizing PCA, volcano plots, heatmaps, and performing GO/KEGG/GSEA enrichment.

📋 Sample Group Metadata (`group.txt`)
This file contains metadata used to define sample groups for differential expression analysis.

Each row corresponds to one sample, including its unique ID (run_accession) and its experimental label (sample_title).

There are 3 conditions:

A: untreated control group

B: treatment group 

C: the other treatment group 
Each group has 3 biological replicates.

Example group.txt:

run_accession    sample_title
A1              XXX_1
A2              XXX_2
A3              XXX_3
B1              YYY_1
B2              YYY_2
B3              YYY_3
C1              ZZZ_1
C2              ZZZ_2
C3              ZZZ_3

### 📁 Project Structure

```text
airway-DEG-analysis/
├── Step00-R packages.R           <- Required R packages
├── Step01-airwayCounts.R         <- Load expression matrix
├── Step02-sampleDistribution.R   <- Plot sample distribution
├── Step03-PCA_Cor.R              <- Perform PCA and sample correlation analysis
├── Step04-edgeR_DEG.R            <- Differential expression analysis with edgeR
├── Step05-DEG_volcano.R          <- Generate volcano plot of DEGs
├── Step05-DGE_heatmap.R          <- Draw heatmap of top DEGs
├── Step06-GO_KEGG_enrich.R       <- GO and KEGG enrichment analysis
├── Step06-GSEA_analysis.R        <- GSEA enrichment analysis
├── data/                         <- Input data (e.g. count matrix, metadata)
├── result/                       <- Output figures and result tables
└── Diff_analysis/                <- Intermediate differential analysis outputs
