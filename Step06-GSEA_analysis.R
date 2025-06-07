# —— Clear workspace & load packages ——
rm(list=ls())
options(stringsAsFactors=FALSE)

library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)

# Optional: fix KEGG download on Windows
options(clusterProfiler.download.method = "wininet")

# 1. Load differential expression results (assuming it contains deg_list: DEG tables for all comparisons)
load("data/Step03-edgeR_nrDEG.Rdata")
comps <- names(deg_list)  # e.g. "MD_vs_MT", "MD_vs_MV", "MT_vs_MV"
print(comps)

# 2. Load Hallmark gene sets
geneset <- read.gmt("data/MsigDB/h.all.v2023.2.Hs.symbols.gmt")

# 3. Ensure output directory exists
if(!dir.exists("result")) dir.create("result")

# 4. Loop through comparisons and run GSEA
for(comp in comps) {
  message(">>> Running GSEA for ", comp)
  df <- deg_list[[comp]]
  
  # If no SYMBOL column, perform ENSEMBL -> SYMBOL mapping
  if(!"SYMBOL" %in% colnames(df)) {
    id2s <- bitr(rownames(df),
                 fromType="ENSEMBL",
                 toType="SYMBOL",
                 OrgDb=org.Hs.eg.db)
    df <- merge(id2s, df, by.x="ENSEMBL", by.y="row.names", all.y=TRUE)
  }
  
  # Key step: remove NA and duplicated SYMBOLs
  df <- df %>%
    filter(!is.na(SYMBOL)) %>%
    distinct(SYMBOL, .keep_all=TRUE)
  
  # Construct ranked gene list
  geneList <- df$logFC
  names(geneList) <- df$SYMBOL
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Run GSEA
  egmt <- GSEA(geneList,
               TERM2GENE   = geneset,
               pvalueCutoff= 1,
               verbose     = FALSE)
  
  # Save GSEA result table
  res <- as.data.frame(egmt@result)
  write.csv(res,
            file      = paste0("result/GSEA_", comp, ".csv"),
            row.names = FALSE)
  
  # Plot and save dotplot
  p <- dotplot(egmt, label_format=100) +
    ggtitle(paste0("GSEA: ", comp))
  ggsave(paste0("result/GSEA_dot_", comp, ".png"),
         plot   = p,
         width  = 8,
         height = 6,
         dpi    = 150)
}

message("✅ Pairwise GSEA completed. CSVs and plots are saved in the 'result/' folder.")
