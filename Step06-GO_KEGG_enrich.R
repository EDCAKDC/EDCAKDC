# Clear environment & load packages
rm(list=ls())
options(stringsAsFactors=FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(ggplot2)
library(tidyverse)
library(patchwork)

# If KEGG download fails, uncomment below (for Windows)
# options(clusterProfiler.download.method = "wininet")

# Load DEG results
load("data/Step03-edgeR_nrDEG.Rdata")
comps <- names(sig_list)

# Load Hallmark gene sets (optional for GSEA)
geneset <- read.gmt("data/MsigDB/v7.4/h.all.v7.4.symbols.gmt")
geneset$term <- str_to_title(sub("HALLMARK_", "", geneset$term))

# Output directory
if(!dir.exists("result")) dir.create("result")

# Loop through each contrast to generate 4-panel enrichment plot
for(comp in comps) {
  message(">>> Drawing combined plot for ", comp)
  
  # Map ENSEMBL -> SYMBOL
  df <- sig_list[[comp]]
  df$ENSEMBL <- rownames(df)
  id2s <- bitr(df$ENSEMBL, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)
  df2 <- merge(df, id2s, by="ENSEMBL", all.x=TRUE)
  
  # Filter significant SYMBOL genes
  genes <- df2 %>%
    filter(regulated != "normal", !is.na(SYMBOL)) %>%
    pull(SYMBOL) %>%
    unique()
  
  if(length(genes) == 0) {
    message("⚠️ No valid SYMBOL genes for ", comp, ". Consider relaxing p-value or fold-change thresholds.")
    next
  }
  
  # GO enrichment
  ego_BP <- enrichGO(gene=genes, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="BP", pvalueCutoff=1, qvalueCutoff=1)
  ego_MF <- enrichGO(gene=genes, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="MF", pvalueCutoff=1, qvalueCutoff=1)
  ego_CC <- enrichGO(gene=genes, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="CC", pvalueCutoff=1, qvalueCutoff=1)
  
  # KEGG enrichment
  ekid <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>%
    pull(ENTREZID) %>% unique()
  ekegg <- enrichKEGG(gene=ekid, organism="hsa", keyType="ncbi-geneid", pvalueCutoff=1, qvalueCutoff=1)
  
  # Individual plots
  p1 <- barplot(ego_BP, showCategory=10, label_format=function(x) x) +
    ggtitle("Biological Process") + theme(axis.text.y=element_text(size=10), plot.title=element_text(hjust=0.5))
  
  p2 <- barplot(ego_MF, showCategory=10, label_format=function(x) x) +
    ggtitle("Molecular Function") + theme(axis.text.y=element_text(size=10), plot.title=element_text(hjust=0.5))
  
  p3 <- barplot(ego_CC, showCategory=10, label_format=function(x) x) +
    ggtitle("Cellular Component") + theme(axis.text.y=element_text(size=10), plot.title=element_text(hjust=0.5))
  
  p4 <- barplot(ekegg, showCategory=10, label_format=function(x) x) +
    ggtitle("KEGG Pathway") + theme(axis.text.y=element_text(size=10), plot.title=element_text(hjust=0.5))
  
  # Merge into 2x2 layout
  title_text <- paste0("Comparison: ", sub("_vs_", " vs ", comp))
  combo <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = title_text,
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
    )
  
  ggsave(
    filename = paste0("result/Combined_", comp, ".png"),
    plot     = combo,
    width    = 12,
    height   = 10,
    dpi      = 150
  )
}
# GO BP Enrichment: upregulated vs downregulated genes (separate)
for(comp in names(sig_list)) {
  message(">>> Processing GO-BP for ", comp)
  df <- sig_list[[comp]]
  df$ENSEMBL <- rownames(df)
  
  id2s <- bitr(df$ENSEMBL, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db) %>% distinct(ENSEMBL, SYMBOL)
  df2 <- merge(df, id2s, by="ENSEMBL", all.x=TRUE)
  
  genes_up <- df2 %>% filter(regulated=="up", !is.na(SYMBOL)) %>% pull(SYMBOL) %>% unique()
  genes_down <- df2 %>% filter(regulated=="down", !is.na(SYMBOL)) %>% pull(SYMBOL) %>% unique()
  
  ego_up <- if(length(genes_up)>0) enrichGO(genes_up, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="BP",
                                            pvalueCutoff=0.05, qvalueCutoff=0.2) else NULL
  ego_dn <- if(length(genes_down)>0) enrichGO(genes_down, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="BP",
                                              pvalueCutoff=0.05, qvalueCutoff=0.2) else NULL
  
  n_up <- if(!is.null(ego_up)) min(10, nrow(ego_up)) else 0
  n_dn <- if(!is.null(ego_dn)) min(10, nrow(ego_dn)) else 0
  
  p_up <- if(n_up > 0) {
    barplot(ego_up, showCategory=n_up) +
      labs(title="GO – BP (Upregulated)") +
      theme(axis.text.y=element_text(size=10), plot.title=element_text(hjust=0.5))
  } else {
    ggplot() + ggtitle("No Upregulated GO-BP Terms") + theme_void()
  }
  
  p_dn <- if(n_dn > 0) {
    barplot(ego_dn, showCategory=n_dn) +
      labs(title="GO – BP (Downregulated)") +
      theme(axis.text.y=element_text(size=10), plot.title=element_text(hjust=0.5))
  } else {
    ggplot() + ggtitle("No Downregulated GO-BP Terms") + theme_void()
  }
  
  combo <- p_up + p_dn +
    plot_annotation(title = paste0("GO – BP: ", sub("_vs_"," vs ", comp)),
                    theme = theme(plot.title=element_text(size=16, face="bold", hjust=0.5)))
  
  ggsave(file.path("result", paste0("GO_", comp, "_combined.png")),
         combo, width=12, height=6, dpi=150)
}
# KEGG Enrichment by Up/Down-regulated genes
for(comp in comps) {
  message(">>> KEGG up/down for ", comp)
  
  df <- sig_list[[comp]]
  df$ENSEMBL <- rownames(df)
  id2s <- bitr(df$ENSEMBL, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db) %>% distinct(ENSEMBL, SYMBOL)
  df2 <- merge(df, id2s, by="ENSEMBL", all.x=TRUE)
  
  genes_up <- df2 %>% filter(regulated=="up", !is.na(SYMBOL)) %>% pull(SYMBOL) %>% unique()
  genes_down <- df2 %>% filter(regulated=="down", !is.na(SYMBOL)) %>% pull(SYMBOL) %>% unique()
  
  ids_up <- if(length(genes_up)>0) bitr(genes_up, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID else character()
  ids_dn <- if(length(genes_down)>0) bitr(genes_down, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID else character()
  
  ek_up <- if(length(ids_up)>0) enrichKEGG(ids_up, organism="hsa", keyType="ncbi-geneid", pvalueCutoff=1, qvalueCutoff=1) else NULL
  ek_dn <- if(length(ids_dn)>0) enrichKEGG(ids_dn, organism="hsa", keyType="ncbi-geneid", pvalueCutoff=1, qvalueCutoff=1) else NULL
  
  p_up <- if(!is.null(ek_up) && nrow(ek_up)>0) {
    barplot(ek_up, showCategory=10) +
      labs(title="KEGG – Upregulated") +
      theme(axis.text.y=element_text(size=10), plot.title=element_text(hjust=0.5))
  } else {
    ggplot() + ggtitle("No Upregulated KEGG Terms") + theme_void()
  }
  
  p_dn <- if(!is.null(ek_dn) && nrow(ek_dn)>0) {
    barplot(ek_dn, showCategory=10) +
      labs(title="KEGG – Downregulated") +
      theme(axis.text.y=element_text(size=10), plot.title=element_text(hjust=0.5))
  } else {
    ggplot() + ggtitle("No Downregulated KEGG Terms") + theme_void()
  }
  
  combo <- (p_up | p_dn) +
    plot_annotation(title = paste0("KEGG Pathway Enrichment: ", sub("_vs_"," vs ", comp)),
                    theme = theme(plot.title=element_text(size=18, face="bold", hjust=0.5)))
  
  ggsave(filename = file.path("result", paste0("KEGG_", comp, "_up_down.png")),
         plot = combo, width = 12, height = 6, dpi = 150)
}
