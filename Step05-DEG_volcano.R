# Step04-DEG_Volcano.R
# -----------------------------------------------------
# This script generates volcano plots from edgeR DEG results.
# It labels top up/down-regulated genes using ENSEMBL -> SYMBOL mapping.
# -----------------------------------------------------

# ———— Clear environment & load libraries ————
rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)
library(ggrepel)
library(org.Hs.eg.db)      # Use corresponding OrgDb if analyzing non-human data
library(AnnotationDbi)

# ———— Load DEG results from Step03 ————
# Each element is a data.frame with gene IDs, logFC, PValue, etc.
load("data/Step03-edgeR_nrDEG.Rdata")

# ———— Global settings ————
fc_cutoff <- 1.5      # Fold-change threshold
p_cutoff  <- 0.05     # P-value threshold
nLabel    <- 6        # Number of top genes (up/down) to label
offset    <- 0.8      # Label horizontal offset (adjust based on logFC range)

# ———— Loop through all comparisons to generate volcano plots ————
for (comp in names(deg_list)) {
  
  # 1. Extract DEG result and prepare data
  df <- deg_list[[comp]] %>%
    tibble::rownames_to_column("ENSEMBL") %>%
    as_tibble()
  
  # 2. Map ENSEMBL IDs to gene SYMBOLs
  map <- AnnotationDbi::select(org.Hs.eg.db,
                               keys    = df$ENSEMBL,
                               columns = "SYMBOL",
                               keytype = "ENSEMBL")
  df <- left_join(df, map, by = "ENSEMBL") %>%
    mutate(Gene = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))
  
  # 3. Annotate regulation status
  df <- df %>%
    mutate(regulated = case_when(
      logFC >  log2(fc_cutoff) & PValue < p_cutoff ~ "Up",
      logFC < -log2(fc_cutoff) & PValue < p_cutoff ~ "Down",
      TRUE                                          ~ "None"
    )) %>%
    mutate(regulated = factor(regulated, levels = c("Down", "None", "Up")))
  
  # 4. Select top significant genes for labeling
  up_lab   <- df %>% filter(regulated == "Up")   %>% arrange(PValue) %>% slice_head(n = nLabel)
  down_lab <- df %>% filter(regulated == "Down") %>% arrange(PValue) %>% slice_head(n = nLabel)
  
  # 5. Extract group names from comparison label
  groups <- strsplit(comp, "_vs_")[[1]]
  g1 <- groups[1]; g2 <- groups[2]
  
  # 6. Create volcano plot
  p <- ggplot(df, aes(x = logFC, y = -log10(PValue))) +
    geom_point(aes(color = regulated), size = 1.8, alpha = 0.6) +
    scale_color_manual(values = c(Down = "dodgerblue", None = "grey70", Up = "red")) +
    geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)),
               linetype = "dashed", color = "orange", size = 0.7) +
    geom_hline(yintercept = -log10(p_cutoff),
               linetype = "dashed", color = "orange", size = 0.7) +
    
    # Label down-regulated genes (left side)
    geom_text_repel(
      data = down_lab,
      aes(label = Gene),
      nudge_x = -offset,
      hjust = 1,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +
    
    # Label up-regulated genes (right side)
    geom_text_repel(
      data = up_lab,
      aes(label = Gene),
      nudge_x = offset,
      hjust = 0,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +
    
    # Theme and labels
    theme_bw(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste0("Volcano: ", g1, " vs ", g2),
      x     = paste0("log2FC (", g1, " vs ", g2, ")"),
      y     = "-log10(PValue)"
    )
  
  # 7. Save plot
  ggsave(
    filename = paste0("result/Volcano_", comp, ".png"),
    plot     = p,
    width    = 7, height = 6, dpi = 150
  )
}
