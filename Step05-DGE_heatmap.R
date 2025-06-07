# Step05-DEG_Heatmap.R
# ----------------------------------------------------------
# For each pairwise comparison in sig_list, this script:
# 1. Selects top N up/down regulated genes by logFC & P-value
# 2. Extracts CPM expression matrix
# 3. Maps ENSEMBL IDs to SYMBOLs
# 4. Plots a heatmap of the selected genes
# 5. Saves corresponding gene lists (up/down separately)
# ----------------------------------------------------------

# Clear environment & load libraries
rm(list = ls())
options(stringsAsFactors = FALSE)

library(dplyr)
library(pheatmap)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)  # For ENSEMBL to SYMBOL mapping

# Load normalized CPM expression data and group info
load("data/Step01-airwayData.Rdata")
rownames(express_cpm) <- sub("\\..*$", "", rownames(express_cpm))  # Remove version suffix from ENSEMBL

# Load significant DEGs from previous step
load("data/Step03-edgeR_nrDEG.Rdata")  # Contains sig_list

# Output directory
outdir <- "result"
if (!dir.exists(outdir)) dir.create(outdir)

# Set number of top up/down genes to plot
topN <- 50

# Loop over all pairwise comparisons
for (comp in names(sig_list)) {
  message("=== Processing ", comp, " ===")
  
  # 1. Retrieve DEG table and calculate ranking score
  df <- sig_list[[comp]]
  df$ENSEMBL <- rownames(df)
  df <- df %>%
    mutate(Score = abs(logFC) * -log10(PValue)) %>%
    arrange(desc(Score))
  
  # 2. Select top N up/down regulated genes
  top_df <- df %>%
    filter(regulated != "normal") %>%
    group_by(regulated) %>%
    slice_max(order_by = Score, n = topN) %>%
    ungroup()
  
  if (nrow(top_df) == 0) {
    warning("No DEGs found for ", comp)
    next
  }
  
  # 3. Get ENSEMBL IDs (remove version suffix if needed)
  ens_keep <- unique(sub("\\..*$", "", top_df$ENSEMBL))
  
  # 4. Identify sample groups based on contrast
  parts <- strsplit(comp, "_vs_")[[1]]  # e.g., c("RD", "RT")
  samp1 <- parts[1]
  samp2 <- parts[2]
  
  # 5. Get sample names for each group
  samples1 <- grep(paste0("^", samp1), colnames(express_cpm), value = TRUE)
  samples2 <- grep(paste0("^", samp2), colnames(express_cpm), value = TRUE)
  samples  <- c(samples1, samples2)
  
  if (length(samples) < 2) {
    warning("No valid samples found for comparison: ", comp)
    next
  }
  
  # 6. Subset CPM expression matrix for selected genes & samples
  mat <- express_cpm[ens_keep, samples, drop = FALSE]
  if (nrow(mat) == 0) {
    warning("No expression data found for selected genes: ", comp)
    next
  }
  
  # 7. Map ENSEMBL IDs to gene SYMBOLs using EnsDb
  gene_ids <- rownames(mat)
  sym_map <- mapIds(
    EnsDb.Hsapiens.v86,
    keys = gene_ids,
    column = "SYMBOL",
    keytype = "GENEID",
    multiVals = "first"
  )
  keep_idx <- which(!is.na(sym_map) & sym_map != "")
  mat <- mat[keep_idx, , drop = FALSE]
  rownames(mat) <- sym_map[keep_idx]
  
  # 8. Save up/down gene lists (ENSEMBL + logFC + PValue)
  write.csv(
    top_df %>% filter(regulated == "up") %>% select(ENSEMBL, logFC, PValue),
    file = file.path(outdir, paste0("UpGenes_", comp, ".csv")),
    row.names = FALSE
  )
  write.csv(
    top_df %>% filter(regulated == "down") %>% select(ENSEMBL, logFC, PValue),
    file = file.path(outdir, paste0("DownGenes_", comp, ".csv")),
    row.names = FALSE
  )
  
  # 9. Prepare sample annotation for heatmap (group label = prefix)
  anno <- data.frame(
    Group = ifelse(grepl(paste0("^", samp1), samples), samp1, samp2),
    row.names = samples
  )
  
  # 10. Draw and save heatmap
  png(file.path(outdir, paste0("Heatmap_", comp, ".png")),
      width = 1500, height = 1800, res = 150)
  
  pheatmap(
    mat,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = anno,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    fontsize_row = 7,
    fontsize_col = 10,
    main = paste("Top DEGs Heatmap -", comp)
  )
  
  dev.off()
}

message("✅ All pairwise heatmaps & gene lists saved in '", outdir, "/'.")
