rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required packages
library(edgeR)
library(ggplot2)

# Load preprocessed count matrix and group info
lname <- load(file = "data/Step01-airwayData.Rdata")
print(lname)

# Preview expression matrix
filter_count[1:4, 1:4]

# Define group labels manually (3 samples per group) A: untreated. B and C: 2 different treated groups. 
group_list <- c(rep("A", 3), rep("B", 3), rep("C", 3))
group_list <- factor(group_list, levels = c("A", "C", "C"))
table(group_list)

# Construct design matrix (no intercept)
design <- model.matrix(~0 + group_list)
rownames(design) <- colnames(filter_count)
colnames(design) <- levels(group_list)
design

# Create DGEList object
DEG <- DGEList(counts = filter_count, group = group_list)

# Normalize library sizes
DEG <- calcNormFactors(DEG)

# Estimate dispersion parameters
DEG <- estimateGLMCommonDisp(DEG, design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

# Fit the GLM model
fit <- glmFit(DEG, design)

# Load limma for contrast matrix creation
library(limma)

# Define pairwise contrasts
cont.matrix <- makeContrasts(
  B_vs_A = B - A,
  C_vs_A = C - A,
  C_vs_B = C - B,
  levels = design
)

# Perform likelihood ratio tests for each contrast
lrt_B_vs_A <- glmLRT(fit, contrast = cont.matrix[,"B_vs_A"])
lrt_C_vs_A <- glmLRT(fit, contrast = cont.matrix[,"C_vs_A"])
lrt_C_vs_B <- glmLRT(fit, contrast = cont.matrix[,"C_vs_B"])

# Store LRT results in a named list
lrt_list <- list(
  C_vs_RD = lrt_C_vs_RD,
  RV_vs_RD = lrt_RV_vs_RD,
  RV_vs_C = lrt_RV_vs_C
)

# Set significance thresholds
fc_cutoff <- 1.5
pvalue <- 0.05

# Prepare lists to store DEG tables
deg_list <- list()
sig_list <- list()

# Loop through each contrast to extract and save results
for (name in names(lrt_list)) {
  # Extract full differential expression table
  deg <- as.data.frame(topTags(lrt_list[[name]], n = Inf, adjust.method = "BH")$table)
  
  # Annotate regulation status
  deg$regulated <- "normal"
  deg$regulated[deg$logFC >  log2(fc_cutoff) & deg$PValue < pvalue] <- "up"
  deg$regulated[deg$logFC < -log2(fc_cutoff) & deg$PValue < pvalue] <- "down"
  
  # Save full DEG results
  write.csv(deg,
            file = paste0("result/DEG_", name, "_all.csv"),
            row.names = FALSE)
  
  # Save significantly differentially expressed genes
  sig <- subset(deg, regulated != "normal")
  write.csv(sig,
            file = paste0("result/DEG_", name, "_sig.csv"),
            row.names = FALSE)
  
  # Store in lists
  deg_list[[name]] <- deg
  sig_list[[name]] <- sig
}

# Save all DEG result lists for downstream analysis
save(
  deg_list,
  sig_list,
  file = "data/Step03-edgeR_nrDEG.Rdata"
)
