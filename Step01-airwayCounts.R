rm(list = ls())
library(stringr)

## ==================== 1. Load Data

# Read raw count expression matrix
rawcount <- read.table("data/gene_counts.xls", row.names = 1, 
                       sep = "\t", header = TRUE)

# Select columns of interest (modify x and y based on your dataset)
rawcount <- rawcount[, x:y]  # x and y are the column indices to keep
colnames(rawcount)

# Preview expression matrix
rawcount[1:6, 1:6]

# Dimensions before filtering
dim(rawcount)

# Read sample group information (treated vs untreated)
group <- data.table::fread("data/group.txt", data.table = FALSE)
group

## ==================== 2. Filter and Normalize Data

# Filter genes: keep those expressed (count > 0) in at least 75% of samples
library(edgeR)
keep <- rowSums(rawcount > 0) >= floor(0.75 * ncol(rawcount))

# Check how many genes are kept
table(keep)

# Apply the filter to the expression matrix
filter_count <- rawcount[keep, ]
filter_count[1:4, 1:4]

# Show filtered matrix dimensions
dim(filter_count)

# Convert raw counts to CPM (Counts Per Million)
express_cpm <- cpm(filter_count)
express_cpm[1:6, 1:6]

# Save filtered count matrix, CPM matrix, and group info for downstream analysis
save(filter_count, express_cpm, group, file = "data/Step01-airwayData.Rdata")
