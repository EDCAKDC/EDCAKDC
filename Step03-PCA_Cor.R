rm(list = ls())
options(stringsAsFactors = FALSE)

# Load libraries
library(FactoMineR)
library(factoextra)
library(corrplot)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)

# Load preprocessed expression data
lname <- load(file = 'data/Step01-airwayData.Rdata')
print(lname)

# 1. Hierarchical clustering of samples (dendrogram)

dat <- log2(express_cpm + 1)
sampleTree <- hclust(dist(t(dat)), method = "average")

pdf("result/2.sample_Treeplot.pdf", width = 7, height = 6)
plot(sampleTree,
     main = "Cluster Dendrogram",
     sub = "", xlab = "", cex = 0.8)
dev.off()

# 2. PCA visualization (centroid version)

dat <- as.data.frame(t(log2(express_cpm + 1)))  # Samples as rows

# Match group labels
sample_names <- rownames(dat)
idx <- match(sample_names, group$run_accession)
group_list <- group[idx, 2]
valid <- !is.na(group_list)
dat <- dat[valid, ]
group_list <- sub("_[0-9]+$", "", group_list[valid])  # Clean suffix
group_list <- factor(group_list)

# Run PCA
res.pca <- PCA(dat, graph = FALSE)

# Extract PC1 & PC2 coordinates
coords <- as.data.frame(res.pca$ind$coord[, 1:2])
colnames(coords) <- c("Dim1", "Dim2")
coords$group <- group_list

# Compute group centroids
centroids <- coords %>%
  group_by(group) %>%
  summarize(Dim1 = mean(Dim1), Dim2 = mean(Dim2))

# Axis labels with % variance explained
pct1 <- round(res.pca$eig[1, 2], 1)
pct2 <- round(res.pca$eig[2, 2], 1)
xlab <- paste0("Dim1 (", pct1, "%)")
ylab <- paste0("Dim2 (", pct2, "%)")

# Colors & shapes
ngrp <- nrow(centroids)
cols <- colorRampPalette(brewer.pal(9, "Set1"))(ngrp)
shps <- c(16, 17, 15)[seq_len(ngrp)]
names(cols) <- names(shps) <- centroids$group

# Plot centroid PCA
p <- ggplot(centroids, aes(x = Dim1, y = Dim2, color = group, shape = group)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 5) +
  geom_text_repel(aes(label = group), show.legend = FALSE) +
  scale_color_manual(name = "Treatment Group", values = cols) +
  scale_shape_manual(name = "Treatment Group", values = shps) +
  labs(title = "Individuals – PCA", x = xlab, y = ylab) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 10),
        legend.position = "right")

ggsave("result/2.sample_PCA_centroids.pdf", p, width = 7.5, height = 6.5)

# 3. PCA visualization (individual points using fviz)

group_list <- group[match(rownames(dat), group$run_accession), 2]
valid_idx <- !is.na(group_list)
dat <- dat[valid_idx, ]
group_list <- factor(group_list[valid_idx])

color_count <- length(unique(group_list))
palette_colors <- colorRampPalette(brewer.pal(9, "Set1"))(color_count)

dat_pca <- PCA(dat, graph = FALSE)
p <- fviz_pca_ind(
  dat_pca,
  geom.ind = "point",
  col.ind = group_list,
  palette = palette_colors,
  addEllipses = FALSE,
  legend.title = "Groups"
) + theme_bw()

pdf(file = "result/2.sample_PCA.pdf", width = 8, height = 6)
print(p)
dev.off()

# 4. Sample-to-sample correlation heatmap

exprSet <- log2(express_cpm + 1)
exprSet <- exprSet[names(sort(apply(exprSet, 1, mad), decreasing = TRUE)[1:800]), ]
M <- cor(exprSet, method = "spearman")

anno <- data.frame(group = group$sample_title,
                   row.names = group$run_accession)

pheatmap(M,
         display_numbers = TRUE,
         annotation_col = anno,
         fontsize = 10,
         cellheight = 20,
         cellwidth = 20,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         filename = "result/2.sample_Cor.pdf",
         width = 10, height = 8)
