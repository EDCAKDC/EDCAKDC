rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required packages and set default ggplot2 theme
library(ggplot2)
library(ggsci)
library(tidyverse)
library(RColorBrewer)

mythe <- theme_bw() + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

# Load preprocessed expression data (filtered and CPM-normalized)
lname <- load(file = "data/Step01-airwayData.Rdata")
lname  # Show names of loaded objects

# Apply log2 transformation to CPM data
exprSet <- log2(as.matrix(express_cpm) + 1)
exprSet[1:6, 1:6]

## 1. Global Expression Distribution Across Samples - Boxplot

# Reshape data for ggplot (long format: sample, expression)
data <- exprSet %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "sample", values_to = "expression")

head(data)

# Create boxplot
p <- ggplot(data = data, aes(x = sample, y = expression, fill = sample))
p1 <- p +
  geom_boxplot() +
  mythe +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab(NULL) + ylab("log2(CPM+1)") +
  scale_fill_lancet()

p1

# Save boxplot
png(file = "result/1.sample_boxplot.png", width = 800, height = 900, res = 150)
print(p1)
dev.off()

## 2. Global Expression Distribution Across Samples - Violin Plot

p2 <- p +
  geom_violin() +
  mythe +
  theme(
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 90)
  ) +
  xlab(NULL) + ylab("log2(CPM+1)") +
  scale_fill_lancet()

p2

# Save violin plot
png(file = "result/1.sample_violin.png", width = 800, height = 900, res = 150)
print(p2)
dev.off()

## 3. Global Expression Distribution Across Samples - Density Plot

m <- ggplot(data = data, aes(x = expression))
p3 <- m +
  geom_density(aes(fill = sample, colour = sample), alpha = 0.1) +
  xlab("log2(CPM+1)") +
  mythe +
  scale_fill_npg()

p3

# Save density plot
png(file = "result/1.sample_density.png", width = 800, height = 700, res = 150)
print(p3)
dev.off()
