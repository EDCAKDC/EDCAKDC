rm(list = ls())

## Installing R packages
bioPackages <-c( "corrplot", "ggrepel", "stringr", "FactoMineR",
  "factoextra", "limma", "pheatmap", "edgeR", "DESeq2", "clusterProfiler",
  "org.Hs.eg.db", "GSEABase", "tidyverse", "GSVA" )


library(limma)
library(edgeR)
library(DESeq2)
library(FactoMineR)
library(factoextra)
library(clusterProfiler)
library(org.Hs.eg.db)
library(FactoMineR)
