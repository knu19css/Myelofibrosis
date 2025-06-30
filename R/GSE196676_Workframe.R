library(Seurat)
library(clustree)
library(usethis)
library(future)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)
library(officer)
library(rvg)

plan(sequential)
plan()
options(future.globals.maxSize = 30*1024^3)
plan(multisession, workers = 8)
plan()

#load

data <- fread("GSE196676/GSE196676_BM_HSPCs_gene_counts.csv")
data <- as.matrix(as.data.frame(data[, -1], row.names = data$V1))
data <- as(data, "dgCMatrix")
data <- CreateSeuratObject(count = data)

Idents(data) <- "Merge"

data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 30000)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

data <- NormalizeData(data)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 20)
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

rm(plot1, plot2)

all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)

data <- RunPCA(data, features = VariableFeatures(object = data))
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")

DimPlot(data, reduction = "pca") + NoLegend()
DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data, dims = 1:20, cells = 500, balanced = TRUE)

ElbowPlot(data)

data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 1.5)
head(Idents(data), 5)

data <- RunUMAP(data, dims = 1:20)
DimPlot(data, reduction = "umap", label = T)
VlnPlot(data, features = top10, stack = T)
