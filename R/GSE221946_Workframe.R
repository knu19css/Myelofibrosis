library(Seurat)
library(future)
library(Matrix)
library(data.table)
library(dplyr)
library(glmGamPoi)
library(ggplot2)
library(officer)
library(rvg)
getwd()

plan(sequential)
options(future.globals.maxSize = 20e9)
plan(multisession)
plan()

#load
data <- fread("GSE221946/GSM6910347_S1_RNA_counts.csv")
data <- as.matrix(as.data.frame(data[, -1], row.names = data$V1))
data <- as(data, "dgCMatrix")
data <- CreateSeuratObject(count = data)


data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- SCTransform(data, vars.to.regress = "percent.mt", verbose = FALSE)

data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
rm(plot1, plot2)

data <- RunPCA(data, verbose = FALSE)
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca") + NoLegend()
DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(data)

data <- FindNeighbors(data, dims = 1:17)
data <- FindClusters(data, resolution = 0.7)
head(Idents(data), 5)

data <- RunUMAP(data, dims = 1:17, verbose = FALSE)
DimPlot(data, reduction = "umap", label = T)

markers <- FindAllMarkers(data, only.pos = F, verbose = T)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


#markers by SS
#####
#megakaryocytes
VlnPlot(data, features = c("GP9","ARTNL","ITGA2B","PLEK","G6B","ITGB3"))
FeaturePlot(data, features = c("GP9","ARTNL","ITGA2B","PLEK","G6B","ITGB3"))

#erythrocyte
VlnPlot(data, features = c("HBG2", "HBB"))
FeaturePlot(data, features = c("HBG2", "HBB"))

#Platelets
VlnPlot(data, features = c("PPBP"))
FeaturePlot(data, features = c("PPBP"))

#eosinophils
VlnPlot(data, features = c("IL3RA", "SIGLEC"))
FeaturePlot(data, features = c("IL3RA", "SIGLEC"))

#basophils
VlnPlot(data, features = c("IL3RA", "ENPP3", "IL5RA"))
FeaturePlot(data, features = c("IL3RA", "ENPP3", "IL5RA"))

#Neutrophils
VlnPlot(data, features = c("CXCR2", "FCGR3B", "CSF3R", "G0S2"))
FeaturePlot(data, features = c("CXCR2", "FCGR3B", "CSF3R", "G0S2"))

#Monocytes
VlnPlot(data, features = c("CD14", "LYZ", "VCAN", "FCN1", "FCGR3A"))
FeaturePlot(data, features = c("CD14", "LYZ", "VCAN", "FCN1", "FCGR3A"))

#Dendritic cells
VlnPlot(data, features = c("FCER1A", "CD1C", "NRP1", "CLEC4C", "IL3RA"))
FeaturePlot(data, features = c("FCER1A", "CD1C", "NRP1", "CLEC4C", "IL3RA"))

#B cell
VlnPlot(data, features = c("MS4A1", "IGHM", "CD79A"))
FeaturePlot(data, features = c("MS4A1", "IGHM", "CD79A"))

#T cell
VlnPlot(data, features = c("CD3E"))
FeaturePlot(data, features = c("CD3E"))

#CD4 T cell
VlnPlot(data, features = c("CD4", "ILR7", "CCR7"))
FeaturePlot(data, features = c("CD4", "ILR7", "CCR7"))

#CD8 T cell
VlnPlot(data, features = c("CD8A"))
FeaturePlot(data, features = c("CD8A"))
#####

#markers by Azimuth
#####
#HSPC
VlnPlot(data, features = unlist(strsplit("ARPP21,CFAP73,AVP,SPTA1,EPCAM,RHCE,AKAP12,RHAG,C1QTNF4,CYTL1", ",")))
FeaturePlot(data, features = unlist(strsplit("ARPP21,CFAP73,AVP,SPTA1,EPCAM,RHCE,AKAP12,RHAG,C1QTNF4,CYTL1", ",")))

#Mono
VlnPlot(data, features = unlist(strsplit("LYPD2,FOLR3,CLEC4E,LILRA1,CDA,RBP7,CD300LF,FPR1,CD93,MTMR11", ",")))
FeaturePlot(data, features = unlist(strsplit("LYPD2,FOLR3,CLEC4E,LILRA1,CDA,RBP7,CD300LF,FPR1,CD93,MTMR11", ",")))

#B
VlnPlot(data, features = unlist(strsplit("FCRL1,FCRL2,CD22,ARHGAP24,BANK1,MS4A1,RALGPS2,CD37,SWAP70,CD79A",",")))
FeaturePlot(data, features = unlist(strsplit("FCRL1,FCRL2,CD22,ARHGAP24,BANK1,MS4A1,RALGPS2,CD37,SWAP70,CD79A",",")))

#CD4 T
VlnPlot(data, features = unlist(strsplit("MAL,LEF1,IL7R,TCF7,LTB,LEPROTL1,CD3E,CD3D,CD3G,CD27", ",")))
FeaturePlot(data, features = unlist(strsplit("MAL,LEF1,IL7R,TCF7,LTB,LEPROTL1,CD3E,CD3D,CD3G,CD27", ",")))

#CD8 T
VlnPlot(data, features = unlist(strsplit("CD8B,CD8A,NELL2,CD3D,CD3E,S100B,GZMH,CD3G,TRGC2,CCL5", ",")))
FeaturePlot(data, features = unlist(strsplit("CD8B,CD8A,NELL2,CD3D,CD3E,S100B,GZMH,CD3G,TRGC2,CCL5", ",")))

#DC
VlnPlot(data, features = unlist(strsplit("CLEC4C,PROC,SCT,SCN9A,SHD,PPM1J,ENHO,CLEC10A,LILRA4,DNASE1L3", ",")))
FeaturePlot(data, features = unlist(strsplit("CLEC4C,PROC,SCT,SCN9A,SHD,PPM1J,ENHO,CLEC10A,LILRA4,DNASE1L3", ",")))

#NK
VlnPlot(data, features = unlist(strsplit("NCR1,SH2D1B,KLRC1,PRSS23,CD160,B3GNT7,IL2RB,IL18RAP,XCL1,KLRF1",",")))
FeaturePlot(data, features = unlist(strsplit("NCR1,SH2D1B,KLRC1,PRSS23,CD160,B3GNT7,IL2RB,IL18RAP,XCL1,KLRF1", ",")))

#Other T
VlnPlot(data, features = unlist(strsplit("SLC4A10,KLRB1,IL7R,TNFRSF25,RORA,SPOCK2,MAF,GZMK,ZFP36L2,LYAR", ",")))
FeaturePlot(data, features = unlist(strsplit("SLC4A10,KLRB1,IL7R,TNFRSF25,RORA,SPOCK2,MAF,GZMK,ZFP36L2,LYAR", ",")))
#####

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(data, features = top10$gene) + NoLegend() +
  theme(axis.text.y = element_text(size = 8))

#celldex, singleR auto-annotaion, #not available result... I'll try again next time with the other data
library(celldex)
library(SingleR)
DefaultAssay(data) <- "RNA"
test.data <- as.SingleCellExperiment(data)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
result <- SingleR(test = test.data, ref = hpca.se, assay.type.test = 1, labels = hpca.se$label.fine)

DimPlot(data, reduction = "umap", label = T, repel = T)
tab <- table(cluster = data$seurat_clusters, label = result$labels)
p <- pheatmap::pheatmap(t(log10(tab+10)), cluster_cols = F)

#Azimuth
#remotes::install_github('satijalab/azimuth', ref = 'master')
library(Azimuth)
#not work yet, need to build-up

#ppt
ppt <- read_pptx()
ppt %>% add_slide(layout = "Blank", master = "Office Theme") %>% 
  ph_with(
  dml(ggobj = p),
  location = ph_location_fullsize()
)
print(ppt, target = "ppt/GSE221946_S1.pptx")

ggsave("Feature.png", dpi = 600)
