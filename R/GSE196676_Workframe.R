library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(usethis)
library(future)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)
library(officer)
library(rvg)
library(qs)
library(hMsig)

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

#Idents(data) <- "Merge"

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

qsave(data, file = "GSE196676/GSE196676_data.qs")
data <- qread("GSE196676/GSE196676_data.qs")
data <- qread("GSE196676/GSE196676_data_MkP.qs")

#Annotation

#split
group_list <- SplitObject(data, split.by = "orig.ident")
HC_data <- group_list$HC
ITP_data <- group_list$ITP
DimPlot(HC_data, reduction = "umap", label = T)
DimPlot(ITP_data, reduction = "umap", label = T)

#subset only MkP cells
Idents(MkP) <- MkP$orig.ident


#Enricment
MkP <- subset(data, idents = c("MkP"))
Idents(MkP)<-MkP$orig.ident
MkP_markers <- FindMarkers(MkP, ident.1="ITP", ident.2="HC")
fwrite(MkP_markers, "GSE196676/ITPvsHC_MkP.csv", row.names = T)

MkP1 <- subset(data, idents = c("MkP"))
Idents(MkP1)<-MkP1$orig.ident

MkP1_markers<-FindMarkers(MkP1, ident.1="ITP", ident.2="HC")
fwrite(MkP1_markers, "GSE196676/ITPvsHC_MkP1.csv", row.names = T)

MkP2 <- subset(data, idents = c("MkP2"))
Idents(MkP2)<-MkP2$orig.ident

MkP2_markers<-FindMarkers(MkP2, ident.1="ITP", ident.2="HC")
fwrite(MkP2_markers, "GSE196676/ITPvsHC_MkP2.csv", row.names = T)

#Enrichment
  
HC_MkP <- subset(HC_data, idents = c("MkP", "MkP2"))
ITP_MkP <- subset(ITP_data, idents = c("MkP", "MkP2"))

DimPlot(HC_MkP, reduction = "umap", label = T)
DimPlot(ITP_MkP, reduction = "umap", label = T)

#extract cell data
ITP_RNA <- as.data.frame(as.matrix(GetAssayData(ITP_MkP, assay = "RNA", layer = "data")))
HC_RNA <- as.data.frame(as.matrix(GetAssayData(HC_MkP, assay = "RNA", layer = "data")))
fwrite(ITP_RNA, "GSE196676/ITP_MkP.csv", row.names = T)
fwrite(HC_RNA, "GSE196676/HC_MkP.csv", row.names = T)

#option, load count data from original data
org_data <- fread("GSE196676/GSE196676_BM_HSPCs_gene_counts.csv")
ITP_RNA_t <- read.csv("GSE196676/ITP_MkP.csv", row.names = 1)
HC_RNA <- read.csv("GSE196676/HC_MkP.csv", row.names = 1)

HC_barcode <- colnames(ITP_RNA)
ITP_barcode <- colnames(HC_RNA)
HC_data <- org_data[, c("V1", HC_barcode), with = F]
ITP_data <- org_data[, c("V1", ITP_barcode), with = F]

#monocle3
cds <- as.cell_data_set(MkP)
cds@clusters$UMAP$clusters <- MkP@active.ident
names(cds@clusters$UMAP$clusters) <- colnames(MkP)

part <- factor(rep(1, length(colnames(MkP))))
names(part) <- colnames(MkP)
cds@clusters$UMAP$partitions <- part
cds@clusters$UMAP$louvain_res <- 1

cds <- learn_graph(cds, use_partition = F)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

qsave(MkP, file = "GSE196676/GSE196676_data_MkP.qs")

#subclustering
data <- qread(file = "GSE196676/GSE196676_data_MkP.qs")
data$seurat_clusters <- seurat_clusters(data$seurat_clusters)
merge_cluster <- unique(data$seurat_clusters)
data$seurat_clusters[data$seurat_clusters %in% merge_cluster] <- "MkP"
data$seurat_clusters <- factor(data$seurat_clusters)

DimPlot(data, reduction = "pca")
ElbowPlot(data, ndims = 50)
data <- FindNeighbors(data, dims = 1:40)
data <- FindClusters(data, resolution = 1)
data <- RunUMAP(data, dims = 1:40)
DimPlot(data, reduction = "umap")
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 13)
feat <- c("PF4", "PPBP", "GYPA", "KLF1", "GATA1", "GATA2", "CD34", "SPINK2", "CDHS", "EPCAM", "AFP")
VlnPlot(data, features = feat, pt.size = 0, stack = T)

qsave(data, file = "GSE196676/GSE196676_MkP_subset.qs")

for (i in 1:7) {
  n <- 1 + 9*(i-1)
  m <- n+8
FeaturePlot(data, features = genesets$gene[-1][n:m])

ggsave(paste0("GSE196676/Mk_feature/Mk", i, ".png"), dpi = 300)
}

Mk <- subset(data, idents = "14")
DimPlot(data, reduction = "umap", label = T,label.size=6)+NoLegend()
qsave(Mk, "GSE196676/GSE196676_MK.qs")
Mk <- qread("GSE196676/GSE196676_MK.qs")

Idents(Mk)<-Mk$orig.ident
Mk_markers<-FindMarkers(Mk, ident.1="ITP", ident.2="HC")
write.csv(Mk_markers, "GSE196676/Mk_ITPvsHC.csv", row.names = T)
DimPlot(Mk, reduction = "umap")
features <- genesets
VlnPlot(Mk, features = genesets$gene[-1][41:63])

genesets_name <- msigdb_search()
genesets <- as.data.frame(msigdb_browse("Hu", genesets_name))

cell_numbers <- table(Idents(data))

clusters_id <- Idents(data)
char_clusters_id <- as.character(Idents(data))
char_clusters_id[clusters_id %in% c("12", "7", "8")] <- "1"
char_clusters_id[clusters_id %in% c("2", "9")] <- "2"
char_clusters_id[clusters_id %in% c("0", "10", "4")] <- "3"
char_clusters_id[clusters_id %in% c("13")] <- "4"
char_clusters_id[clusters_id %in% c("1")] <- "5"
char_clusters_id[clusters_id %in% c("11")] <- "6"
char_clusters_id[clusters_id %in% c("5", "6")] <- "7"
char_clusters_id[clusters_id %in% c("3")] <- "8"
char_clusters_id[clusters_id %in% c("14")] <- "Platelet"

char_clusters_id <- factor(char_clusters_id)
data_merged <- data
Idents(data_merged) <- char_clusters_id
DimPlot(data, reduction = "umap")

qsave(data, "GSE196676/GSE196676_MkP_merged.qs")


#CytoTRACE2#####


library(CytoTRACE2)
library(CytoTRACE)

splited <- SplitObject(data, split.by = "seurat_clusters")
Mk <- splited$`2`
expr_mat <- as.matrix(GetAssayData(Mk, slot = "count", assay= "RNA"))
cyto_result <- cytotrace2(expr_mat, species = "human")
ann <- as.data.frame(data$seurat_clusters)
plots <- plotData(cytotrace2_result = cyto_result, annotation = ann, expression_data = expr_mat)

for (i in seq_along(plots)) {
  ggsave(plots[[i]], file = paste0("GSE196676/CytoTrace2/GSE196676_MkP_", names(plots[i]), ".png"), dpi = 300)
}

cyto_result <- CytoTRACE(expr_mat)
plotCytoGenes(cyto_result, numOfGenes = 10, outputDir = "./GSE196676/")
#####
#GSEA
if (!require("tidyverse", character.only = TRUE) | 
    !require("cowplot", character.only = TRUE) |
    !require("tcltk", character.only = TRUE)) {
  install.packages(c("tidyverse", "cowplot", "tcltk"))
  lapply(c("tidyverse", "cowplot", "tcltk"), library, character.only = TRUE)
} else {
  lapply(c("tidyverse", "cowplot", "tcltk"), library, character.only = TRUE)
}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("enrichplot", character.only = TRUE) | 
    !require("fgsea", character.only = TRUE) | 
    !require("msigdbr", character.only = TRUE) | 
    !require("clusterProfiler", character.only = TRUE) | 
    !require("DOSE", character.only = TRUE) |
    !require("org.Hs.eg.db", character.only = TRUE)) {
  
  BiocManager::install(c("enrichplot", "fgsea", "msigdbr", "clusterProfiler", "DOSE", "org.Hs.eg.db"))
  lapply(c("enrichplot", "fgsea", "msigdbr", "clusterProfiler", "DOSE", "org.Hs.eg.db"), library, character.only = TRUE)
} else {
  lapply(c("enrichplot", "fgsea", "msigdbr", "clusterProfiler", "DOSE", "org.Hs.eg.db"), library, character.only = TRUE)
}


#data loading
df <- MkP1_markers#...load own data, It must have geneID and logFC.
names(df) <- c("geneSymbol", "logFC") #...Rename col

#personal function of extract abs/max
max_abs_pos <- function(x) {
  max_abs <- max(abs(x))
  candidates <- x[abs(x) == max_abs]
  pos_max <- candidates[candidates > 0]
  if (length(pos_max) > 0) {
    return(pos_max[1])
  } else {
    return(candidates[1])
  }
}

#Abs_max_probe
df <- aggregate(logFC ~ geneSymbol, df, max_abs_pos)

#num list
row_names <- df$geneSymbol
geneList <- list(df$logFC)
geneList <- geneList[[1]]
names(geneList) <- row_names
geneList <- sort(geneList, decreasing = TRUE)

# H, C1, C2, C4, C7, C8
T2G <- msigdbr(species = "Homo sapiens", category = "H")  %>% 
  dplyr::select(gs_name, gene_symbol)

edo2 <- GSEA(geneList,
             exponent = 0,
             minGSSize = 1,
             maxGSSize = 3000,
             pvalueCutoff  = 1,
             pAdjustMethod = "BH",
             TERM2GENE = T2G,
             nPermSimple = 10000,
             by = "fgsea",
             eps = 0)

#log -> Sig
edo2log <- edo2
pvalues <- -log10(edo2@result[["pvalue"]])
pvalues[pvalues > 6] <- 6 #10을 원하는 숫자로 변경 시 최대값 설정 가능
edo2log@result[["pvalue"]] <- pvalues
names(edo2log@result)[names(edo2log@result) == 'pvalue'] <- "Sig"

#dotplot
dotplot(edo2log, showCategory=50, font.size = 10, label_format = 60, color = "NES", size = "Sig") +
  scale_fill_gradient(low="blue", high="red", name = 'NES') +
  scale_size(name = "Sig", breaks = c(2.5, 5, 7.5, 10), limits = c(0, 10)) +
  ggtitle("Hallmark") +
  theme_test()

dotplot(edo2, showCategory=50, font.size = 10, label_format = 60, color = "NES", size = "pvalue") +
  scale_fill_gradient(low="blue", high="red", name = 'NES') +
  scale_size(name = "p-value", breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  ggtitle("Hallmark") +
  theme_test()

#GO
library(clusterProfiler)
library(org.Hs.eg.db)

gene_df <- bitr(names(geneList)[abs(geneList) > 2],
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)


gene <- names(geneList)[abs(geneList) > 2]

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                keyType = "SYMBOL",
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
goplot(ego)

edo <- enrichDGN(gene_df$ENTREZID)
barplot(edo, showCategory=20) 



mito.genes<-grep(pattern="^MT-",x=rownames(data),value=T)
