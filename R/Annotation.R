scTypedb <- fread("R/GSE196676_markergenes.csv")
scTypedb_ic <- subset(scTypedb, scTypedb$tissueType == "Immune system")

celltype <- scTypedb_ic$cellName

is <- setNames(
  lapply(strsplit(scTypedb_ic$geneSymbolmore1, ","), unlist),
  scTypedb_ic$cellName)

names(is) <- gsub("/", "_", names(is))

p <- DimPlot(data, reduction = "umap", label = T)
ggsave(p, filename = paste0("GSE196676/features2/0_Umap.png"), dpi = 300)

for (i in 1:length(is)) {
  temp1 <- VlnPlot(data, features = is[[i]], stack = T)
  ggsave(temp1, filename = paste0("GSE196676/features2/", i, "_VlnPlot_", names(is)[i], ".png"), dpi = 300)
  
  # expr <- FetchData(data, vars = is[[i]])
  # avg_expr <- colMeans(expr)
  # top9_genes <- names(sort(avg_expr, decreasing = TRUE))[1:9]
  temp2 <- FeaturePlot(data, features = is[[i]])
  ggsave(temp2, filename = paste0("GSE196676/features2/", i, "_Ft_", names(is)[i], ".png"), dpi = 300)
}

library(officer)
image_list <- list.files("GSE196676/features3", full.names = T)
image_name <- list.files("GSE196676/features3")

ppt <- read_pptx()

for (i in seq_along(image_list)) {
  txt <- fpar(
    ftext(image_name[i], fp_text(font.size = 10, color = "black"))
  )
  
  ppt %>% add_slide(layout = "Blank", master = "Office Theme") %>% 
    ph_with(
      external_img(image_list[i], width = 6, height = 4),
      location = ph_location_fullsize()
    ) %>%
    ph_with(
      value = txt,
      location = ph_location(left = 0, top = 0, width = 6, height = 1)
    )
}

print(ppt, target = "ppt/GSE196676_feature_auto2.pptx")


markergenes <- unlist(is)
DoHeatmap(data@meta.data$seurat_clusters[1:10], features = markergenes)

for (i in seq_along(is)) {
p <- DotPlot(data, is[[i]]) + 
  ggtitle(toupper(names(is)[i])) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("GSE196676/dotplot/", i, "_", names(is)[i], ".png"), plot = p, dpi = 300)
}

library(officer)
image_list <- list.files("GSE196676/dotplot", full.names = T)
image_name <- list.files("GSE196676/dotplot")

image_list <- image_list[order(as.numeric(gsub("\\D", "", basename(image_list))))]
image_name <- image_name[order(as.numeric(gsub("\\D", "", basename(image_name))))]

ppt <- read_pptx()

for (i in seq_along(image_list)) {
  txt <- fpar(
    ftext(image_name[i], fp_text(font.size = 10, color = "black"))
  )
  
  ppt %>% add_slide(layout = "Blank", master = "Office Theme") %>% 
    ph_with(
      external_img(image_list[i], width = 6, height = 4),
      location = ph_location_fullsize()
    ) %>%
    ph_with(
      value = txt,
      location = ph_location(left = 0, top = 0, width = 6, height = 1)
    )
}

print(ppt, target = "ppt/GSE196676_dotplot.pptx")

cluster_name <- fread("GSE196676/feature.csv")

new.cluster.ids <- cluster_name$V2
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)

DimPlot(data, reduction = "umap", label = T)
