library(Seurat)
library(data.table)

sub1 <- ReadMtx(
  mtx = "GSE212092/GSM6509143_Subject1_Day0_BM_matrix.mtx.gz",
  features = "GSE212092/GSM6509143_Subject1_Day0_BM_features.tsv.gz",
  cells = "GSE212092/GSM6509143_Subject1_Day0_BM_barcodes.tsv.gz"
)
sub2 <- ReadMtx(
  mtx = "GSE212092/GSM6509146_Subject2_Day0_BM_matrix.mtx.gz",
  features = "GSE212092/GSM6509146_Subject2_Day0_BM_features.tsv.gz",
  cells = "GSE212092/GSM6509146_Subject2_Day0_BM_barcodes.tsv.gz"
)
sub4 <- ReadMtx(
  mtx = "GSE212092/GSM6509149_Subject4_Day0_BM_matrix.mtx.gz",
  features = "GSE212092/GSM6509149_Subject4_Day0_BM_features.tsv.gz",
  cells = "GSE212092/GSM6509149_Subject4_Day0_BM_barcodes.tsv.gz"
)

S_sub1 <- CreateSeuratObject(counts = sub1)
S_sub2 <- CreateSeuratObject(counts = sub2)
S_sub4 <- CreateSeuratObject(counts = sub4)

data <- merge(x = S_sub1,
              y = list(S_sub2, S_sub4),
              add.cell.ids = c("Sub1", "Sub2", "Sub4"))

#QC
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 30000)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)


