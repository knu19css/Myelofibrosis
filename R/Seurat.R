#Seurat

Data<-read.csv(filename,"",row.names=1)
Seuratobject<-CreateSeuratObject(counts=Data)