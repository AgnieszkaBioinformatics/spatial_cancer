```{r}
library(Seurat)
library(SeuratDisk)
library(data.table)
library(ggplot2)
```

```{r, setup, include=FALSE}
knitr::opts_knit$set(root.dir = '/mnt/Data/Projects/EyeAtlas/atlas_diseases/')
```

```{r}
setwd('/mnt/Data/Projects/EyeAtlas/atlas_diseases/count_matrix/')
datasets <- list.dirs(full.names = FALSE, recursive = FALSE)
```

```{r}
read_dataset <- function(index){
  name <- datasets[index]
  print(name)
  matrix_path <- paste0("count_matrix/", name, "/matrix_", name, ".mtx")
  features_path <- paste0("count_matrix/", name, "/genes_", name, ".tsv")
  cells_path <- paste0("count_matrix/", name, "/barcodes_", name, ".tsv")
  metadata_path <- paste0("metadata/", name, "_metadata.csv")
  
  count_matrix <- ReadMtx(mtx = matrix_path,
                          features = features_path,
                          cells = cells_path,
                          feature.column = 1)
  metadata <- fread(metadata_path)
  metadata <- tibble::column_to_rownames(metadata, "V1")
  seurat <- CreateSeuratObject(counts=count_matrix,
                               meta.data = metadata,
                               min.cells = 10,
                               min.features = 100)
  return(seurat)
}


list_of_seurat_objects<- lapply(1:length(datasets), read_dataset)


integrated <- list_of_seurat_objects[[1]]
for (i in 2:length(list_of_seurat_objects)){
  integrated <- merge(x=integrated, y=list_of_seurat_objects[[i]])
}


integrated@meta.data[["Cell.Annotation"]]  <- integrated@meta.data[["Cell.Anatomy"]]
integrated@meta.data <- subset(integrated@meta.data, select= -Cell.Anatomy)


```

Kanały po których był puszczany Scrublet. Scrublet był robiony w jupyterlabie (brak możliwości instalacji pythonowej biblioteki bezpośrednio na serwerze)

```{r}
zliczenie_sampli <- table(integrated@meta.data$GSM)
zliacznie_sampli <- as.data.frame(zliczenie_sampli)
setwd('/mnt/Data/Projects/EyeAtlas/atlas_diseases/scrublet')
write.csv(zliczenie_sampli, "liczebność_sampli.csv")
```

```{r}
scrublet_results <- fread("scrublet/scrublet.tsv")
rownames(scrublet_results) <- scrublet_results$CellIDs
integrated <- AddMetaData(integrated, scrublet_results[,c(2,3)])

thresh_doublets<-0.2
isDoublet<-rep("FALSE",nrow(integrated@meta.data))
isDoublet[integrated@meta.data$ScrubletScores>=thresh_doublets]<-T
isDoublet<-factor(isDoublet,levels=c(F,T))
integrated = AddMetaData(integrated,isDoublet,col.name="isDoublet")

plot(density(log10(integrated@meta.data$ScrubletScores)),main="Scublet's Doublet Scores");abline(v=log10(thresh_doublets))

#ggplot(integrated@meta.data, aes(x = ScrubletScores, fill = isDoublet)) +
#    geom_histogram(binwidth = 0.01, color = "grey", fill = "grey") +
#    theme_minimal() +
#    geom_vline(xintercept = thresh_doublets) +
#    facet_grid(~GSM) #TO DO, z kodu Marcina


#VlnPlot(object = integrated,
#          features = "ScrubletScores",
#          pt.size = 0.1,
#          group.by = "GSE")
#po czym grupować? Czym było Workflow we wzorcowym kodzie? Podstawiłam na GSE
```

```{r}
indeksy_dublety <- which(integrated@meta.data[["isDoublet"]]==TRUE)
barkody_dublety <- rownames(integrated@meta.data)[indeksy_dublety]
integrated <- subset(integrated, cells = barkody_dublety, invert = TRUE)
```

```{r}
integrated <- NormalizeData(integrated)
integrated <- FindVariableFeatures(integrated, selection.method="vst", nfeatures=5000)
top10 <- head(VariableFeatures(integrated), 10)
plot1 <- VariableFeaturePlot(integrated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

```

```{r}
all.genes <- rownames(integrated)
integrated <- ScaleData(integrated, features = all.genes)
integrated <- RunPCA(integrated,  npcs=200, features = VariableFeatures(object = integrated))
```

```{r}
ElbowPlot(integrated)
```

```{r}
options(future.globals.maxSize = 8000 * 1024^2)
integrated <- FindNeighbors(integrated, dims = 1:100)
integrated <- FindClusters(integrated, resolution = 0.5, k.param = 40)
integrated <- RunUMAP(integrated, dims = 1:100, n.neighbors = 25, n.components = 2, min.dist = 0.35)

```


```{r, fig.width = 12, fig.height = 10}
DimPlot(integrated, reduction="umap", label="TRUE", group.by="Patient")
#ggplot2::ggsave(file="umap/umap.png") 
```

```{r}
FeaturePlot(integrated, "ScrubletScores", cols = c("grey", "red"), order = T)

DimPlot(integrated, reduction = "umap", label=F, group.by = "isDoublet",cols = c("grey80","red"))
```


```{r}
FeaturePlot(integrated, features = "DoubletScores", pt.size = 0.01)
```



```{r}
markery_oka <- c("PDE6A", "ABCA8", "GRM6", "OPN1MW", "GAD1", "SLC6A9", "GRIK1", "PECAM1", "ONECUT1", "POU4F1", "C1QA", "ACTA2", "PDGFRA", "RPE65")
for (gene in markery_oka){
  print(gene)
  FeaturePlot(integrated, gene)
  file_path <- paste0("markers/", gene, ".png")
  ggplot2::ggsave(file=file_path)

}
```

```{r}
SaveH5Seurat(integrated, "integrated.h5Seurat")

```
