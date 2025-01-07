library(Seurat, lib.loc = "/mnt/Work/software/rlibraries")
library(data.table)
library(reticulate)
library(Matrix)
library(SeuratData)
library(SeuratDisk)

setwd("/mnt/Work/asiwiak/retina/ver4/")
files_objects <- list.files(path = "Seurat_orto/", pattern = "*.h5Seurat", full.names = TRUE)
seurat <- lapply(files_objects, LoadH5Seurat)

meta <- rep("chicken", nrow(seurat[[1]]))
seurat[[1]] <- AddMetaData(seurat[[1]], meta, col.name = "Species")

meta <- rep("cow", nrow(seurat[[2]]))
seurat[[2]] <- AddMetaData(seurat[[2]], meta, col.name = "Species")

meta <- rep("ferret", nrow(seurat[[3]]))
seurat[[3]] <- AddMetaData(seurat[[3]], meta, col.name = "Species")

meta <- rep("human", nrow(seurat[[4]]))
seurat[[4]] <- AddMetaData(seurat[[4]], meta, col.name = "Species")

meta <- rep("lamprey", nrow(seurat[[5]]))
seurat[[5]] <- AddMetaData(seurat[[5]], meta, col.name = "Species")

meta <- rep("lizard", nrow(seurat[[6]]))
seurat[[6]] <- AddMetaData(seurat[[6]], meta, col.name = "Species")

meta <- rep("macaque", nrow(seurat[[7]]))
seurat[[7]] <- AddMetaData(seurat[[7]], meta, col.name = "Species")

meta <- rep("marmoset", nrow(seurat[[8]]))
seurat[[8]] <- AddMetaData(seurat[[8]], meta, col.name = "Species")

meta <- rep("mouse", nrow(seurat[[9]]))
seurat[[9]] <- AddMetaData(seurat[[9]], meta, col.name = "Species")

meta <- rep("opossum", nrow(seurat[[10]]))
seurat[[10]] <- AddMetaData(seurat[[10]], meta, col.name = "Species")

meta <- rep("peromyscus", nrow(seurat[[11]]))
seurat[[11]] <- AddMetaData(seurat[[11]], meta, col.name = "Species")

meta <- rep("pig", nrow(seurat[[12]]))
seurat[[12]] <- AddMetaData(seurat[[12]], meta, col.name = "Species")

meta <- rep("rhabdomys", nrow(seurat[[13]]))
seurat[[13]] <- AddMetaData(seurat[[13]], meta, col.name = "Species")

meta <- rep("sheep", nrow(seurat[[14]]))
seurat[[14]] <- AddMetaData(seurat[[14]], meta, col.name = "Species")

meta <- rep("squirrel", nrow(seurat[[15]]))
seurat[[15]] <- AddMetaData(seurat[[15]], meta, col.name = "Species")

meta <- rep("tree shrew", nrow(seurat[[16]]))
seurat[[16]] <- AddMetaData(seurat[[16]], meta, col.name = "Species")

meta <- rep("zebrafish", nrow(seurat[[17]]))
seurat[[17]] <- AddMetaData(seurat[[17]], meta, col.name = "Species")

setwd('/mnt/Work/asiwiak/retina/ver4/scrublet/')

for (i in 1:length(seurat)) {
    samples <- unique(seurat[[i]]@meta.data$orig.ident)
    genes <- rownames(seurat[[i]])
    scores <- list()
    pred <- list()
    cells <- list()
    l = 1
    for (j in samples) {
        cellIDs <- rownames(subset(seurat[[i]]@meta.data, orig.ident == j))
        if (length(cellIDs) < 30) {
            next
        }
        mat <- seurat[[i]][, cellIDs]
        mat <- mat@assays$RNA$counts
        print("Writing data for Scrublet")
        writeMM(mat, file = "expression_matrix.mtx")

        scr <- reticulate::import("scrublet", convert = FALSE)    # for using python code
        np <- reticulate::import("numpy", convert = FALSE)
        sc <- reticulate::import("scipy.io", convert = FALSE)

        print("Reading data")
        ds <- sc$mmread("expression_matrix.mtx")
        ds <- ds$T$tocsc()    # sparse column 
        print("Running scrublet")
        scrub <- scr$Scrublet(ds, expected_doublet_rate = 0.05)
        res <- scrub$scrub_doublets(min_counts = 2L,         # dataframe with True and False values
                                    min_cells = 3L, 
                                    min_gene_variability_pctl = 85L, 
                                    n_prin_comps = 30L)
        scores[[l]] <- py_to_r(res)[[1]]    # convert pandas dataframe to r dataframe
        pred[[l]] <- py_to_r(res)[[2]]
        cells[[l]] <- cellIDs
        l = l + 1
        rm(mat); rm(ds); gc()
    }
    df <- data.frame(CellIDs = unlist(cells), 
                     ScrubletScores = unlist(scores), 
                     ScrubletPredictedDoublets = unlist(pred))
    filename <- paste0(i, "_scrublet_seurat.tsv")
    fwrite(as.data.frame(df), file = filename)
}
