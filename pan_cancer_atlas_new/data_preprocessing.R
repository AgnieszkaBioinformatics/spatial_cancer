library(data.table)
library(dplyr)
library(Matrix)
library(glue)
library(readr)
#library(GEOquery)
setwd("/mnt/Work/amichalak")



                              ######## 1st publication ########

## count matrix
mtx_files_breast <- list.files(path = "pub1/breast/counts", 
                               pattern = "\\.mtx$", full.names = TRUE)

read_matrix <- function(index) {
  name <- mtx_files_breast[index]
  file_path <- paste0(glue('{getwd()}/'), name)
  print(name)
  print(file_path)
  matrix <- readMM(file_path)
  return(matrix)
}


####### breast
list_breast_matrices <- lapply(1:length(mtx_files_breast), read_matrix)

## GSE110686
mtx_list <- list()

for (channel in c(1, 2)) {
  # Construct file paths
  base_path <- paste0(getwd(), "/pub1/breast/GSE110686/GSM3011853_tils20_channel", channel)
  
  # Read matrix market files
  mtx <- readMM(paste0(base_path, "_matrix.mtx"))
  genes <- read_tsv(paste0(base_path, "_genes.tsv"),
                    col_names = c("gene_id", "gene_name"),
                    show_col_types = FALSE)
  cell_ids <- read_tsv(paste0(base_path, "_barcodes.tsv"),
                       col_names = "cell_id",
                       show_col_types = FALSE)
  
  # Set matrix dimensions
  rownames(mtx) <- genes$gene_id
  colnames(mtx) <- cell_ids$cell_id
  
  # Store in list
  mtx_list[[paste0("channel", channel)]] <- mtx
}


  