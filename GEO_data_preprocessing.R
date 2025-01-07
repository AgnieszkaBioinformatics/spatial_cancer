```{r}
library(data.table)
library(dplyr)
library(Matrix)
library(GEOquery)
```

```{r, setup, include=FALSE}
knitr::opts_knit$set(root.dir = '/mnt/Data/Projects/EyeAtlas/')
```

```{r}
## count data
file_names <- c("macula_donor_1.csv", "macula_donor_2.csv", "macula_donor_3.csv", "macula_donor_4.csv", "macula_donor_5.csv", "macula_donor_6.csv", "macula_donor_7.csv",
                "peripheral_donor_1.csv", "peripheral_donor_2.csv", "peripheral_donor_3.csv", "peripheral_donor_4.csv", "peripheral_donor_5.csv", "peripheral_donor_6.csv", "peripheral_donor_7.csv" )

## metadata
gsm <- c("GSM4037981", "GSM4037982", "GSM4037983", "GSM4037987", "GSM4037988", "GSM4037989", "GSM4037990",
         "GSM4037984", "GSM4037985", "GSM4037986", "GSM4037991", "GSM4037992", "GSM4037993", "GSM4037994")

read_count_matrix <- function(index){
  name <- file_names[index]
  file_path <- paste0("healthy_eye/GSE135922/", name)
  count_matrix <- fread(file_path)
  #count_matrix <- tibble::column_to_rownames(count_matrix, colnames(count_matrix)[1])
  colnames(count_matrix) <- paste("GSE135922", gsm[index], colnames(count_matrix), sep="_")
  colnames(count_matrix)[1]<-paste0("genes", index)
  return(count_matrix)
}

## read in matrices -> list of matrices
list_of_count_matrixes <- lapply(1:length(file_names), read_count_matrix)

```


```{r}
## merge matrices by gene
merged_count_matrix <- list_of_count_matrixes[[1]]
for (i in 2:length(file_names)){
  nazwa_kolumny <- paste0("genes", i)
  merged_count_matrix <- full_join(merged_count_matrix, list_of_count_matrixes[[i]], by = c("genes1" = nazwa_kolumny))
}

merged_count_matrix[is.na(merged_count_matrix)] <- 0
merged_count_matrix <- tibble::column_to_rownames(merged_count_matrix, "genes1")
matrix_merged <- as.matrix(merged_count_matrix)
sparse_merged <- as(matrix_merged, "sparseMatrix")

## save files
writeMM(obj = sparse_merged, file="atlas_healthy_eye/count_matrix/GSE135922/matrix_GSE135922.mtx")
write(x = rownames(merged_count_matrix), file = "atlas_healthy_eye/count_matrix/GSE135922/genes_GSE135922.tsv")
write(x = colnames(sparse_merged), file = "atlas_healthy_eye/count_matrix/GSE135922/barcodes_GSE135922.tsv")
```

```{r}
## read in metadata - only for one record???
geo_metadata <- getGEO("GSE135922")
metadata <- matrix(,length(colnames(merged_count_matrix)), 11)   # empty matrix at size colnames(merged_counts) x 11
metadata <- data.frame(metadata)
row.names(metadata) <- colnames(sparse_merged)
colnames(metadata) <- c("GSE", "GSM", "PMID", "Tissue", "Cell Annotation", "Region", "Patient", "Sex", "Age", "Ethnicity", "Dissected Tissue")

metadata[,"GSE"] <- "GSE135922"
metadata[,"Tissue"] <- "RPE-choroid"
metadata[,"PMID"] <- "31712411, 35705048"
for (i in 1:nrow(metadata)){
  gsm <- rownames(metadata)[i]  #get the row name
  gsm <- substr(gsm, 11, 20)  #extract substring of gsm
  index <- which(geo_metadata[["GSE135922_series_matrix.txt.gz"]]@phenoData@data[["geo_accession"]]==gsm)  #find gsm in geo_metadata
  metadata[i, "GSM"] <- gsm  # fill GSM column with extracted gsm
  patient_number <- geo_metadata[["GSE135922_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.1"]][index]  # extract the patient number
  patient_number <- sub("donor: ", "", patient_number)
  metadata[i, "Patient"] <-paste("patient", patient_number, sep=" ")
  metadata[i, "Region"] <- geo_metadata[["GSE135922_series_matrix.txt.gz"]]@phenoData@data[["location:ch1"]][index]  # extract the chromosome location and fill the column
  if (patient_number=="1"){
    metadata[i, "Age"] <- "54"
    metadata[i, "Sex"] <- "male"
  } else if (patient_number=="2"){
    metadata[i, "Age"] <- "82"
    metadata[i, "Sex"] <- "female"
  }  else if (patient_number=="3"){
    metadata[i, "Age"] <- "79"
    metadata[i, "Sex"] <- "male"
  }  else if (patient_number=="4"){
    metadata[i, "Age"] <- "92"
    metadata[i, "Sex"] <- "female"
  }  else if (patient_number=="5"){
    metadata[i, "Age"] <- "80"
    metadata[i, "Sex"] <- "male"
  }  else if (patient_number=="6"){
    metadata[i, "Age"] <- "75"
    metadata[i, "Sex"] <- "female"
  }  else if (patient_number=="7"){
    metadata[i, "Age"] <- "77"
    metadata[i, "Sex"] <- "male"
  }
}

write.csv(metadata, "atlas_healthy_eye/metadata/GSE135922_metadata.csv")
