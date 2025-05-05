library(Seurat)
library(dplyr)

DCIS_file <- read.table("/./DCIS_all.txt")
DCIS_seurat <- list()

for(i in c(1:length(DCIS_file$V1))){
    sample <- as.character(DCIS_file$V1[i])
    sample_sc <- CreateSeuratObject(counts=Read10X(data.dir = sample))
    sample_sc[["percent.mt"]] <- PercentageFeatureSet(sample_sc, pattern = "^MT-")
    sample_sc[["percent.rp"]] <- PercentageFeatureSet(sample_sc, pattern = "^RP[L|S]")
    sample_sc_sub <- subset(sample_sc, subset= percent.mt < 20 & nFeature_RNA > 200)
    sample_name <- basename(dirname(dirname(sample)))
    sample_sc_sub$orig.ident <- sample_name
    sample_sc_sub$project <- "DCIS"
    DCIS_seurat[[i]] <- sample_sc_sub
}

DCIS_seurat_merged <- merge(x = DCIS_seurat[[1]], y = DCIS_seurat[c(2: length(DCIS_seurat))])

saveRDS(DCIS_seurat_merged, "/./DCIS_all.rds")
