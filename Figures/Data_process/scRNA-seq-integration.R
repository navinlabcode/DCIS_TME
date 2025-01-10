library(Seurat)
library(dplyr)
library(future)

plan("multiprocess", workers = 5)

options(future.globals.maxSize = 50000 * 1024^2)

##normal
normal <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/normal/Normal_all.rds")
DefaultAssay(normal) <- "RNA"

##DCIS
dcis1 <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/DCIS/ECIS_all.rds")
DefaultAssay(dcis1) <- "RNA"

##IDC
idc1 <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/IDC/idc_er_ng_seurat.rds")
DefaultAssay(idc1) <- "RNA"

hbca_dcis_ibc <- merge(normal_sub, c(dcis1, idc1))

message("normalizing......")
hbca_dcis_ibc.list <- SplitObject(hbca_dcis_ibc, split.by = "project")
hbca_dcis_ibc.list <- lapply(X = hbca_dcis_ibc.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = hbca_dcis_ibc.list)

k.filter <- min(200, min(sapply(hbca_dcis_ibc.list, ncol)))

hbca_dcis_ibc.anchors <- FindIntegrationAnchors(object.list = hbca_dcis_ibc.list, anchor.features = features, k.filter = k.filter)

hbca_dcis_ibc.combined <- IntegrateData(anchorset = hbca_dcis_ibc.anchors)
hbca_dcis_ibc.combined <- ScaleData(hbca_dcis_ibc.combined, verbose = FALSE)
hbca_dcis_ibc.combined <- RunPCA(hbca_dcis_ibc.combined, npcs = 30, verbose = FALSE)
hbca_dcis_ibc.combined <- RunUMAP(hbca_dcis_ibc.combined, reduction = "pca", dims = 1:20)
hbca_dcis_ibc.combined <- FindNeighbors(hbca_dcis_ibc.combined, reduction = "pca", dims = 1:20)
hbca_dcis_ibc.combined <- FindClusters(hbca_dcis_ibc.combined, resolution = 0.5)

all_cell <- hbca_dcis_ibc.combined

saveRDS(all_cell, "all_cell.rds")

