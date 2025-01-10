library(Seurat)
library(dplyr)
library(future)

plan("multiprocess", workers = 5)

setwd("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/integrated")
options(future.globals.maxSize = 50000 * 1024^2)

##normal
normal <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/normal/Normal_all.rds")
DefaultAssay(normal) <- "RNA"
normal_sub <- subset(normal, orig.ident %in% c("BCMHBCA83L_3h", "BCMHBCA15L_3h", "BCMHBCA82R_24hTis_4h", "BCMHBCA85L_3h", "BCMHBCA12L_3h", "BCMHBCA81R_3h", "BCMHBCA02_RIM", "BCMHBCA32L", "BCMHBCA10L", "BCMHBCA22R_4h", "BCMHBCA35R", "BCMHBCA20L_3h", "BCMHBCA21L_24hTis_3h", "BCMHBCA09L", "BCMHBCA72R_3h"))
dim(normal_sub)

##DCIS
dcis1 <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/DCIS/ECIS_all.rds")
DefaultAssay(dcis1) <- "RNA"
dcis1_sub <- subset(dcis1, orig.ident %in% c("ECIS01T", "ECIS02T", "ECIS05T", "ECIS06T", "ECIS07T", "ECIS09T", "ECIS11T", "ECIS12T", "ECIS14T", "ECIS15T", "ECIS39T", "ECIS40T", "ECIS41T", "ECIS43T", "ECIS44T", "ECIS45T", "ECIS46T", "ECIS48T", "ECIS49T_cryo", "ECIS50T", "ECIS51T", "ECIS52T"))
dim(dcis1_sub)

dcis2 <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/DCIS/BCMDCIS_all.rds")
DefaultAssay(dcis2) <- "RNA"
dcis2_sub <- subset(dcis2, orig.ident %in% c("BCMDCIS01T_pc1", "BCMDCIS02T", "BCMDCIS03T", "BCMDCIS04T", "BCMDCIS05T", "BCMDCIS06T", "BCMDCIS07T", "BCMDCIS08T", "BCMDCIS09T", "BCMDCIS10T", "BCMDCIS35T", "BCMDCIS36T", "BCMDCIS37T_cryo", "BCMDCIS38T_24hTis", "BCMDCIS39T", "BCMDCIS40T_24hTis", "BCMDCIS41T", "BCMDCIS42T_24htis", "BCMDCIS43T", "BCMDCIS44T", "BCMDCIS46T", "BCMDCIS47T", "BCMDCIS48T_24htis", "BCMDCIS49T_24htis", "BCMDCIS50T", "BCMDCIS51T", "BCMDCIS52T"))
dim(dcis2_sub)

##IDC
idc1 <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/IDC/idc_er_ng_seurat.rds")
DefaultAssay(idc1) <- "RNA"
idc1_sub <- subset(idc1, orig.ident %in% c("CID4461", "CID4463", "CID4530N", "CID4040"))
dim(idc1_sub)

idc2 <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/IDC/idc_er_nm_seurat_updated.rds")
DefaultAssay(idc2) <- "RNA"
idc2_sub <- subset(idc2, orig.ident %in% c("BIOKEY_3", "BIOKEY_5", "BIOKEY_30", "BIOKEY_12", "BIOKEY_20", "BIOKEY_22", "BIOKEY_21"))
dim(idc2_sub)

idc5 <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/IDC/idc_er_nc_seurat_updated.rds")
DefaultAssay(idc5) <- "RNA"
idc5_sub <- subset(idc5, orig.ident %in% c("TBB011", "TBB075", "TBB129"))
dim(idc5_sub)

idc3 <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/IDC/idc_er_hbca.rds")
DefaultAssay(idc3) <- "RNA"
idc3_sub <- subset(idc3, orig.ident %in% c("HBCA07T", "HBCA16T"))
dim(idc3_sub)

idc4 <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/IDC/idc_er_kh.rds")
DefaultAssay(idc4) <- "RNA"
idc4_sub <- subset(idc4, orig.ident %in% c("ER10", "ER11", "ER13", "ER14T", "ER15T"))
dim(idc4_sub)

rm(normal, dcis1, dcis2, idc1, idc2, idc3, idc4, idc5)

#hbca_dcis_ibc <- merge(normal, c(dcis1, dcis2, idc1, idc2, idc3, idc4, idc5))
message("merging samples......")
hbca_dcis_ibc <- merge(normal_sub, c(dcis1_sub, dcis2_sub, idc1_sub, idc2_sub, idc3_sub, idc4_sub, idc5_sub))

rm(normal_sub, dcis1_sub, dcis2_sub, idc1_sub, idc2_sub, idc3_sub, idc4_sub, idc5_sub)

gc()

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

saveRDS(hbca_dcis_ibc.combined, "hbca_dcis_ibc.integrated.b1.rds")

hbca_dcis_ibc_sub <- hbca_dcis_ibc.combined[, sample(ncol(hbca_dcis_ibc.combined), 26000)]
saveRDS(hbca_dcis_ibc_sub, "hbca_dcis_ibc.integrated.b1.subsample.rds")

hbca_dcis_ibc.combined.sub <- subset(hbca_dcis_ibc.combined, downsample = 4000)
saveRDS(hbca_dcis_ibc.combined, "hbca_dcis_ibc.integrated.b1.subclu.rds")

rm(list=ls())
gc()

