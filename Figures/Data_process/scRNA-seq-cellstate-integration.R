library(Seurat)
library(dplyr)
library(future)

plan("multiprocess", workers = 4)

setwd("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/integrated")
options(future.globals.maxSize = 50000 * 1024^2)

hbca_dcis_ibc.combined.b1 <- readRDS("hbca_dcis_ibc.integrated.b1.res0.23.rds")
b1_obj <- subset(hbca_dcis_ibc.combined.b1, idents = c(9))
DefaultAssay(b1_obj) <- "RNA"

hbca_dcis_ibc.combined.b2 <- readRDS("hbca_dcis_ibc.integrated.b2.res0.23.rds")
b2_obj <- subset(hbca_dcis_ibc.combined.b2, idents = c(6))
DefaultAssay(b2_obj) <- "RNA"

all_obj <- merge(b1_obj, b2_obj)

DefaultAssay(all_obj) <- "RNA"

obj.list <- SplitObject(all_obj, split.by = "project")

obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = obj.list)
k.filter <- min(200, min(sapply(obj.list, ncol)))
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, k.filter = k.filter)
obj.combined <- IntegrateData(anchorset = obj.anchors)
obj.combined <- ScaleData(obj.combined, verbose = FALSE)
obj.combined <- RunPCA(obj.combined, npcs = 30, verbose = FALSE)


saveRDS(obj.combined, "cellstate.seurat.inte.proj.all.upd.rds")
