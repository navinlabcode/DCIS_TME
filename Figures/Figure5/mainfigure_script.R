
##Some NK/T cells have low transcripts (CD4, CD8 gene expression). Here we are using random forest to predict the cell types for these cells.
nkt_cell_anno <- rbind(data.frame(cellstate = cd8_all$cellstate), data.frame(cellstate = cd4_all$cellstate), data.frame(cellstate = nk_all_reint$cellstate))

t_cell_anno <- nkt_cell_anno1 %>% dplyr::filter(!cellstate %in% c('NK-GZMK', 'NK-FCGR3A', 'ILC1', 'ILC3'))
t_cell_anno$celltype <- ifelse(t_cell_anno$cellstate %in% c('CD8-EM', 'CD8-RM', 'CD8-Naive', 'gdT', 'CD8-Cytotoxic','CD8-stress', 'CD8-Exhausted','MAIT',  'CD8-Cycling', 'CD8-Interferon', 'pt-specific'), 'CD8T', "CD4T")

t_inte <- RunUMAP(t_inte, reduction = "pca", dims = 1:25)

DefaultAssay(t_inte) <- "integrated"

##first predict CD4 and CD8 cells (include CD4 and CD8 genes in PCs)
t_inte <- AddMetaData(t_inte, metadata = t_cell_anno, col.name = "celltype")
t_inte_anno <- subset(t_inte, celltype != "NA")

DefaultAssay(t_inte_anno) <- "RNA" 

t_rf <- randomForest::randomForest(x = as.matrix(Embeddings(t_inte_anno, 'pca')[, 1:30]), y = factor(t_inte_anno$celltype))
t_rf_errate <- t_rf$err.rate[, -1] %>% melt %>% set_colnames(c('Trees', 'CellType', 'ErrorRate'))
t_inte_anno$celltype_rf <- t_rf$predicted

t_inte_anno$celltype_rf_dif <- ifelse(t_inte_anno$celltype == t_inte_anno$celltype_rf, "yes", "no")
t_inte_anno_conf <- subset(t_inte_anno, celltype_rf_dif %in% c('yes'))
t_inte_anno_conf_id_label <- data.frame(anno = t_inte_anno_conf$celltype)

##predict unannotated cells
t_inte_unanno <- subset(t_inte, cells = c(rownames(t_inte_anno@meta.data)), invert = T)
DefaultAssay(t_inte_unanno) <- "RNA"
t_rf_pred <- predict(t_rf, as.matrix(Embeddings(t_inte_unanno, 'pca')[, 1:30]), type = c('prob'))

t_rf_pred_rf_raw <- t_rf_pred %>% melt %>% set_colnames(c('Cell', 'CellType', 'Prob'))
t_rf_pred_rf_raw_max_prob <- t_rf_pred_rf_raw %>% group_by(Cell) %>% slice(which.max(Prob))

t_rf_pred_rf_raw_max_prob_fil <- t_rf_pred_rf_raw_max_prob %>% filter(Prob > 0.6)
t_inte_unanno_conf_id_label <- data.frame(anno = t_rf_pred_rf_raw_max_prob_fil$CellType)
rownames(t_inte_unanno_conf_id_label) <- t_rf_pred_rf_raw_max_prob_fil$Cell

##only keep anno (with high conf) and predict cells with score > 0.6
t_inte_all_conf_id_label <- rbind(t_inte_unanno_conf_id_label, t_inte_anno_conf_id_label)
t_inte <- AddMetaData(t_inte, metadata = t_inte_all_conf_id_label, col.name = "celltype_anno_conf")
t_inte_fil <- subset(t_inte, celltype_anno_conf %in% c('CD4T', 'CD8T'))

t_inte_fil <- ScaleData(t_inte_fil, verbose = FALSE)
t_inte_fil <- RunPCA(t_inte_fil, npcs = 30, verbose = FALSE)
t_inte_fil <- RunUMAP(t_inte_fil, reduction = "pca", dims = 1:25)

