
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


# CD4 T cells ---------------------------------------------------
cd4_all_final <- subset(t_inte_fil, celltype_anno_conf %in% c('CD4T'))
DefaultAssay(cd4_all_final) <- "integrated"

cd4_all_final <- ScaleData(cd4_all_final, verbose = FALSE)
cd4_all_final <- RunPCA(cd4_all_final, npcs = 30, verbose = FALSE)
cd4_all_final <- RunUMAP(cd4_all_final, reduction = "pca", dims = 1:30)
cd4_all_final <- FindNeighbors(cd4_all_final, reduction = "pca", dims = 1:30) #k.para = 10
cd4_all_final <- FindClusters(cd4_all_final, resolution = 0.29)

DimPlot(cd4_all_final, label = T, label.size = 5, cols = cd4_col_temp)
VlnPlot(cd4_all_final, c("nCount_RNA", "nFeature_RNA"),log = T, ncol = 1, pt.size = 0, cols=cd4_col_temp)

DimPlot(cd4_all_final, split.by = 'project', ncol = 3, label = T)

DefaultAssay(cd4_all_final) <- "RNA"
cd4_all_final_marker <- FindAllMarkers(subset(cd4_all_final, downsample = 900), logfc.threshold = 0.4, only.pos = T, slot = "data") #for cluster number iden
cd4_all_final_marker_top10 <- cd4_all_final_marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=8, wt=avg_log2FC)
cd4_all_final <- ScaleData(cd4_all_final, features = rownames(cd4_all_final))



#CD4 T cell state comparison -----------------------------------------------------
cell_stat_plot <- c('cd4-cyc', 'cd4-ifn', 'cd4-tfh', 'treg', 'cd4-mem')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% 
  rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + 
  geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + 
  rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


# CD8 T cells----------------------------------------------------
cd8_all_final <- subset(t_inte_fil, celltype_anno_conf %in% c('CD8T'))
DefaultAssay(cd8_all_final) <- "integrated"

cd8_all_final <- ScaleData(cd8_all_final, verbose = FALSE)
cd8_all_final <- RunPCA(cd8_all_final, npcs = 30, verbose = FALSE)
cd8_all_final <- RunUMAP(cd8_all_final, reduction = "pca", dims = 1:20)
cd8_all_final <- FindNeighbors(cd8_all_final, reduction = "pca", dims = 1:20, k.param = 10)
cd8_all_final <- FindClusters(cd8_all_final, resolution = 0.47)

Idents(cd8_all_final) <- cd8_all_final$integrated_snn_res.0.47

DimPlot(cd8_all_final, label = T, label.size = 5, cols = cd8_col_temp)
VlnPlot(cd8_all_final, c("nCount_RNA", "nFeature_RNA"),log = T, ncol = 1, pt.size = 0, cols=cd8_col_temp)

DefaultAssay(cd8_all_final) <- "RNA"
cd8_all_final_marker <- FindAllMarkers(subset(cd8_all_final, downsample = 900), logfc.threshold = 0.4, only.pos = T, slot = "data") #for cluster number iden
cd8_all_final_marker_top10 <- cd8_all_final_marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=8, wt=avg_log2FC)
cd8_all_final <- ScaleData(cd8_all_final, features = rownames(cd8_all_final))
DoHeatmap(subset(cd8_all_final, downsample = 600), features = cd8_all_final_marker_top10$gene, size=4, group.colors = cd8_col_temp)



#CD8 T cell state comparison -----------------------------------------------------
cell_stat_plot <- c('cd8-ifn', 'cd8-trm-aust2')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% 
  rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + 
  geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + 
  rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")



# T cell spatial ----------------------------------------------------------
##proportion comparison between normal and dcis
dcis_tcell_count <- data.frame(table(dcis_300_xenium_srt_12cons_tcell$sample, dcis_300_xenium_srt_12cons_tcell$cellstate_ml_com))
dcis_tcell_count$type <- 'DCIS'
dcis_tcell_count_prop <- dcis_tcell_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))


hbca_tcell_count <- data.frame(table(hbca_300_xenium_srt_6cons_tcell$sample, hbca_300_xenium_srt_6cons_tcell$cellstate_ml_com))
hbca_tcell_count$type <- 'HBCA'
hbca_tcell_count_prop <- hbca_tcell_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))

hbca_dcis_tcell_prop <- rbind(dcis_tcell_count_prop, hbca_tcell_count_prop)
hbca_dcis_tcell_prop$type <- factor(hbca_dcis_tcell_prop$type, levels = c('HBCA', 'DCIS'))

ggboxplot(hbca_dcis_tcell_prop %>% filter(Var2 %in% c('IgA')), x = "Var2", y = "Prop", color = "type", palette = tissue_type_col,#x = "supp", y = "len", fill = "#00AFBB", 
          add = "jitter") + stat_compare_means(aes(group = type), label = "p.format", method = "wilcox.test") + rotate_x_text(angle = 45, hjust = 1, vjust = 1)




# B cells ------------------------------------------------------------
bcell_all_unfiltered <- readRDS("bcell.seurat.inte.proj.all.upd.rds")
DefaultAssay(bcell_all_unfiltered) <- "integrated"
bcell_all_unfiltered <- ScaleData(bcell_all_unfiltered, verbose = FALSE)
bcell_all_unfiltered <- RunPCA(bcell_all_unfiltered, npcs = 30, verbose = FALSE)
ElbowPlot(bcell_all_unfiltered, ndims = 30)
bcell_all_unfiltered <- RunUMAP(bcell_all_unfiltered, reduction = "pca", dims = 1:20)
#bcell_all_unfiltered <- RunUMAP(bcell_all_unfiltered, reduction = "pca", dims = 1:20)
DimPlot(bcell_all_unfiltered, label = T, label.size = 5)
DefaultAssay(bcell_all_unfiltered) <- "RNA"

DimPlot(bcell_all, group.by = 'project', cols = proj_col)
DimPlot(bcell_all, group.by = 'project', cols = proj_col, split.by = 'project', ncol = 3)
DimPlot(bcell_all, cols = bcell_col, split.by = 'project', ncol = 3)

DimPlot(subset(bcell_all, project == "NG_ER_IBC"), group.by = "celltype_subset", cols = sample(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(bcell_all$celltype_subset)))))
DimPlot(subset(bcell_all, project == "NM_ER_IBC"), group.by = "cellSubType",  cols = sample(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(bcell_all$cellSubType)))))
DimPlot(subset(bcell_all, project == "NC_ER_IBC"), group.by = "bcell_cluster", cols = sample(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(bcell_all$bcell_cluster)))))

bcell_all <- subset(bcell_all_unfiltered, nFeature_RNA > 500)
DefaultAssay(bcell_all) <- "integrated"
bcell_all <- ScaleData(bcell_all, verbose = FALSE)
bcell_all <- RunPCA(bcell_all, npcs = 30, verbose = FALSE)
ElbowPlot(bcell_all, ndims = 30)
bcell_all <- RunUMAP(bcell_all, reduction = "pca", dims = 1:20)
#bcell_all <- RunUMAP(bcell_all, reduction = "pca", dims = 1:20)
DimPlot(bcell_all, label = T, label.size = 5)

bcell_all <- FindNeighbors(bcell_all, dims = 1:20)

bcell_all <- FindClusters(bcell_all, resolution = 0.45) ##
DimPlot(bcell_all, label = T, label.size = 5)

VlnPlot(bcell_all, c("nCount_RNA", "nFeature_RNA"),log = T, ncol = 1, pt.size = 0, cols=bcell_col_temp)

DefaultAssay(bcell_all) <- "RNA"
bcell_all_marker <- FindAllMarkers(subset(bcell_all, downsample = 3000), logfc.threshold = 0.4, only.pos = T, slot = "data") #for cluster number iden
#bcell_all_marker <- FindAllMarkers(bcell_all, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden
bcell_all_marker_top10 <- bcell_all_marker %>% filter(p_val_adj < 0.05) %>% filter(pct.1 > 0.3) %>% filter(pct.2 > 0.5) %>% group_by(cluster) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=12, wt=avg_logFC)
bcell_all <- ScaleData(bcell_all, features = rownames(bcell_all))
DoHeatmap(subset(bcell_all, downsample = 500), features = bcell_all_marker_top10$gene, size=4, assay = "RNA", group.colors = bcell_col_temp)

DimPlot(bcell_all, split.by = 'project', label = T, cols = bcell_col_temp, ncol = 3)

bcell_all.rm <- subset(bcell_all, idents = c(0:2,4))
bcell_all.rm <- subset(bcell_all, idents = c(0:7))#1 round
DimPlot(bcell_all.rm, cols=colorRampPalette(brewer.pal(12,"Paired"))(length(unique(bcell_all.rm@active.ident))), label = T)
bcell_all <- bcell_all.rm

DefaultAssay(bcell_all) <- "integrated"
bcell_all <- ScaleData(bcell_all, verbose = FALSE)
bcell_all <- RunPCA(bcell_all, npcs = 30, verbose = FALSE)
bcell_all <- RunUMAP(bcell_all, reduction = "pca", dims = 1:20)
DimPlot(bcell_all, label = T, label.size = 5, cols = bcell_col_temp)


#B cell state comparison -----------------------------------------------------
cell_stat_plot <- c('memoryB', 'plasma-igA')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% 
  rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + 
  geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + 
  rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


# B cell spatial ----------------------------------------------------------
##proportion comparison between normal and dcis
dcis_bcell_count <- data.frame(table(dcis_300_xenium_srt_12cons_bcell$sample, dcis_300_xenium_srt_12cons_bcell$cellstate_ml_com))
dcis_bcell_count$type <- 'DCIS'
dcis_bcell_count_prop <- dcis_bcell_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))


hbca_bcell_count <- data.frame(table(hbca_300_xenium_srt_6cons_bcell$sample, hbca_300_xenium_srt_6cons_bcell$cellstate_ml_com))
hbca_bcell_count$type <- 'HBCA'
hbca_bcell_count_prop <- hbca_bcell_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))

hbca_dcis_bcell_prop <- rbind(dcis_bcell_count_prop, hbca_bcell_count_prop)
hbca_dcis_bcell_prop$type <- factor(hbca_dcis_bcell_prop$type, levels = c('HBCA', 'DCIS'))

ggboxplot(hbca_dcis_bcell_prop %>% filter(Var2 %in% c('IgA')), x = "Var2", y = "Prop", color = "type", palette = tissue_type_col,#x = "supp", y = "len", fill = "#00AFBB", 
          add = "jitter") + stat_compare_means(aes(group = type), label = "p.format", method = "wilcox.test") + rotate_x_text(angle = 45, hjust = 1, vjust = 1)








