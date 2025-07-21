

# NK cell (final version)--------------------------------------------------------------
nk_all_final <- subset(tcell_all, idents = c(8,13)) 
DefaultAssay(nk_all_final) <- "integrated"

nk_all_final <- ScaleData(nk_all_final, verbose = FALSE)
nk_all_final <- RunPCA(nk_all_final, npcs = 30, verbose = FALSE)
nk_all_final <- RunUMAP(nk_all_final, reduction = "pca", dims = 1:10)
nk_all_final <- FindNeighbors(nk_all_final, reduction = "pca", dims = 1:10)
#nk_all_final <- FindNeighbors(nk_all_final, reduction = "pca", dims = 1:30, k.param = 10)
nk_all_final <- FindClusters(nk_all_final, resolution = 0.3)

nk_col_temp <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(nk_all_final@active.ident)))
DimPlot(nk_all_final, label = T, label.size = 5, cols = nk_col_temp)
VlnPlot(nk_all_final, c("nCount_RNA", "nFeature_RNA"),log = T, ncol = 1, pt.size = 0, cols=nk_col_temp)

DefaultAssay(nk_all_final) <- "RNA"
nk_all_final_marker <- FindAllMarkers(subset(nk_all_final, downsample = 3000), logfc.threshold = 0.4, only.pos = T, slot = "data") #for cluster number iden
#nk_all_final_marker <- FindAllMarkers(nk_all_final, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden
nk_all_final_marker_top10 <- nk_all_final_marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=8, wt=avg_log2FC)
nk_all_final <- ScaleData(nk_all_final, features = rownames(nk_all_final))
DoHeatmap(subset(nk_all_final, downsample = 500), features = nk_all_final_marker_top10$gene, size=4, assay = "RNA", group.colors = nk_col_temp)

DimPlot(nk_all_final, ncol = 3, split.by = 'project', cols = nk_col_temp, label = T)

nk_all_final.rm <- subset(nk_all_final, idents = c(1,2,5))
#nk_all_final.rm <- subset(nk_all_final, idents = c(0:14,16,18))
nk_all_final <- nk_all_final.rm

nk_all_final <- AddMetaData(nk_all_final, metadata = nkt_sc.prediction2)
nk_all_final_ze_col <- sample(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(nk_all_final$pruned.labels.y))))
#DimPlot(nk_all_final, group.by = 'pruned.labels.x', label = T, cols = nk_all_final_ze_col) 
DimPlot(nk_all_final, group.by = 'pruned.labels.y', cols = nk_all_final_ze_col)
DimPlot(nk_all_final, group.by = 'orig.ident', cols = sample(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(nk_all_final$orig.ident)))), label = T)

saveRDS(nk_all_final, "nk_all_final_raw.rds")

###Re-integrate the NK data
nk_all_final_reint <- readRDS("nk.seurat.inte.proj.all.final.rds")
DefaultAssay(nk_all_final_reint) <- "integrated"
nk_all_final_reint <- RunUMAP(nk_all_final_reint, reduction = "pca", dims = 1:10) 
nk_all_final_reint <- FindNeighbors(nk_all_final_reint, reduction = "pca", dims = 1:10)
#nk_all_final <- FindNeighbors(nk_all_final, reduction = "pca", dims = 1:30, k.param = 10)
nk_all_final_reint <- FindClusters(nk_all_final_reint, resolution = 0.25)
nk_col_temp <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(nk_all_final_reint@active.ident)))
DimPlot(nk_all_final_reint, label = T, label.size = 5, cols = nk_col_temp)
VlnPlot(nk_all_final_reint, c("nCount_RNA", "nFeature_RNA"),log = T, ncol = 1, pt.size = 0, cols=nk_col_temp)
DimPlot(nk_all_final_reint, ncol = 3, split.by = 'project', cols = nk_col_temp, label = T)
DimPlot(nk_all_final_reint, group.by = 'project')

DefaultAssay(nk_all_final_reint) <- "RNA"
FeaturePlot(nk_all_final_reint, c('rna_SELL', 'rna_IL7R'), order = T)
FeaturePlot(nk_all_final_reint, c('RORC', 'KIT', 'IL23R', 'CXCR3', 'IFNG', 'GATA3', 'RORA','SELL', 'IL7R','ITGA1', 'ITGAE', 'GZMA', 'GZMB', 'PRF1', 'GNLY', 'XCL1'), order = T)
FeaturePlot(nk_all_final_reint, c('CD3D', 'FCGR3A', 'NCAM1', 'KLRB1', 'KLRD1'))
nk_all_final_reint_marker <- FindAllMarkers(subset(nk_all_final_reint, downsample = 3000), logfc.threshold = 0.4, only.pos = T, slot = "data") #for cluster number iden
#nk_all_final_reint_marker <- FindAllMarkers(nk_all_final, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden
nk_all_final_reint_marker_top10 <- nk_all_final_reint_marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene) & !grepl("^HSP",gene)) %>% top_n(n=8, wt=avg_log2FC)
nk_all_final_reint_marker_top10 <- nk_all_final_reint_marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=8, wt=avg_log2FC)
nk_all_final_reint <- ScaleData(nk_all_final_reint, features = rownames(nk_all_final_reint))
DoHeatmap(subset(nk_all_final_reint, downsample = 1000), features = nk_all_final_reint_marker_top10$gene, size=4, assay = "RNA", group.colors = nk_col_temp)

nk_all_final_reint.rm <- subset(nk_all_final_reint, idents = c(0:4,7,8))
nk_all_final_reint <- nk_all_final_reint.rm

##final
{
  nk_all_final_reint$cellstate <- "NULL"
  nk_all_final_reint$cellstate[nk_all_final_reint$integrated_snn_res.0.25 %in% c("0")] = "NK-FCGR3A"
  nk_all_final_reint$cellstate[nk_all_final_reint$integrated_snn_res.0.25 %in% c("1")] = "ILC1"
  nk_all_final_reint$cellstate[nk_all_final_reint$integrated_snn_res.0.25 %in% c("2")] = "NK-GZMK"
  nk_all_final_reint$cellstate[nk_all_final_reint$integrated_snn_res.0.25 %in% c("3")] = "ILC3"
}

nk_all_final_reint$cellstate <- factor(nk_all_final_reint$cellstate, levels = c("NK-GZMK", "NK-FCGR3A", "ILC1", "ILC3"))

Idents(nk_all_final_reint) <- nk_all_final_reint$cellstate

DimPlot(nk_all_final_reint, cols = nk_col, label = T)

VlnPlot(nk_all_final_reint, c("nCount_RNA", "nFeature_RNA"),log = T, ncol = 1, pt.size = 0, cols=nk_col)


##Run DE genes
DefaultAssay(nk_all_final_reint) <- "RNA"

nk_all_final_reint_marker <- FindAllMarkers(nk_all_final_reint, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden
#nk_all_final_reint_marker <- FindAllMarkers(nk_all_final_reint, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden
nk_all_final_reint_marker_top10 <- nk_all_final_reint_marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% filter(pct.2 < 0.35) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene) & !grepl("^HSP",gene)) %>% top_n(n=10, wt=avg_log2FC)
nk_all_final_reint <- ScaleData(nk_all_final_reint, features = rownames(nk_all_final_reint))

DoHeatmap(subset(nk_all_final_reint, downsample = 800), features = nk_all_final_reint_marker_top10$gene, size=4, assay = "RNA", group.colors = nk_col)

DimPlot(nk_all_final_reint, group.by = 'project', cols = proj_col)

nk_all_final_reint$celltype <- "Tcells"

DimPlot(nk_all_final_reint, group.by = 'orig.ident', cols = sample(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(nk_all_final_reint$orig.ident)))), label = T)
DimPlot(nk_all_final_reint, cols = nk_col, split.by = 'project', ncol = 3)


#NK cell state comparison -----------------------------------------------------
cell_stat_plot <- c('ilc1', 'ilc3', 'nk-fcgr3a', 'nk-gzmk')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% 
  rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + 
  geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + 
  rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")














