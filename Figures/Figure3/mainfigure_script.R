library(Seurat)
library(dplyr)
library(ggpubr)


# figure 3a basal/myoepithelial analysis -----------------------------------------------------
basal_all <- readRDS("basal.seurat.inte.proj.all.upd.rds")

basal_all_unfiltered <- readRDS("basal.seurat.inte.proj.all.upd.rds")
DefaultAssay(basal_all_unfiltered) <- "integrated"

basal_all <- subset(basal_all_unfiltered, nFeature_RNA > 500)
DefaultAssay(basal_all) <- "integrated"
basal_all <- ScaleData(basal_all, verbose = FALSE)
basal_all <- RunPCA(basal_all, npcs = 30, verbose = FALSE)
basal_all <- RunUMAP(basal_all, reduction = "pca", dims = 1:25)
#basal_all <- RunUMAP(basal_all, reduction = "pca", dims = 1:20)
DimPlot(basal_all, label = T, label.size = 5)

basal_all <- FindNeighbors(basal_all, dims = 1:20)

basal_all <- FindClusters(basal_all, resolution = 0.12)

DimPlot(basal_all, label = T, label.size = 5, cols = basal_col_temp)
#FeaturePlot(basal_all, c("nCount_RNA", "nFeature_RNA"))
VlnPlot(basal_all, c("nCount_RNA", "nFeature_RNA"),log = T, ncol = 1, pt.size = 0, cols=basal_col_temp)

DefaultAssay(basal_all) <- "RNA"
basal_all_clu_marker <- FindAllMarkers(subset(basal_all, downsample = 3000), logfc.threshold = 0.4, only.pos = T, slot = "data") #for cluster number iden
basal_all_marker_top10 <- basal_all_clu_marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=9, wt=avg_log2FC)

basal_all <- ScaleData(basal_all, features = rownames(basal_all))
DoHeatmap(subset(basal_all, downsample = 1000), features = basal_all_marker_top10$gene, size=4, assay = "RNA", group.colors = basal_col_temp)

DimPlot(basal_all, split.by = 'project', label = T, cols = basal_col_temp, ncol = 3)

FeaturePlot(basal_all, basal_all_marker_top10 %>% filter(cluster == 1) %>% pull(gene))


# figure 3b basal top markers -----------------------------------------------------
DotPlot(basal_all, basal_all_marker_top10)


# figure 3c basal cell state comparison -----------------------------------------------------
cell_stat_plot <- c('basal-major', 'cab')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% 
  rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + 
  geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + 
  rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


# figure 3d enrichment analysis -----------------------------------------------------
basal_signature <- basal_all_marker %>% group_by(cluster) %>% filter(avg_log2FC > 0.5 & pct.2 < 0.7 & p_val_adj < 0.05) %>% top_n(n = 50, wt = avg_log2FC)
basal_df <- basal_signature[,6:7]
basal_dfsample <- split(basal_df$gene, basal_df$cluster)
{
  basal1_gene <- bitr(as.character(basal_dfsample$`Basal-major`), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  basal2_gene <- bitr(as.character(basal_dfsample$`Basal-interferon`), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  basal3_gene <- bitr(as.character(basal_dfsample$CAB), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
}
basalph_genes <- list(major = basal1_gene$ENTREZID, int = basal2_gene$ENTREZID, cab = basal3_gene$ENTREZID)
basalph_GOclusterplot <- compareCluster(geneCluster = basalph_genes, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP")
dotplot(basalph_GOclusterplot, showCategory = 8) + scale_y_discrete(labels=function(x) str_wrap(x, width=85)) +
  scale_colour_gradientn(colours = (paletteer::paletteer_c("grDevices::Blue-Red 2", 50)))

basalph_GOclusterplot1 <- simplify(basalph_GOclusterplot, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(basalph_GOclusterplot1, showCategory = 7) + scale_y_discrete(labels=function(x) str_wrap(x, width=80)) + 
  scale_colour_gradientn(colours = (paletteer::paletteer_c("grDevices::Blue-Red 2", 50)))


# figure 3e integrin markers -----------------------------------------------------
VlnPlot(basal_all, c('ITGB1', 'ITGA2', 'ITGA6', 'ITGA3', 'ITGB4', 'ITGBL1', 'ITGAV', 'ITGB6'), pt.size = 0, cols = basal_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()


# figure 3 basal spatial -----------------------------------------------------
## basal spatial plot
dcis_300_xenium_srt_12cons_tumor_basal <- merge(dcis_300_xenium_srt_12cons_basal_rm, dcis_300_xenium_srt_12cons_lumhr, merge.dr = 'SP')
dcis_300_xenium_srt_12cons_tumor_basal$cellstate_ml_merge <- ifelse(dcis_300_xenium_srt_12cons_tumor_basal$cellstate_ml_merge %in% c("AMS_others", "others"), 'Basal_others', dcis_300_xenium_srt_12cons_tumor_basal$cellstate_ml_merge)
sample_check <- unique(dcis_300_xenium_srt_12cons_tumor_basal$sample)[2]

{
  basal_plot <- subset(dcis_300_xenium_srt_12cons_tumor_basal, sample == sample_check) 
  basal_plot <- subset(basal_plot, cellstate_ml_merge == 'Normal_like_lumhr', invert = T) 
  basal_plot_df <- data.frame(basal_plot@reductions$SP@cell.embeddings, basal_plot$cellstate_ml_merge)
  colnames(basal_plot_df) <- c('X', 'Y', 'cellanno')
  basal_plot_df$point_size <- ifelse(basal_plot_df$cellanno %in% c("Tumor_lumhr"), 0.01, 0.3)
  #basal_plot_df$point_size <- ifelse(basal_plot_df$cellanno %in% c("vein", "artery", "capillary"), 0.4, basal_plot_df$point_size)
  ggplot(basal_plot_df, aes(x = X, y = Y, color = cellanno)) +
    geom_point(aes(size = point_size)) + #scale_size_continuous(range = c(0.2, 1)) + 
    theme_classic() + theme(panel.background = element_rect(fill = '#000000FF'))+
    guides(colour = guide_legend(override.aes = list(size=3))) + scale_color_manual(values = xenium_basal_col) + scale_size_identity()   
}



## normal vs dcis basal comparison
dcis_basal_count <- data.frame(table(dcis_300_xenium_srt_12cons_basal_rm$sample, dcis_300_xenium_srt_12cons_basal_rm$cellstate_ml_merge))
dcis_basal_count$Var2 <- as.character(dcis_basal_count$Var2)
dcis_basal_count$Var2 <- ifelse(dcis_basal_count$Var2 %in% c('AMS_others', 'others'), 'basal_others', dcis_basal_count$Var2)
dcis_basal_count$type <- 'DCIS'
dcis_basal_count_prop <- dcis_basal_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))

hbca_basal_count <- data.frame(table(hbca_300_xenium_srt_6cons_basal$sample, hbca_300_xenium_srt_6cons_basal$cellstate_ml_merge))
hbca_basal_count$Var2 <- as.character(hbca_basal_count$Var2)
hbca_basal_count$Var2 <- ifelse(hbca_basal_count$Var2 %in% c('AMS_others', 'others'), 'basal_others', hbca_basal_count$Var2)
hbca_basal_count$type <- 'HBCA'
hbca_basal_count_prop <- hbca_basal_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))

hbca_dcis_basal_prop <- rbind(dcis_basal_count_prop, hbca_basal_count_prop)
hbca_dcis_basal_prop$type <- factor(hbca_dcis_basal_prop$type, levels = c('HBCA', 'DCIS'))

ggboxplot(hbca_dcis_basal_prop %>% filter(Var2 %in% c('cab')), x = "Var2", y = "Prop", color = "type", palette = tissue_type_col,
          add = "jitter") + rotate_x_text(angle = 45, hjust = 1, vjust = 1) + ylim(0, 0.2)






