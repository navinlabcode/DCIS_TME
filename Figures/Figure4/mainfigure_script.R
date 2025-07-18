
###--------Figure 4a--------###
myeloid_all_filtered <- readRDS("myeloid.seurat.inte.proj.all.upd.rds")
myeloid_all_filtered <- subset(myeloid_all_filtered, nFeature_RNA > 500)
DefaultAssay(myeloid_all_filtered) <- "integrated"

myeloid_all_filtered <- ScaleData(myeloid_all_filtered, verbose = FALSE)
myeloid_all_filtered <- RunPCA(myeloid_all_filtered, npcs = 30, verbose = FALSE)
myeloid_all_filtered <- RunUMAP(myeloid_all_filtered, reduction = "pca", dims = 1:25)
#myeloid_all_filtered <- RunUMAP(myeloid_all_filtered, reduction = "pca", dims = 1:20)
DimPlot(myeloid_all_filtered, label = T, label.size = 5)

myeloid_all_filtered <- FindNeighbors(myeloid_all_filtered, dims = 1:20)

myeloid_all_filtered <- FindClusters(myeloid_all_filtered, resolution = 0.12)

DimPlot(myeloid_all_filtered, label = T, label.size = 5, cols = basal_col_temp)
#FeaturePlot(myeloid_all_filtered, c("nCount_RNA", "nFeature_RNA"))
VlnPlot(myeloid_all_filtered, c("nCount_RNA", "nFeature_RNA"),log = T, ncol = 1, pt.size = 0, cols=basal_col_temp)

DefaultAssay(myeloid_all_filtered) <- "RNA"
myeloid_all_filtered_clu_marker <- FindAllMarkers(subset(myeloid_all_filtered, downsample = 3000), logfc.threshold = 0.4, only.pos = T, slot = "data") #for cluster number iden
myeloid_all_filtered_marker_top10 <- myeloid_all_filtered_clu_marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=9, wt=avg_log2FC)

myeloid_all_filtered <- ScaleData(myeloid_all_filtered, features = rownames(myeloid_all_filtered))
DoHeatmap(subset(myeloid_all_filtered, downsample = 1000), features = myeloid_all_filtered_marker_top10$gene, size=4, assay = "RNA", group.colors = basal_col_temp)

FeaturePlot(myeloid_all_filtered, myeloid_all_filtered_marker_top10 %>% filter(cluster == 1) %>% pull(gene))

DimPlot(myeloid_all_filtered, label = T, label.size = 3.5, cols = mye_col, group.by = "cellstate")


###--------Figure 4b--------###
cell_stat_plot <- c('mye-cycling', 'macro-lipo', 'macro-c3', 'macro-lyve1')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% 
  rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + 
  geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + 
  rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")



###--------Figure 4c--------###
macrophage_all <- subset(myeloid_all_filtered, cellstate %in% c('Macro_C3', 'Macro_LYVE1', 'Macro_CCL4',  'Macro_APOC1', 'Macro_ISG15'))

VlnPlot(macrophage_all, c('APOC1','SPP1','APOE', 'LIPA' ,'FABP5', 'LPL', 'TREM2'), pt.size = 0, cols = mye_col[c(3:7)], stack = T, flip = T, fill.by = 'ident') & NoLegend()

VlnPlot(macrophage_all, c('FOLR2', 'SLC40A1', 'SELENOP', 'LYVE1', 'MRC1', 'CD163', 'MAF'), pt.size = 0, cols = mye_col[c(3:7)], stack = T, flip = T, fill.by = 'ident') & NoLegend()


###--------Figure 4d--------###
macrophage_signature <- myeloid_all_filtered_marker %>% filter(cluster %in% c('Macro_C3', 'Macro_LYVE1', 'Macro_CCL4',  'Macro_APOC1', 'Macro_ISG15')) %>% group_by(cluster) %>% filter(avg_log2FC > 0.5 & pct.2 < 0.7 & p_val_adj < 0.05) %>% top_n(n = 50, wt = avg_log2FC)

macrophage_df <- macrophage_signature[,6:7]
macrophage_dfsample <- split(macrophage_df$gene, macrophage_df$cluster)
{
  macrophage0_gene <- bitr(as.character(macrophage_dfsample$Macro_C3), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  macrophage1_gene <- bitr(as.character(macrophage_dfsample$Macro_LYVE1), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  macrophage2_gene <- bitr(as.character(macrophage_dfsample$Macro_CCL4), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  macrophage3_gene <- bitr(as.character(macrophage_dfsample$Macro_APOC1), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  macrophage4_gene <- bitr(as.character(macrophage_dfsample$Macro_ISG15), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
}
macrophageph_genes <- list(C3 = macrophage0_gene$ENTREZID, LYVE1 = macrophage1_gene$ENTREZID, CCL4 = macrophage2_gene$ENTREZID, APOC1 = macrophage3_gene$ENTREZID, ISG15 = macrophage4_gene$ENTREZID)
macrophageph_GOclusterplot <- compareCluster(geneCluster = macrophageph_genes, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP")
dotplot(macrophageph_GOclusterplot, showCategory = 5) + scale_y_discrete(labels=function(x) str_wrap(x, width=85)) + 
  scale_colour_gradientn(colours = (paletteer::paletteer_c("grDevices::Blue-Red 2", 50)))

macrophageph_GOclusterplot1 <- simplify(macrophageph_GOclusterplot, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(macrophageph_GOclusterplot1, showCategory = 5) + scale_y_discrete(labels=function(x) str_wrap(x, width=85)) + 
  scale_colour_gradientn(colours = (paletteer::paletteer_c("grDevices::Blue-Red 2", 50)))


###--------Figure 4 macrophage spatial--------###
dcis_300_xenium_srt_12cons_macrophage <- merge(dcis_300_xenium_srt_12cons_basal_rm, dcis_300_xenium_srt_12cons_lumhr, merge.dr = 'SP')
dcis_300_xenium_srt_12cons_macrophage$cellstate_ml_merge <- ifelse(dcis_300_xenium_srt_12cons_macrophage$cellstate_ml_merge %in% c("AMS_others", "others"), 'Basal_others', dcis_300_xenium_srt_12cons_macrophage$cellstate_ml_merge)
sample_check <- unique(dcis_300_xenium_srt_12cons_macrophage$sample)[2]

{
  macrophage_plot <- subset(dcis_300_xenium_srt_12cons_macrophage, sample == sample_check) 
  macrophage_plot <- subset(macrophage_plot, cellstate_ml_merge == 'Normal_like_lumhr', invert = T) 
  macrophage_plot_df <- data.frame(macrophage_plot@reductions$SP@cell.embeddings, macrophage_plot$cellstate_ml_merge)
  colnames(macrophage_plot_df) <- c('X', 'Y', 'cellanno')
  macrophage_plot_df$point_size <- ifelse(macrophage_plot_df$cellanno %in% c("Tumor_lumhr"), 0.01, 0.3)
  #macrophage_plot_df$point_size <- ifelse(macrophage_plot_df$cellanno %in% c("vein", "artery", "capillary"), 0.4, macrophage_plot_df$point_size)
  ggplot(macrophage_plot_df, aes(x = X, y = Y, color = cellanno)) +
    geom_point(aes(size = point_size)) + #scale_size_continuous(range = c(0.2, 1)) + 
    theme_classic() + theme(panel.background = element_rect(fill = '#000000FF'))+
    guides(colour = guide_legend(override.aes = list(size=3))) + scale_color_manual(values = xenium_basal_col) + scale_size_identity()   
}


##normal vs dcis macrophage comparison
dcis_myeloid_count <- data.frame(table(dcis_300_xenium_srt_12cons_myeloid$sample, dcis_300_xenium_srt_12cons_myeloid$cellstate_ml_com))
dcis_myeloid_count$type <- 'DCIS'
dcis_myeloid_count_prop <- dcis_myeloid_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))
dcis_myeloid_count_prop <- dcis_myeloid_count_prop %>% filter(!Var1 %in% c('bcmdcis49',  'bcmdcis62', 'ECIS51T'))

hbca_myeloid_count <- data.frame(table(hbca_300_xenium_srt_6cons_myeloid$sample, hbca_300_xenium_srt_6cons_myeloid$cellstate_ml_com))
hbca_myeloid_count$type <- 'HBCA'
hbca_myeloid_count_prop <- hbca_myeloid_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))

hbca_dcis_myeloid_prop <- rbind(dcis_myeloid_count_prop, hbca_myeloid_count_prop)
hbca_dcis_myeloid_prop$type <- factor(hbca_dcis_myeloid_prop$type, levels = c('HBCA', 'DCIS'))

ggboxplot(hbca_dcis_myeloid_prop %>% filter(Var2 %in% c('Macro_APOC1', 'Macro_LYVE1', 'Macro_C3')), x = "Var2", y = "Prop", color = "type", palette = tissue_type_col,#x = "supp", y = "len", fill = "#00AFBB", 
          add = "jitter") + stat_compare_means(aes(group = type), label = "p.format", method = "wilcox.test") + rotate_x_text(angle = 45, hjust = 1, vjust = 1)


##CellTrek##
DCIS_st_image <- png::readPNG("./spatial/tissue_hires_image.png")
DCIS.st <- readRDS("./DCIS.rds")

DefaultAssay(DCIS.st) <- "Spatial"
DCIS.st <- NormalizeData(DCIS.st) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)

dcis.combined.sub2 <- readRDS("./all_cell_merged_sub40_diet.rds")
dcis.combined.sub2@meta.data <- dcis.combined.sub2@meta.data[,-grep("^integrated_snn_res",colnames(dcis.combined.sub2@meta.data))]
DimPlot(dcis.combined.sub2, group.by = "cellStatedcis", cols = sample(colorRampPalette(brewer.pal(12,"Paired"))(51)))

DCIS.tumor.sub <- readRDS("./DCIST.rds")
DCIS.tumor.sub$cellstate <- 'tumor'
DCIS.tumor.sub$celltype <- 'tumor'

dcis.combined.sub3 <- merge(dcis.combined.sub2, DCIS.tumor.sub)
Idents(dcis.combined.sub3) <- dcis.combined.sub3$celltype
dcis.combined.sub3 <- NormalizeData(dcis.combined.sub3) %>% FindVariableFeatures() %>% ScaleData() %>% 
                        RunPCA(verbose = FALSE)
dcis.combined.sub3$id <- names(dcis.combined.sub3$orig.ident)

DCIS_traint <- CellTrek::traint(st_data=DCIS.st, sc_data=dcis.combined.sub3,
                                st_assay='Spatial', sc_assay='RNA', cell_names='celltype')
DimPlot(DCIS_traint, group.by='type')+DimPlot(DCIS_traint, group.by = "celltype", label = T)

DCIS_schart <- CellTrek::celltrek(st_sc_int=DCIS_traint, int_assay='traint', 
                                  reduction='pca', intp=T, intp_pnt=10000, intp_lin=F, nPCs=30,
                                  ntree=1000, dist_thresh=0.55, top_spot=5, spot_n=10, repel_r=5, repel_iter=10, keep_model=T)$celltrek

SpatialDimPlot(DCIS_schart, group.by = 'celltype', crop = F, pt.size.factor = 1)













