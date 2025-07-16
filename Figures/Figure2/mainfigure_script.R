library(Seurat)
library(dplyr)


# fig2a fibro UMAP-------------------------------------------------------------------
fibro_all_unfiltered <- readRDS("fibro.seurat.inte.proj.all.upd.rds")
DefaultAssay(fibro_all_unfiltered) <- "integrated"

fibro_all <- subset(fibro_all_unfiltered, nFeature_RNA > 500)
DefaultAssay(fibro_all) <- "integrated"
fibro_all <- ScaleData(fibro_all, verbose = FALSE)
fibro_all <- RunPCA(fibro_all, npcs = 30, verbose = FALSE)
fibro_all <- RunUMAP(fibro_all, reduction = "pca", dims = 1:20)
fibro_all <- FindNeighbors(fibro_all, dims = 1:20)
fibro_all <- FindClusters(fibro_all, resolution = 0.2)
                                                                      

DimPlot(fibro_all)

# fig2b fibroblast cell state comparison---------------------------------------------
cell_stat_plot <- c('Fibro_FN1', 'Fibro_IGF1')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% 
            rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + 
      geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
      facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + 
      rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


# fig2c fibro markers-------------------------------------------------------------------
VlnPlot(fibro_all_diet, c('COL15A1', 'IGF1', 'DLK1', 'IGFBP2', 'ACTA2', 'FAP', 'TAGLN', 'POSTN'), pt.size = 0, 
        cols = fibro_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()



# fig2 fibroblast spatial-------------------------------------------------------------------
dcis_300_xenium_srt_12cons_fibroblast_basal <- merge(dcis_300_xenium_srt_12cons_basal_rm, c(dcis_300_xenium_srt_12cons_fibroblast, dcis_300_xenium_srt_12cons_lumhr), merge.dr = 'SP')

fibro_plot <- subset(dcis_300_xenium_srt_12cons_fibroblast_basal, cellstate_ml_com %in% c("basal", "fibro_sfrp4", "fibro_cfd",  "fibro_fn1",  "fibro_igf1",       
                                                                                          "fibro_cycling", "Tumor_lumhr") & sample == sample_check)  
fibro_plot_df <- data.frame(fibro_plot@reductions$SP@cell.embeddings, fibro_plot$cellstate_ml_com)
colnames(fibro_plot_df) <- c('X', 'Y', 'cellanno')
fibro_plot_df$point_size <- ifelse(fibro_plot_df$cellanno %in% c("basal", "Tumor_lumhr"), 0.06, 0.6)

ggplot(fibro_plot_df, aes(x = X, y = Y, color = cellanno)) +
  geom_point(aes(size = point_size)) + #scale_size_continuous(range = c(0.2, 1)) + 
  theme_classic() + theme(panel.background = element_rect(fill = '#000000FF'))+
  guides(colour = guide_legend(override.aes = list(size=3))) + scale_color_manual(values = fibro_xenium_col) + scale_size_identity()   


# normal vs dcis fibroblast
dcis_fibro_count <- data.frame(table(dcis_300_xenium_srt_12cons_fibroblast$sample, dcis_300_xenium_srt_12cons_fibroblast$cellstate_ml_com))
dcis_fibro_count$type <- 'DCIS'
dcis_fibro_count_prop <- dcis_fibro_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))

hbca_fibro_count <- data.frame(table(hbca_300_xenium_srt_6cons_fibro$sample, hbca_300_xenium_srt_6cons_fibro$cellstate_ml_com))
hbca_fibro_count$type <- 'HBCA'
hbca_fibro_count_prop <- hbca_fibro_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))

hbca_dcis_fibro_prop <- rbind(dcis_fibro_count_prop, hbca_fibro_count_prop)
hbca_dcis_fibro_prop$type <- factor(hbca_dcis_fibro_prop$type, levels = c('HBCA', 'DCIS'))

ggboxplot(hbca_dcis_fibro_prop %>% filter(Var2 %in% c('fibro_fn1', 'fibro_igf1')), x = "Var2", y = "Prop", color = "type", palette = tissue_type_col,#x = "supp", y = "len", fill = "#00AFBB", 
          add = "jitter") + stat_compare_means(aes(group = type), label = "p.format", method = "wilcox.test") + rotate_x_text(angle = 45, hjust = 1, vjust = 1)





# vascular analysis -----------------------------------------------------------
# fig2g vascular state umap-------------------------------------------------------------------
vas_all_unfiltered <- readRDS("vas.seurat.inte.proj.all.upd.rds")
DefaultAssay(vas_all_unfiltered) <- "integrated"

vas_all <- subset(vas_all_unfiltered, nFeature_RNA > 500)
DefaultAssay(vas_all) <- "integrated"
vas_all <- ScaleData(vas_all, verbose = FALSE)
vas_all <- RunPCA(vas_all, npcs = 30, verbose = FALSE)
vas_all <- RunUMAP(vas_all, reduction = "pca", dims = 1:20)
vas_all <- FindNeighbors(vas_all, dims = 1:20)
vas_all <- FindClusters(vas_all, resolution = 0.2)
                                                                      
DefaultAssay(vas_all) <- "RNA"
vas_all_marker <- FindAllMarkers(vas_all, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden

DimPlot(vas_all)


# fig2h vascular cell state comparison---------------------------------------------
cell_stat_plot <- c('capillary', 'tec', 'vas-cycling')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% 
  rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + 
  geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + 
  rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


# fig2i vascular state umap-------------------------------------------------------------------
VlnPlot(vas_all, c('APLNR', 'INSR', 'FLT1', 'KDR', 'NRP1', 'VWA1', 'COL4A1', 'COL4A2', 'COL6A2','LAMA4','LAMB1'), 
        pt.size = 0, cols = vas_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()



# fig2j vascular spatial-------------------------------------------------------------------
dcis_300_xenium_srt_12cons_vascular_basal <- merge(dcis_300_xenium_srt_12cons_basal_rm, c(dcis_300_xenium_srt_12cons_vascular, dcis_300_xenium_srt_12cons_lumhr), merge.dr = 'SP')

vascular_plot <- subset(dcis_300_xenium_srt_12cons_vascular_basal, cellstate_ml_com %in% c("basal", "artery", "vein", "tec", "capillary", "Tumor_lumhr") & sample == sample_check)  
vascular_plot_df <- data.frame(vascular_plot@reductions$SP@cell.embeddings, vascular_plot$cellstate_ml_com)
colnames(vascular_plot_df) <- c('X', 'Y', 'cellanno')
vascular_plot_df$point_size <- ifelse(vascular_plot_df$cellanno %in% c("basal", "Tumor_lumhr"), 0.06, 0.6)

ggplot(vascular_plot_df, aes(x = X, y = Y, color = cellanno)) +
  geom_point(aes(size = point_size)) + #scale_size_continuous(range = c(0.2, 1)) + 
  theme_classic() + theme(panel.background = element_rect(fill = '#000000FF'))+
  guides(colour = guide_legend(override.aes = list(size=3))) + scale_color_manual(values = vascular_xenium_col) + scale_size_identity()   



