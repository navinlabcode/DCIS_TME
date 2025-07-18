


###--------Supplementary figure 5 macrophage top marker plot--------###
DefaultAssay(macrophage_all) <- "RNA"
macrophage_all_marker <- FindAllMarkers(macrophage_all, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden

macrophage_all_filtered_marker_top5 <- macrophage_all_marker %>% filter(p_val_adj < 0.05) %>% filter(pct.2 < 0.5) %>% group_by(cluster) %>% 
  filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object = macrophage_all_filtered, features = unique(macrophage_all_filtered_marker_top5$gene))+ RotatedAxis() + scale_color_gradientn(colours = rev(paletteer_c("ggthemes::Red-Blue Diverging", 100)))


###--------Supplementary figure 5 macrophage proportion comparison--------###
cell_stat_plot <- c('mye-cycling', 'macro-lipo', 'macro-c3', 'macro-lyve1')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% 
  rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


###--------Supplementary figure 5 macrophage signature analysis--------###
macrophage_all[['signature']] <- CreateAssayObject(data = t(x = FetchData(object = macrophage_all, vars = c('Phago_signature1',"M1_signature1", "M2_signature1"))))
macrophage_all_avg <- AverageExpression(macrophage_all, return.seurat = T, assays = "signature")
macrophage_all_avg$cellstate <- rownames(macrophage_all_avg@meta.data)
macrophage_all_avg$cellstate <- factor(macrophage_all_avg$cellstate, levels = c("Macro_C3", "Macro_LYVE1", "Macro_CCL4","Macro_APOC1",  "Macro_ISG15"))
DoHeatmap(macrophage_all_avg, features = c("M1-signature1", "M2-signature1", "Phago-signature1"),  draw.lines = F, size = 3.6, raster = F, group.by = 'cellstate',  group.colors = mye_col[c(3:7)], slot = 'data')+ 
  scale_fill_gradientn(colours = paletteer::paletteer_c("grDevices::Red-Green", 30))


###--------Supplementary figure 5 macrophage spatial analysis--------###
dcis_300_xenium_srt_12cons_macrophage <- merge(dcis_300_xenium_srt_12cons_basal_rm, dcis_300_xenium_srt_12cons_lumhr, merge.dr = 'SP')
dcis_300_xenium_srt_12cons_macrophage$cellstate_ml_merge <- ifelse(dcis_300_xenium_srt_12cons_macrophage$cellstate_ml_merge %in% c("AMS_others", "others"), 'Basal_others', dcis_300_xenium_srt_12cons_macrophage$cellstate_ml_merge)
sample_check <- unique(dcis_300_xenium_srt_12cons_macrophage$sample)[2]

{
  macrophage_plot <- subset(dcis_300_xenium_srt_12cons_macrophage, sample == sample_check) 
  macrophage_plot <- subset(macrophage_plot, cellstate_ml_merge == 'Normal_like_lumhr', invert = T) 
  macrophage_plot_df <- data.frame(macrophage_plot@reductions$SP@cell.embeddings, macrophage_plot$cellstate_ml_merge)
  colnames(macrophage_plot_df) <- c('X', 'Y', 'cellanno')
  #macrophage_plot_df$point_size <- ifelse(macrophage_plot_df$cellanno %in% c("vein", "artery", "capillary"), 0.4, macrophage_plot_df$point_size)
  ggplot(macrophage_plot_df, aes(x = X, y = Y, color = cellanno)) +
    geom_point(aes(size = point_size)) + #scale_size_continuous(range = c(0.2, 1)) + 
    theme_classic() + theme(panel.background = element_rect(fill = '#000000FF'))+
    guides(colour = guide_legend(override.aes = list(size=3))) + scale_color_manual(values = xenium_basal_col) + scale_size_identity()   
}
















