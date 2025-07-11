library(Seurat)
library(ggplot2)
library(dplyr)

dcis_300_xenium_srt_9cons_fibroblast_basal <- merge(dcis_300_xenium_srt_9cons_basal_rm, c(dcis_300_xenium_srt_9cons_fibroblast, dcis_300_xenium_srt_9cons_lumhr), merge.dr = 'SP')

{
  fibro_plot <- subset(dcis_300_xenium_srt_9cons_fibroblast_basal, cellstate_ml_com %in% c("basal", "fibro_sfrp4", "fibro_cfd",  "fibro_fn1",  "fibro_igf1",       
                                                                                            "fibro_cycling", "Tumor_lumhr") & sample == sample_check)  
  fibro_plot_df <- data.frame(fibro_plot@reductions$SP@cell.embeddings, fibro_plot$cellstate_ml_com)
  colnames(fibro_plot_df) <- c('X', 'Y', 'cellanno')
  fibro_plot_df$point_size <- ifelse(fibro_plot_df$cellanno %in% c("basal", "Tumor_lumhr"), 0.06, 0.6)
  ggplot(fibro_plot_df, aes(x = X, y = Y, color = cellanno)) +
    geom_point(aes(size = point_size)) + #scale_size_continuous(range = c(0.2, 1)) + 
    theme_classic() + theme(panel.background = element_rect(fill = '#000000FF'))+
    guides(colour = guide_legend(override.aes = list(size=3))) + scale_color_manual(values = fibro_xenium_col) + scale_size_identity()   
}
