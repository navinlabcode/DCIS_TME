
er_tumor_obj <- readRDS("./er_tumor_obj.rds")

###--------Figure 6b--------###
DotPlot(object = er_tumor_obj, features = c('TMC5', hbca_epi_final_use_sub_marker_top9$gene[!hbca_epi_final_use_sub_marker_top9$gene == 'TCIM']), group.by = "orig.ident", dot.scale = 3, dot.min = 0.01)+ RotatedAxis() + 
  scale_color_gradientn(colours = viridis::viridis(100)) & coord_flip()

###--------Figure 6c--------###
normal_epi_hvg_24 <- read.table('/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/tumor/cell_origin/diffusion_map/normal_epi_hvg_24.txt', header = T)
hbca_epi_count_dm_24 <- readRDS("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/tumor/cell_origin/diffusion_map/normal_epi_sub_count_dm.rds")
hbca_epi_anno_24 <- read.table("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/tumor/cell_origin/diffusion_map/normal_epi_sub_anno.txt", header = T)
hbca_epi_dm_tmp_24 <- data.frame(DC1 = eigenvectors(hbca_epi_count_dm_24)[,1], DC2 = eigenvectors(hbca_epi_count_dm_24)[,2], CellType = hbca_epi_anno_24$x)

er_tumor_obj_diet <- readRDS('/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/tumor/ER_pos_tumor_diet_0924.rds')
DefaultAssay(er_tumor_obj_diet) <- 'RNA'
er_tumor_obj_diet <- NormalizeData(er_tumor_obj_diet) 

tumor_df <- as.data.frame(er_tumor_obj_diet@assays$RNA@data)
tumor_count <- tumor_df %>% filter(rownames(tumor_df) %in% normal_epi_hvg_24$x)
tumor_count_df <- as.data.frame(t(tumor_count))
tumor_count_dm_pred <- dm_predict(hbca_epi_count_dm_24, tumor_count_df)

tumor_anno <- data.frame(anno = er_tumor_obj_diet$tissue_upd_0924)
tumor_count_df$id <- 1:nrow(tumor_count_df)
tumor_count_df_anno <- merge(tumor_count_df, tumor_anno, by = 0)
tumor_count_df_anno <- tumor_count_df_anno[order(tumor_count_df_anno$id), ]

tumor_count_dm_pred_tmp <- data.frame(DC1 = tumor_count_dm_pred[,1], DC2 = tumor_count_dm_pred[,2], CellType = paste0("tumor_", tumor_count_df_anno$anno))

tumor_dif_pred <- rbind(hbca_epi_dm_tmp_24, tumor_count_dm_pred_tmp)

ggplot()+geom_point(data = tumor_dif_pred %>% filter(!CellType == "tumor_synch_yes"), aes(x = DC1, y = DC2, colour = CellType), size = 0.4)+
  geom_point(data = tumor_dif_pred %>% filter(CellType == "tumor_synch_yes"), aes(x = DC1, y = DC2, colour = CellType), size = 0.4)+ggpubr::theme_pubr()+
  theme(axis.title = element_text(size = 16), axis.text=ggplot2::element_text(size=12), plot.margin = ggplot2::unit(c(0.5,1,0.5,1), "cm"), legend.title = element_text(size = 14),legend.text = element_text(size = 12))+
  guides(colour = guide_legend(override.aes = list(size=3))) + scale_colour_manual(values = c('#FC7EADFF', "#41B6E6FF", '#31D64DFF', "#F9B90AFF",  "#F24C3DFF", '#9063CDFF', '#2CCCD3FF')) 
