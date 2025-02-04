
er_tumor_obj <- readRDS("./er_tumor_obj.rds")

###--------Figure 6b--------###
DotPlot(object = er_tumor_obj, features = c('TMC5', hbca_epi_final_use_sub_marker_top9$gene[!hbca_epi_final_use_sub_marker_top9$gene == 'TCIM']), group.by = "orig.ident", dot.scale = 3, dot.min = 0.01)+ RotatedAxis() + 
  scale_color_gradientn(colours = viridis::viridis(100)) & coord_flip()

###--------Figure 6c--------###
normal_epi_hvg_24 <- read.table('./diffusion_map/normal_epi_hvg_24.txt', header = T)
hbca_epi_count_dm_24 <- readRDS("./normal_epi_sub_count_dm.rds")
hbca_epi_anno_24 <- read.table("./diffusion_map/normal_epi_sub_anno.txt", header = T)
hbca_epi_dm_tmp_24 <- data.frame(DC1 = eigenvectors(hbca_epi_count_dm_24)[,1], DC2 = eigenvectors(hbca_epi_count_dm_24)[,2], CellType = hbca_epi_anno_24$x)

er_tumor_obj_diet <- er_tumor_obj
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

###--------Figure 6e--------###
nmf_programs_genes_tumor <- readRDS('tumor_w_list_all95.rds')
tissue_type_upd_0709$sample_upd <- sub('_pc1|_24hTis|_cryo|_24htis', '', tissue_type_upd_0709$sample)
nmf_programs_genes_tumor <- nmf_programs_genes_tumor[intersect(tissue_type_upd_0709$sample_upd, names(nmf_programs_genes_tumor))] 

tumor_nmf_programs_top50_ori <- lapply(nmf_programs_genes_tumor, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
tumor_nmf_filter <- robust_nmf_programs(tumor_nmf_programs_top50_ori, intra_min=35, intra_max=10, inter_filter=T, inter_min=10)
tumor_nmf_programs_top50 <- lapply(tumor_nmf_programs_top50_ori, function(x) x[, is.element(colnames(x), tumor_nmf_filter), drop=F])
tumor_nmf_programs_top50 <- do.call(cbind, tumor_nmf_programs_top50)

# calculate similarity between NMF programs
nmf_programs_sig_tumor_rshp <- tumor_nmf_programs_top50 %>% melt %>% set_colnames(c('X1', 'X2', 'gene')) %>% dplyr::select(-X1) %>% mutate(exst=1) %>% 
  tidyr::pivot_wider(id_cols=X2, names_from=gene, values_from=exst, values_fill=0)
nmf_programs_sig_tumor_rshp_mat <- as.matrix(nmf_programs_sig_tumor_rshp[, -1])
tumor_nmf_intersect <- nmf_programs_sig_tumor_rshp_mat %*% t(nmf_programs_sig_tumor_rshp_mat) %>% set_colnames(nmf_programs_sig_tumor_rshp$X2) %>% set_rownames(nmf_programs_sig_tumor_rshp$X2)
tumor_nmf_intersect_rmrg <- remerge(as.dist((50-tumor_nmf_intersect)/50), intra_thr=0.8, inter_thr=.72, min_size=4, method='ward.D2', init_k=23, max_iter=20) #intra_thr smaller, between cluster are less similar
tumor_nmf_intersect_rmrg_df <- as.matrix(tumor_nmf_intersect_rmrg$dist_out)
tumor_nmf_intersect_rmrg_group <- data.frame(grp=tumor_nmf_intersect_rmrg$grp_out) %>% set_rownames(names(tumor_nmf_intersect_rmrg$grp_out))

rec_break <- cumsum(table(tumor_nmf_intersect_rmrg_group$grp)[unique(tumor_nmf_intersect_rmrg_group$grp)])
                                   
tumor_nmf_intersect_rmrg_df_reord_melt <- reshape2::melt(tumor_nmf_intersect_rmrg_df)
ggplot(data = tumor_nmf_intersect_rmrg_df_reord_melt, aes(x=Var1, y=Var2, fill=1-value)) + 
    geom_tile() +
    scale_fill_gradientn(limits=c(0, 1), colors = c("white", rev(paletteer::paletteer_c("grDevices::Reds 3", n = 100))), name="Similarity score") +
    scale_x_discrete(name="\nPrograms", breaks=unique(tumor_nmf_intersect_rmrg_df_reord_melt$Var1)[seq(50, length(unique(tumor_nmf_intersect_rmrg_df_reord_melt$Var1)), by=50)], labels= seq(50, length(unique(tumor_nmf_intersect_rmrg_df_reord_melt$Var1)), by=50)) + 
    scale_y_discrete(name="\nPrograms", breaks=unique(tumor_nmf_intersect_rmrg_df_reord_melt$Var2)[seq(50, length(unique(tumor_nmf_intersect_rmrg_df_reord_melt$Var2)), by=50)], labels= seq(50, length(unique(tumor_nmf_intersect_rmrg_df_reord_melt$Var2)), by=50)) +
    geom_rect(aes(xmin = 0.5, xmax = rec_break[1] + 0.5, ymin = 0.5, ymax = rec_break[1] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[1] + 0.5, xmax = rec_break[2] + 0.5, ymin = rec_break[1] + 0.5, ymax = rec_break[2] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[2] + 0.5, xmax = rec_break[3] + 0.5, ymin = rec_break[2] + 0.5, ymax = rec_break[3] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[3] + 0.5, xmax = rec_break[4] + 0.5, ymin = rec_break[3] + 0.5, ymax = rec_break[4] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[4] + 0.5, xmax = rec_break[5] + 0.5, ymin = rec_break[4] + 0.5, ymax = rec_break[5] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[5] + 0.5, xmax = rec_break[6] + 0.5, ymin = rec_break[5] + 0.5, ymax = rec_break[6] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[6] + 0.5, xmax = rec_break[7] + 0.5, ymin = rec_break[6] + 0.5, ymax = rec_break[7] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[7] + 0.5, xmax = rec_break[8] + 0.5, ymin = rec_break[7] + 0.5, ymax = rec_break[8] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[8] + 0.5, xmax = rec_break[9] + 0.5, ymin = rec_break[8] + 0.5, ymax = rec_break[9] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[9] + 0.5, xmax = rec_break[10] + 0.5, ymin = rec_break[9] + 0.5, ymax = rec_break[10] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[10] + 0.5, xmax = rec_break[11] + 0.5, ymin = rec_break[10] + 0.5, ymax = rec_break[11] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[11] + 0.5, xmax = rec_break[12] + 0.5, ymin = rec_break[11] + 0.5, ymax = rec_break[12] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[12] + 0.5, xmax = rec_break[13] + 0.5, ymin = rec_break[12] + 0.5, ymax = rec_break[13] + 0.5), fill = NA, color = "black") +
    geom_rect(aes(xmin = rec_break[13] + 0.5, xmax = rec_break[14] + 0.5, ymin = rec_break[13] + 0.5, ymax = rec_break[14] + 0.5), fill = NA, color = "black") 

###--------Figure 6g--------###
normal_epi_hr <- readRDS(./normal_epi_hr.rds)
hbca_tumor_obj_com <- merge(er_tumor_obj_diet, normal_epi_hr)
hbca_tumor_obj_com@meta.data[, names(module_gene_top50)] <- NULL
hbca_tumor_obj_com <- AddModuleScore(hbca_tumor_obj_com, features = module_gene_top50)
colnames(hbca_tumor_obj_com@meta.data)[grepl('Cluster', colnames(hbca_tumor_obj_com@meta.data))] <- names(module_gene_top50)

tumor_obj_dat <- data.frame(sample = hbca_tumor_obj_com$orig.ident, hbca_tumor_obj_com@meta.data[, names(module_gene_top50)]) #tumor_obj$type
tumor_obj_dat_type <- merge(tumor_obj_dat, dcis_clinical_0924[, c(1,5)], by = "sample")

###binarize to proportion
tumor_obj_dat_type_trunc <- data.frame(sample = tumor_obj_dat_type$sample, type = tumor_obj_dat_type$tissue_upd_0924, ifelse(tumor_obj_dat_type[, names(module_gene_top50)] > 0.25, 1, 0)) #test more cutoff
tumor_obj_dat_type_trunc_sum_upd <- tumor_obj_dat_type_trunc %>%
    group_by(sample, type) %>%
    summarize(
      across(1:14, ~sum(. == 1) / n(), .names = "{.col}")
    )
tumor_obj_dat_type_trunc_sum_ref_upd <- tumor_obj_dat_type_trunc_sum_upd %>% pivot_longer(!c(sample, type), names_to = "Module", values_to = "Proportion")
tumor_obj_dat_type_trunc_sum_ref_upd$Module <- factor(tumor_obj_dat_type_trunc_sum_ref_upd$Module, levels = names(module_gene_top50))
tumor_obj_dat_type_trunc_sum_ref_upd$type_upd <- tumor_obj_dat_type_trunc_sum_ref_upd$type
tumor_obj_dat_type_trunc_sum_ref_upd$type_upd <- factor(tumor_obj_dat_type_trunc_sum_ref_upd$type_upd, levels = c("Normal" , "DCIS_yes", "synch_yes", "IDC_yes"))
facet(ggboxplot(tumor_obj_dat_type_trunc_sum_ref_upd, x = "type_upd", y = "Proportion", color = "type_upd", palette = tissue_type_col, add = "jitter")+
          stat_compare_means(aes(group = type_upd), label = "p.format" , method = 'kruskal.test'), facet.by = "Module", nrow = 3, scales = 'free_y')+theme(axis.text.x=element_text(angle=50,hjust=1,vjust=1))




