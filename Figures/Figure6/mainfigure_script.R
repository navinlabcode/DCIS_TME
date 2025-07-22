
er_tumor_obj <- readRDS("./er_tumor_obj.rds")

###--------Figure 6b--------###
DimPlot(er_tumor_obj, group.by = "patient")


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


###--------Figure 6d--------###
normal_tumor_type1 <- dcis_clinical_0924[c(1,5)]
colnames(normal_tumor_type1) <- c('sample', 'tumor_path')
hbca_tumor_obj_com$tissue_type <- NULL
{
  hbca_tumor_obj_com$tissue_type <- ifelse(hbca_tumor_obj_com$orig.ident %in% (normal_tumor_type1 %>% filter(tumor_path == "Normal") %>% pull(sample)), "Normal", "na" )
  hbca_tumor_obj_com$tissue_type <- ifelse(hbca_tumor_obj_com$orig.ident %in% (normal_tumor_type1 %>% filter(tumor_path %in% c("DCIS_no", "DCIS_yes")) %>% pull(sample)), "DCIS", hbca_tumor_obj_com$tissue_type)
  hbca_tumor_obj_com$tissue_type <- ifelse(hbca_tumor_obj_com$orig.ident %in% (normal_tumor_type1 %>% filter(tumor_path %in% c("synch_yes", "synch_no")) %>% pull(sample)), "synch", hbca_tumor_obj_com$tissue_type)
  hbca_tumor_obj_com$tissue_type <- ifelse(hbca_tumor_obj_com$orig.ident %in% (normal_tumor_type1 %>% filter(tumor_path %in% c("IDC_no", "IDC_yes")) %>% pull(sample)), "IDC", hbca_tumor_obj_com$tissue_type)
}
hbca_tumor_obj_com$tissue_type <- factor(hbca_tumor_obj_com$tissue_type, levels = c("Normal", "DCIS", "synch", "IDC"))

hbca_tumor_obj_com_diet <- hbca_tumor_obj_com


##DE genes
Idents(hbca_tumor_obj_com) <- hbca_tumor_obj_com$tissue_type
normal_dcis_idc_upd_marker_0924 <- FindAllMarkers(hbca_tumor_obj_com, logfc.threshold = 0.25, only.pos = T, slot = "data") #for cluster number iden

write.csv(normal_dcis_idc_upd_marker_0924, './tumor/normal_dcis_idc_upd_marker_0924.csv', quote = F)

normal_dcis_idc_upd_marker_0924_top30 <- normal_dcis_idc_upd_marker_0924 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

#combine cosmic genes
normal_dcis_idc_upd_marker_0924_top20 <- normal_dcis_idc_upd_marker_0924 %>% group_by(cluster) %>% filter(avg_log2FC > 0.25) %>% top_n(n = 20, wt = avg_log2FC) 

normal_dcis_idc_upd_maker_cancergene <- normal_dcis_idc_upd_marker_0924 %>% group_by(cluster) %>% filter(avg_log2FC > 0.25) %>% filter(gene %in% cosmic_updated$Gene) %>% top_n(n = 20, wt = avg_log2FC) 
normal_dcis_idc_upd_maker_cancergene_info <- merge(normal_dcis_idc_upd_maker_cancergene, cosmic_updated, by.x = 'gene', by.y = 'Gene')

hbca_tumor_obj_com$tissue_type1 <- factor(hbca_tumor_obj_com$tissue_type, levels = c("IDC", "synch", "DCIS", "Normal"))

DotPlot(hbca_tumor_obj_com, features = c('PGR', 'AR', 'SMAD3', 'SBDS', 'DICER1', 'CHD2', 'KLF6', 'CDH1', 'ERBB2', 'LASP1', 'GATA3', 'CD74', 'KLF4', 'ESR1'),  group.by = "tissue_type1", dot.min = 0.01)+ RotatedAxis() + scale_color_gradientn(colours = c(paletteer::paletteer_c("grDevices::Viridis", n = 100))) +
  theme(axis.text.x=element_text(size=10))


##Known genes
cosmic_updated <- read.table('./cosmic_updated.txt', fill = TRUE, header = T)
normal_dcis_idc_upd_maker_cancergene <- normal_dcis_idc_upd_marker %>% filter(avg_log2FC > 0.3) %>% filter(gene %in% cosmic_updated$Gene)
normal_dcis_idc_upd_maker_cancergene <- normal_dcis_idc_upd_marker %>% group_by(cluster) %>% filter(avg_log2FC > 0.25) %>% filter(gene %in% cosmic_updated$Gene) %>% top_n(n = 8, wt = avg_log2FC) 
genes_to_check <- normal_dcis_idc_upd_maker_cancergene %>% pull(gene) #c('GATA3', 'DNAJB1', 'CCND1')
genes_to_check <- c(genes_to_check[c(1:5)], 'CDH1', genes_to_check[c(6:9)] , 'KLF4')

breast_gene <- read.table('./breast_cancer_uniq_sorted.txt', header = F) 
normal_dcis_idc_upd_maker_breastgene <- normal_dcis_idc_upd_marker %>% filter(avg_log2FC > 0.3) %>% filter(gene %in% breast_gene$V1)
genes_to_check <- normal_dcis_idc_upd_maker_breastgene %>% pull(gene) #c('GATA3', 'DNAJB1', 'CCND1')

##plot genes
hbca_tumor_obj_com_diet$tissue_type <- factor(hbca_tumor_obj_com_diet$tissue_type, levels = c("IDC", "DCIS", "Normal"))

normal_dcis_idc_upd_topmarker <- normal_dcis_idc_upd_marker %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) 
normal_dcis_idc_upd_topmarker_pct <- normal_dcis_idc_upd_marker %>% group_by(cluster) %>% filter(pct.1 > 0.3) %>% top_n(n = 3, wt = avg_log2FC) 
DotPlot(hbca_tumor_obj_com_diet, features = normal_dcis_idc_upd_topmarker$gene,  group.by = "tissue_type", dot.min = 0.01)+ RotatedAxis() + scale_color_gradientn(colours = c(paletteer::paletteer_c("grDevices::Viridis", n = 100))) +
  theme(axis.text.x=element_text(size=10))

normal_dcis_idc_upd_maker_cosmic <- left_join(normal_dcis_idc_upd_marker, cosmic_updated, by = c("gene" = "Gene"))
cancer_gene <- c('SMAD3', 'CDH1', 'SBDS', 'DICER1', 'CHD2', 'KLF6', 'COX6C', 'CD74', 'CCND1', 'GATA3', 'KLF4')

genes_to_check <- c(cancer_gene, normal_dcis_idc_upd_topmarker$gene)
genes_to_check <- c(cancer_gene, normal_dcis_idc_upd_topmarker_pct$gene)

genes_to_check <- c('ESR1', 'PGR', 'AR', 'ERBB2', 'MKI67', 'SMAD3', 'CDH1', 'SBDS', 'DICER1', 'CHD2', 'KLF6', 'COX6C', 'CD74', 'CCND1', 'GATA3', 'KLF4') #previous version

DotPlot(hbca_tumor_obj_com_diet, features = unique(genes_to_check),  group.by = "tissue_type")+ RotatedAxis() + scale_color_gradientn(colours = c(paletteer::paletteer_c("grDevices::Viridis", n = 100))) +
  theme(axis.text.x=element_text(size=10))

DotPlot(hbca_tumor_obj_com_diet, features = unique(c('ESR1', 'PGR', 'AR', 'ERBB2', 'ITGAV', 'SMAD3', 'SFPQ', 'MYH9', 'CDH1', 'GATA3', 'DNAJB1', 'CCND1', 'KLF4')),  group.by = "tissue_type")+ RotatedAxis() + scale_color_gradientn(colours = c(paletteer::paletteer_c("grDevices::Viridis", n = 100))) +
  theme(axis.text.x=element_text(size=10))
DotPlot(hbca_tumor_obj_com_diet, features = unique(c('ESR1', 'PGR', 'AR', 'ERBB2')),  group.by = "tissue_type", scale = F)+ RotatedAxis() + scale_color_gradientn(colours = c(paletteer::paletteer_c("grDevices::Viridis", n = 100))) +
  theme(axis.text.x=element_text(size=10))
VlnPlot(hbca_tumor_obj_com_diet, features = c('ESR1', 'PGR', 'AR', 'ERBB2'),  group.by = "tissue_type", pt.size = 0)+
  theme(axis.text.x=element_text(size=10))
VlnPlot(hbca_tumor_obj_com_diet, c('ESR1', 'PGR', 'AR', 'ERBB2'), pt.size = 0, cols = tissue_type_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()

VlnPlot(hbca_tumor_obj_com_diet, features = 'PALLD', pt.size = 0,  group.by = "tissue_type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(comparisons = list( c("Normal", "DCIS"), c("DCIS", "IDC"), c("Normal", "IDC"))) + ylim(-2, 10)

write.csv(normal_dcis_idc_upd_marker, './tumor/normal_dcis_idc_upd_marker.csv', quote = F)



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


###--------Figure 6 spatial--------###
#normal vs DCIS lumhr proportion comp
dcis_lumhr_count <- data.frame(table(dcis_300_xenium_srt_12cons_lumhr$sample, dcis_300_xenium_srt_12cons_lumhr$cellstate_ml_com))
dcis_lumhr_count$type <- 'DCIS'
dcis_lumhr_count_prop <- dcis_lumhr_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))

hbca_lumhr_count <- data.frame(table(hbca_300_xenium_srt_6cons_lumhr$sample, hbca_300_xenium_srt_6cons_lumhr$cellstate_ml_com))
hbca_lumhr_count$type <- 'HBCA'
hbca_lumhr_count_prop <- hbca_lumhr_count %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))

hbca_dcis_lumhr_prop <- rbind(dcis_lumhr_count_prop, hbca_lumhr_count_prop)
hbca_dcis_lumhr_prop$type <- factor(hbca_dcis_lumhr_prop$type, levels = c('HBCA', 'DCIS'))

ggboxplot(hbca_dcis_lumhr_prop %>% filter(!Var2 %in% c('AMS_others', 'others', 'lumhr_others', 'Normal_like_lumhr')), x = "Var2", y = "Prop", color = "type", palette = tissue_type_col,#x = "supp", y = "len", fill = "#00AFBB", 
          add = "jitter") + stat_compare_means(aes(group = type), label = "p.format", method = "wilcox.test") + rotate_x_text(angle = 45, hjust = 1, vjust = 1) + scale_x_discrete(limits = c('lumhr_cycling', 'lumhr_intef', 'lumhr_hla'))


