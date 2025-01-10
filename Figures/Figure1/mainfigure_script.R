## Note: the analysis below is based on the integrated object ##

all_cell <- readRDS("~/all_cell.rds")

###--------Figure 1b--------###
DimPlot(all_cell, group.by = "celltype", cols = celltype_col, label = T, raster = F)

###--------Figure 1c--------###
Idents(all_cell) <- all_cell$celltype
DefaultAssay(all_cell) <- "RNA"
hbca_dcis_ibc.combined_marker <- FindAllMarkers(all_cell, logfc.threshold = 0.3, only.pos = T, slot = "data")

normal_dcis_idc_deg_top8 <- normal_dcis_idc_deg %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% filter(pct.1 > 0.5) %>% filter(pct.2 < 0.5) %>% top_n(n=7, wt=avg_log2FC) #%>% filter(pct.1 > 0.8)
normal_dcis_idc_deg_top8 <- normal_dcis_idc_deg_top8[with(normal_dcis_idc_deg_top8, order(ave(seq_along(cluster), cluster,  FUN = function(x) match(cluster[x], c("Basal", "lumhr", "lumsec",  "Bcell", "NK_T", "Myeloid", "Fibroblast", "Vascular", "Lymphatic", "PeriVas"))))), ]
DotPlot(object = all_cell, group.by = 'celltype', features = unique(normal_dcis_idc_deg_top8$gene))+ RotatedAxis() + scale_color_gradientn(colours = rev(paletteer_c("ggthemes::Red-Blue Diverging", 100)))

###--------Figure 1d--------###
DimPlot(all_cell, group.by = 'tumor_normal', cols = c('#DC863BFF', '#a4ceb7'), raster = F)

###--------Figure 1e--------###                                                                                          
all_cell$tissue_type_mer <- all_cell$tissue_type
all_cell$tissue_type_mer <- ifelse(all_cell$tissue_type_mer %in% c('DCIS_no', 'DCIS_yes'), 'DCIS', all_cell$tissue_type_mer)
all_cell$tissue_type_mer <- ifelse(all_cell$tissue_type_mer %in% c('synch_yes', 'synch_no'), 'synchr', all_cell$tissue_type_mer)
all_cell$tissue_type_mer <- ifelse(all_cell$tissue_type_mer %in% c('IDC_no', 'IDC_yes'), 'IDC', all_cell$tissue_type_mer)
all_cell$tissue_type_mer <- factor(all_cell$tissue_type_mer, levels = c("Normal", "DCIS", "synchr", "IDC"))

DimPlot(all_cell, group.by = 'tissue_type_mer', cols = tissue_type_col, raster = F, shuffle = T)
                                                                                              
###--------Figure 1f--------### 
hbca_dcis_ibc.integrated.all.meta <- data.frame(sample = all_cell$orig.ident, proj_ori = all_cell$project,
                                                celltype = all_cell$celltype_com, ploidy = all_cell$tumor_normal)
hbca_dcis_ibc.integrated.all.meta.tis <- merge(hbca_dcis_ibc.integrated.all.meta, dcis_clinical, by = 'sample', all.x = T)
all_tme_cell_count2_cli_upd <- data.frame(table(hbca_dcis_ibc.integrated.all.meta.tis$tissue_upd_0924, hbca_dcis_ibc.integrated.all.meta.tis$celltype))
allcells_kept_sum1_proj_upd_mel_sel <- all_tme_cell_count2_cli_upd %>% filter(Var1 %in% c("Normal", "DCIS_yes", "synch_yes", "IDC_yes"))

ggplot(allcells_kept_sum1_proj_upd_mel_sel %>% filter(!Var2 %in% c('lumhr', 'lumsec')), aes(fill=Var2, y=Freq, x=Var1)) + ylab("Proportion")+
  geom_bar(position="fill", stat="identity", colour = "black")+ggpubr::theme_pubr(legend = "right")+
  scale_fill_manual(values = c(celltype_col[1], celltype_col[4:10]))+font("y.title", size = 14)+
  scale_x_discrete(limits = c("Normal", "DCIS_yes", "synch_yes", "IDC_yes"))+theme(axis.text.x=element_text(angle=50,hjust=1,vjust=1))
                                                                                              
###--------Figure 1g--------### 
all_tme_sample_cell_count <- data.frame(table(hbca_dcis_ibc.integrated.all.meta$sample, hbca_dcis_ibc.integrated.all.meta$celltype))
all_tme_sample_cell_count_tis <- merge(all_tme_sample_cell_count, dcis_clinical, by.x = 'Var1', by.y = 'sample')
all_tme_sample_cell_count_tis_sel <- all_tme_sample_cell_count_tis %>% filter(tissue_upd_0924 %in% c("Normal", "DCIS_yes", "synch_yes", "IDC_yes"))
all_tme_sample_cell_count_tis_sel_prop <- all_tme_sample_cell_count_tis_sel %>% dplyr::filter(!Var2 %in% c('lumhr', 'lumsec')) %>% group_by(Var1) %>% mutate(Prop = Freq/sum(Freq))
all_tme_sample_cell_count_tis_sel_prop$tissue_upd_0924 <- factor(all_tme_sample_cell_count_tis_sel_prop$tissue_upd_0924, levels = c("Normal", "DCIS_yes",  "synch_yes", "IDC_yes"))

pwc_celltype <- all_tme_sample_cell_count_tis_sel_prop %>% filter(!Var2 %in% c('Vascular',  'Lymphatic', 'Bcell', 'PeriVas')) %>% group_by(Var2) %>% 
  rstatix::dunn_test(Prop ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")
ggplot(all_tme_sample_cell_count_tis_sel_prop %>% filter(!Var2 %in% c('Vascular', 'Lymphatic', 'Bcell', 'PeriVas')), aes(x = tissue_upd_0924, y = Prop)) + geom_boxplot(outliers = F) + geom_jitter(aes(colour = tissue_upd_0924), size = 0.9) +
  facet_wrap(facets = vars(Var2), nrow = 1) + theme_bw()+rremove("x.ticks") + rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + stat_pvalue_manual(pwc_celltype, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")

                                                                                              
