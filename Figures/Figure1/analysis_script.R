



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
all_cell$tissue_type_mer <- factor(all_cell$tissue_type_mer, levels = c("Normal", "DCIS", "synchr", "IDC", "mucinous"))

DimPlot(subset(all_cell, tissue_type_mer == 'mucinous', invert = T), group.by = 'tissue_type_mer', cols = tissue_type_col, raster = F, shuffle = T) #

                                                                                              


                                                                                              
