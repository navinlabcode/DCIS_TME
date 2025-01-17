myeloid_all_filtered <- subset(all_cell, celltype == "myeloid")

###--------Figure 4a--------###
DimPlot(myeloid_all_filtered, label = T, label.size = 3.5, cols = mye_col, group.by = "cellstate")

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
