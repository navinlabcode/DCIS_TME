
basal_all <- readRDS("")

# figure 3d enrichment analysis -----------------------------------------------------
basal_signature <- basal_all_marker %>% group_by(cluster) %>% filter(avg_log2FC > 0.5 & pct.2 < 0.7 & p_val_adj < 0.05) %>% top_n(n = 50, wt = avg_log2FC)
basal_df <- basal_signature[,6:7]
basal_dfsample <- split(basal_df$gene, basal_df$cluster)
{
  basal1_gene <- bitr(as.character(basal_dfsample$`Basal-major`), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  basal2_gene <- bitr(as.character(basal_dfsample$`Basal-interferon`), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  basal3_gene <- bitr(as.character(basal_dfsample$CAB), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
}
basalph_genes <- list(major = basal1_gene$ENTREZID, int = basal2_gene$ENTREZID, cab = basal3_gene$ENTREZID)
basalph_GOclusterplot <- compareCluster(geneCluster = basalph_genes, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP")
dotplot(basalph_GOclusterplot, showCategory = 8) + scale_y_discrete(labels=function(x) str_wrap(x, width=85)) +
  scale_colour_gradientn(colours = (paletteer::paletteer_c("grDevices::Blue-Red 2", 50)))

basalph_GOclusterplot1 <- simplify(basalph_GOclusterplot, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(basalph_GOclusterplot1, showCategory = 7) + scale_y_discrete(labels=function(x) str_wrap(x, width=80)) + 
  scale_colour_gradientn(colours = (paletteer::paletteer_c("grDevices::Blue-Red 2", 50)))

# figure 3e integrin markers -----------------------------------------------------
VlnPlot(basal_all, c('ITGB1', 'ITGA2', 'ITGA6', 'ITGA3', 'ITGB4', 'ITGBL1', 'ITGAV', 'ITGB6'), pt.size = 0, cols = basal_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()

