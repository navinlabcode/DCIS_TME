
fibro_all_filtered_marker_top5 <- fibro_all_filtered_marker %>% filter(p_val_adj < 0.05) %>% filter(pct.2 < 0.5) %>% group_by(cluster) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=5, wt=avg_log2FC)
DotPlot(object = fibro_all_filtered, features = unique(fibro_all_filtered_marker_top5$gene))+ RotatedAxis() + scale_color_gradientn(colours = rev(paletteer_c("ggthemes::Red-Blue Diverging", 100)))



# Supplementary Fig2c fibroblast signature analysis ---------------------------------------------------------------
fibro_signature <- fibro_all_marker %>% group_by(cluster) %>% filter(avg_log2FC > 0.5 & pct.2 < 0.7 & p_val_adj < 0.05) %>% top_n(n = 50, wt = avg_log2FC)
fibro_df <- fibro_signature[,6:7]

fibro_dfsample <- split(fibro_df$gene, fibro_df$cluster)
{
  fibro0_gene <- bitr(as.character(fibro_dfsample$Fibro_CXCL1), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  fibro1_gene <- bitr(as.character(fibro_dfsample$Fibro_CFD), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  fibro2_gene <- bitr(as.character(fibro_dfsample$Fibro_FN1), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  fibro3_gene <- bitr(as.character(fibro_dfsample$Fibro_IGF1), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  fibro4_gene <- bitr(as.character(fibro_dfsample$Fibro_SFRP4), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
}

fibroph_genes <- list(CXCL1 = fibro0_gene$ENTREZID, CFD = fibro1_gene$ENTREZID, FN1 = fibro2_gene$ENTREZID, IGF1 = fibro3_gene$ENTREZID, SFRP4 = fibro4_gene$ENTREZID) #, stress = fibro5_gene$ENTREZID, cycling = fibro6_gene$ENTREZID)
fibroph_GOclusterplot <- compareCluster(geneCluster = fibroph_genes, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP")
fibroph_GOclusterplot1 <- simplify(fibroph_GOclusterplot, cutoff=0.5, by="p.adjust", select_fun=min)

dotplot(fibroph_GOclusterplot1, showCategory = 12) + scale_y_discrete(labels=function(x) str_wrap(x, width=85)) + 
  scale_colour_gradientn(colours = (paletteer::paletteer_c("grDevices::Blue-Red 2", 50)))
