library(dplyr)
library(clusterProfiler)
library(Seurat)
library(org.Hs.eg.db)


###--------Supplementary figure 2a--------###
DefaultAssay(fibro_all) <- "RNA"
fibro_all_marker <- FindAllMarkers(fibro_all, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden

fibro_all_filtered_marker_top5 <- fibro_all_marker %>% filter(p_val_adj < 0.05) %>% filter(pct.2 < 0.5) %>% group_by(cluster) %>% 
                                    filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object = fibro_all_filtered, features = unique(fibro_all_filtered_marker_top5$gene))+ RotatedAxis() + scale_color_gradientn(colours = rev(paletteer_c("ggthemes::Red-Blue Diverging", 100)))


###--------Supplementary figure 2b--------###
cell_stat_plot <- c('Fibro_FN1', 'Fibro_IGF1')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% 
              rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


###--------Supplementary figure 2c--------###
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

fibroph_genes <- list(CXCL1 = fibro0_gene$ENTREZID, CFD = fibro1_gene$ENTREZID, FN1 = fibro2_gene$ENTREZID, IGF1 = fibro3_gene$ENTREZID, SFRP4 = fibro4_gene$ENTREZID)
fibroph_GOclusterplot <- compareCluster(geneCluster = fibroph_genes, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP")
fibroph_GOclusterplot1 <- simplify(fibroph_GOclusterplot, cutoff=0.5, by="p.adjust", select_fun=min)

dotplot(fibroph_GOclusterplot1, showCategory = 12) + scale_y_discrete(labels=function(x) str_wrap(x, width=85)) + 
  scale_colour_gradientn(colours = (paletteer::paletteer_c("grDevices::Blue-Red 2", 50)))


###--------Supplementary figure 2d--------###
vas_all_filtered_marker_top5 <- vas_all_filtered_marker %>% filter(p_val_adj < 0.05) %>% filter(pct.2 < 0.5) %>% group_by(cluster) %>%
                                  filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object = vas_all_filtered, features = unique(vas_all_filtered_marker_top5$gene))+ RotatedAxis() + scale_color_gradientn(colours = rev(paletteer_c("ggthemes::Red-Blue Diverging", 100)))


###--------Supplementary figure 2e--------###
cell_stat_plot <- c('capillary', 'tec', 'vas-cycling')

supp_pwc <- cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% 
  rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + 
  geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + 
  rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


###--------Supplementary figure 2g--------###
perivas_all_unfiltered <- readRDS("perivas.seurat.inte.proj.all.upd.rds")
DefaultAssay(perivas_all_unfiltered) <- "integrated"

perivas_all <- subset(perivas_all_unfiltered, nFeature_RNA > 500)
DefaultAssay(perivas_all) <- "integrated"
perivas_all <- ScaleData(perivas_all, verbose = FALSE)
perivas_all <- RunPCA(perivas_all, npcs = 30, verbose = FALSE)
perivas_all <- RunUMAP(perivas_all, reduction = "pca", dims = 1:20)
perivas_all <- FindNeighbors(perivas_all, dims = 1:20)
perivas_all <- FindClusters(perivas_all, resolution = 0.2)
                                                                      
DefaultAssay(perivas_all) <- "RNA"
perivas_all_marker <- FindAllMarkers(perivas_all, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden

DimPlot(perivas_all)


###--------Supplementary figure 2h--------###
peri_all_filtered_marker_top5 <- perivas_all_marker %>% filter(p_val_adj < 0.05) %>% filter(pct.2 < 0.5) %>% group_by(cluster) %>% filter(!grepl("^MT-",gene) & !grepl("^RPL",gene) & !grepl("^RPS",gene)) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object = peri_all_filtered, features = unique(peri_all_filtered_marker_top5$gene))+ RotatedAxis() + scale_color_gradientn(colours = rev(paletteer_c("ggthemes::Red-Blue Diverging", 100)))


###--------Supplementary figure 2i--------###
supp_pwc <- cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% 
  rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + 
  geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + 
  rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


                                                                      
