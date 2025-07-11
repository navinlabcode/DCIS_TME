library(Seurat)
library(dplyr)


# fig2a fibro UMAP-------------------------------------------------------------------
fibro_all_unfiltered <- readRDS("fibro.seurat.inte.proj.all.upd.rds")
DefaultAssay(fibro_all_unfiltered) <- "integrated"

fibro_all <- subset(fibro_all_unfiltered, nFeature_RNA > 500)
DefaultAssay(fibro_all) <- "integrated"
fibro_all <- ScaleData(fibro_all, verbose = FALSE)
fibro_all <- RunPCA(fibro_all, npcs = 30, verbose = FALSE)
fibro_all <- RunUMAP(fibro_all, reduction = "pca", dims = 1:20)
fibro_all <- FindNeighbors(fibro_all, dims = 1:20)
fibro_all <- FindClusters(fibro_all, resolution = 0.2)
                                                                      
DefaultAssay(fibro_all) <- "RNA"
fibro_all_marker <- FindAllMarkers(fibro_all, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden

DimPlot(fibro_all)

# fig2b fibroblast cell state comparison---------------------------------------------
supp_pwc <- cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")

ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


# fig2c fibro markers-------------------------------------------------------------------
VlnPlot(fibro_all_diet, c('COL15A1', 'IGF1', 'DLK1', 'IGFBP2', 'ACTA2', 'FAP', 'TAGLN', 'POSTN'), pt.size = 0, cols = fibro_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()

fibro_mmp_marker <- fibro_all_marker %>% filter(grepl("^MMP", gene)) %>% pull(gene)
VlnPlot(fibro_all_diet, fibro_mmp_marker, pt.size = 0, cols = fibro_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()
DotPlot(fibro_all_diet, features = unique(fibro_mmp_marker)) + RotatedAxis() + scale_color_gradientn(colours = viridis::viridis(20)) & coord_flip()





# vascular analysis -----------------------------------------------------------
# fig2g vascular state umap-------------------------------------------------------------------
vas_all_unfiltered <- readRDS("vas.seurat.inte.proj.all.upd.rds")
DefaultAssay(vas_all_unfiltered) <- "integrated"

vas_all <- subset(vas_all_unfiltered, nFeature_RNA > 500)
DefaultAssay(vas_all) <- "integrated"
vas_all <- ScaleData(vas_all, verbose = FALSE)
vas_all <- RunPCA(vas_all, npcs = 30, verbose = FALSE)
vas_all <- RunUMAP(vas_all, reduction = "pca", dims = 1:20)
vas_all <- FindNeighbors(vas_all, dims = 1:20)
vas_all <- FindClusters(vas_all, resolution = 0.2)
                                                                      
DefaultAssay(vas_all) <- "RNA"
vas_all_marker <- FindAllMarkers(vas_all, logfc.threshold = 0.3, only.pos = T, slot = "data") #for cluster number iden

DimPlot(vas_all)

# fig2i vascular state umap-------------------------------------------------------------------
VlnPlot(vas_all, c('APLNR', 'INSR', 'FLT1', 'KDR', 'NRP1', 'VWA1', 'COL4A1', 'COL4A2', 'COL6A2','LAMA4','LAMB1'), 
        pt.size = 0, cols = vas_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()
