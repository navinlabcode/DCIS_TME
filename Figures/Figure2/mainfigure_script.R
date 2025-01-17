
fibro_all_diet <- readRDS("")

# fig2a fibro UMAP-------------------------------------------------------------------
DimPlot(fibro_all_diet, raster = F, cols = fibro_col, label = F)


# fig2c fibro markers-------------------------------------------------------------------
VlnPlot(fibro_all_diet, c('COL15A1', 'IGF1', 'DLK1', 'IGFBP2', 'ACTA2', 'FAP', 'TAGLN', 'POSTN'), pt.size = 0, cols = fibro_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()

fibro_mmp_marker <- fibro_all_marker %>% filter(grepl("^MMP", gene)) %>% pull(gene)
VlnPlot(fibro_all_diet, fibro_mmp_marker, pt.size = 0, cols = fibro_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()
DotPlot(fibro_all_diet, features = unique(fibro_mmp_marker)) + RotatedAxis() + scale_color_gradientn(colours = viridis::viridis(20)) & coord_flip()

 

vascular_all_diet <- readRDS("")

# fig2i vasuclar markers -----------------------------------------------------------
VlnPlot(vascular_all_diet, c('APLNR', 'INSR', 'FLT1', 'KDR', 'NRP1', 'VWA1', 'COL4A1', 'COL4A2', 'COL6A2','LAMA4','LAMB1'), pt.size = 0, cols = vas_col, stack = T, flip = T, fill.by = 'ident') & NoLegend()
