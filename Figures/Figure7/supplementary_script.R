
####################spatial analysis####################
dcis_300_xenium_srt_12cons_all_rm <- readRDS('./hbca_300_xenium_srt_6cons_downsample_092824.rds')

xenium_dcis_srt_merge_list <- SplitObject(dcis_300_xenium_srt_12cons_all_rm, split.by = "sample")

dcis_300_xenium_srt_12cons_merge_sp_df <- data.frame(Embeddings(xenium_dcis_srt_merge_list[[1]], reduction = 'SP'), cell_anno=xenium_dcis_srt_merge_list[[1]]$cellstate_ml_merge, sample='BCMHBCA59R2')
for (i in 2:length(xenium_dcis_srt_merge_list)) {
  print(i)
  srt_i <- xenium_dcis_srt_merge_list[[i]]
  name_i <- names(xenium_dcis_srt_merge_list)[i]
  temp_df <- data.frame(Embeddings(srt_i, reduction = 'SP'), cell_anno=srt_i$cellstate_ml_merge, sample=name_i)
  temp_df$SP_1 <- temp_df$SP_1 + max(dcis_300_xenium_srt_12cons_merge_sp_df$SP_1) + 2000
  dcis_300_xenium_srt_12cons_merge_sp_df <- rbind(dcis_300_xenium_srt_12cons_merge_sp_df, temp_df)
}

message('CellTrek colocalizing ...')

dcis_300_xenium_srt_12cons_merge_sp_df_dt <- CellTrek:::DT_boot_mst(dcis_300_xenium_srt_12cons_merge_sp_df[, c(3, 4)], 
                                                                    coord_df=dcis_300_xenium_srt_12cons_merge_sp_df[, c(1, 2)], col_cell='cell_anno', boot_n=5, dist_cutoff=1000, prop = 0.6)


hbca_300_xenium_srt_6cons_merge_sp_df_dt <- readRDS('./hbca_300_xenium_srt_merge_sp_df_dt_0928.rds')
hbca_300_xenium_srt_6cons_merge_sp_df_dt_boot_dis <- apply(hbca_300_xenium_srt_6cons_merge_sp_df_dt$boot_array, c(1, 2), mean) #, na.rm = T

hbca_300_xenium_srt_6cons_merge_sp_df_dt_boot_dis_median <- apply(hbca_300_xenium_srt_6cons_merge_sp_df_dt$boot_array, c(1, 2), median)

#normalization
hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_dis <- (hbca_300_xenium_srt_6cons_merge_sp_df_dt_boot_dis - min(hbca_300_xenium_srt_6cons_merge_sp_df_dt_boot_dis)) / (max(hbca_300_xenium_srt_6cons_merge_sp_df_dt_boot_dis) - min(hbca_300_xenium_srt_6cons_merge_sp_df_dt_boot_dis))
hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox <- 1 - hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_dis

hbca_dcis_xenium_clu_forhbca2 <- hbca_dcis_xenium_clu_forhbca %>% filter(id %in% colnames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox))

hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox
#hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm[hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm < 0.76] <- 0

hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm[hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm < 0.8] <- 0

# hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm[, colnames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm) != "fibro_cycling"]
# hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm[rownames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm) != "fibro_cycling", ]
# hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm[, colnames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm) != "cab"]
# hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm[rownames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm) != "cab", ]

scoloc_vis_upd(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm, hbca_dcis_xenium_clu_forhbca2)
# scoloc_vis_upd(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm, hbca_xenium_clu, col_pal = c("1" = "#49A85FFF", "2" = "#894FC6FF", "3" = "#1A91EBFF", "4" = "#894FC6FF",
#                                                                                                         "5" = "#81E74AFF", "6" = "#FFACAAFF", "7" = "#894FC6FF"))

hbca_xenium_clu <- data.frame(id = names(membership(hbca_xenium_community_obj)), clu = as.character(membership(hbca_xenium_community_obj)), 
                              freq=c(table(hbca_300_xenium_srt_6cons_all_merge_rm$cellstate_ml_merge)[names(membership(hbca_xenium_community_obj))]))

# hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm1 <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm[, colnames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm) != "fibro_cycling"]
# hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm1 <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm1[rownames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm1) != "fibro_cycling", ]
# 
# scoloc_vis_upd(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm1, hbca_xenium_clu)

hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rm
hbca_old_row_names <- c("adipocyte", "artery", "basal_intef", "Basal_others", "capillary", "CD4_Intef", "fibro_cfd", 
                        "fibro_igf1", "fibro_sfrp4", "IgA", "IgG", "lumsec", "lymphatic", "Macro_APOC1", "Macro_ISG15", "Normal_like_lumhr", "vein")
hbca_new_row_names <- c("Adipocyte", "Artery", "Basal-IFN", "Basal-major", "Capillary", "CD4-IFN",  "Fibro-prematrix",
                        "Fibro-matrix", "Fibro_sfrp4", "Plasma-IgA", "Plasma-IgG", "Lumsec", "Lymphatic", "Macro-lipo", "Macro-IFN", "Lumhr", "Vein")

hbca_old_col_names <- c("adipocyte", "artery", "basal_intef", "Basal_others", "capillary", "CD4_Intef", "fibro_cfd", 
                        "fibro_igf1", "fibro_sfrp4", "IgA", "IgG", "lumsec", "lymphatic", "Macro_APOC1", "Macro_ISG15", "Normal_like_lumhr", "vein")
hbca_new_col_names <- c("Adipocyte", "Artery", "Basal-IFN", "Basal-major", "Capillary", "CD4-IFN",  "Fibro-prematrix",
                        "Fibro-matrix", "Fibro_sfrp4", "Plasma-IgA", "Plasma-IgG", "Lumsec", "Lymphatic", "Macro-lipo", "Macro-IFN", "Lumhr", "Vein" )

# Replace row names
rownames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename)[rownames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename) %in% hbca_old_row_names] <- hbca_new_row_names

# Replace column names
colnames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename)[colnames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename) %in% hbca_old_col_names] <- hbca_new_col_names

hbca_xenium_clu_25 <- hbca_xenium_clu
hbca_xenium_clu_25$id_ori <- hbca_xenium_clu_25$id
hbca_xenium_clu_25$id[hbca_xenium_clu_25$id %in% hbca_old_col_names] <- hbca_new_col_names
rownames(hbca_xenium_clu_25) <- hbca_xenium_clu_25$id


# hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename[, colnames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename) != "Basal-IFN"]
# hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename[rownames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename) != "Basal-IFN", ]
# hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename[, colnames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename) != "CD4-IFN"]
# hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename <- hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename[rownames(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename) != "CD4-IFN", ]

scoloc_vis_upd(hbca_300_xenium_srt_6cons_sp_boot_logORdist_scaled_prox_rename, hbca_xenium_clu_25, col_pal = c("1" = "#49A85FFF", "2" = "#FFB900FF", "3" = "#1A91EBFF", "4" = "#9B7362FF",
                                                                                                               "5" = "#D34817FF", "6" = "#FFACAAFF", "7" = "#C78EF0FF"))

