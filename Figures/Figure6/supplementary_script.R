library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(future)


# label transfer proportion analysis--------------------------------------------------------------
hbca_epi_15k <- readRDS("./hbca_epi_upd_srtint.rds")

Idents(hbca_epi_15k) <- hbca_epi_15k$celltype

all_tumor_obj <- readRDS("./tumor_obj_rm_sub_diet.rds")
dcis_tumor <- all_tumor_obj
DefaultAssay(dcis_tumor) <- "RNA"

##without removing proliferating lumsec cells to predict
cell_num <- ncol(dcis_tumor)
k.filter <- min(200, cell_num)

##remove proliferating lumsec cells to predict
hbca_epi_15k_sub <- subset(hbca_epi_15k, celltype != "Lumsec_Prolif")

###
DefaultAssay(hbca_epi_15k_sub) <- "integrated"
dcis_tumor_sc.anchors1 <- FindTransferAnchors(reference = hbca_epi_15k_sub, query = dcis_tumor, dims = 1:30, k.filter = k.filter) ##this function doesn't work well on seurat 3.2
dcis_tumor_sc.prediction1_1 <- TransferData(anchorset = dcis_tumor_sc.anchors1, refdata = hbca_epi_15k_sub$celltype, dims = 1:30)

write.table(dcis_tumor_sc.prediction1_1, "dcis_tumor_sc.prediction_rm_proliflumsec_upd_integ_1112.txt", quote = F)

dcis_tumor_label_trans_inte <- dcis_tumor_sc.prediction1_1

all_tumor_obj <- readRDS("./tumor/tumor_obj_rm_sub_diet.rds")
all_tumor_obj@meta.data[, grepl('^sheet_', colnames(all_tumor_obj@meta.data))] <- NULL
all_tumor_obj_meta <- data.frame(sample_id = all_tumor_obj$orig.ident)
dcis_tumor_labtrans_sample <- merge(dcis_tumor_label_trans_inte[, c(1,5)], all_tumor_obj_meta, by = 0)

dcis_tumor_labtrans_sample_tis_er <- dcis_tumor_labtrans_sample_tis %>% filter(!ER == 'neg') %>% filter(!tissue_upd_2024_final %in% c('DCIS_no', 'mucinous'))
dcis_tumor_lt_inte_freq_24 <- data.frame(table(dcis_tumor_labtrans_sample_tis_er$sample_id, dcis_tumor_labtrans_sample_tis_er$predicted.id))

dcis_tumor_lt_inte_freq_reord_24 <- dcis_tumor_lt_inte_freq_24 %>% group_by(Var1) %>% mutate(lumsec_prop = sum(Freq[Var2=="Lumsec"], na.rm=T)/sum(Freq, na.rm=T))
#dcis_tumor_lt_inte_freq_reord_24_sample <- merge(dcis_tumor_lt_inte_freq_reord_24, tissue_type_upd_0725, by.x = 'Var1', by.y = 'sample')
dcis_tumor_lt_inte_freq_reord_24_sample <- merge(dcis_tumor_lt_inte_freq_reord_24, dcis_clinical_0924[, c(1,2,3,5)], by.x = 'Var1', by.y = 'sample')
dcis_tumor_lt_inte_freq_reord_24_sample$tissue_upd_0924 <- factor(dcis_tumor_lt_inte_freq_reord_24_sample$tissue_upd_0924, levels = c("DCIS_yes",  "synch_yes",  "IDC_yes"))
dcis_tumor_lt_inte_freq_reord_24_sample_ord <- dcis_tumor_lt_inte_freq_reord_24_sample %>% arrange(tissue_upd_0924, lumsec_prop)
ggplot() +
  geom_col(data = dcis_tumor_lt_inte_freq_reord_24_sample_ord, aes(x = Var1, y = Freq, fill = Var2), position = "fill") + 
  theme_bw() +
  scale_fill_manual(values = c('#69B3E7FF', '#31D64DFF')) + theme(axis.text.x=element_text(angle=50,hjust=1,vjust=1)) + 
  scale_x_discrete(limits = unique(dcis_tumor_lt_inte_freq_reord_24_sample_ord$Var1))

tissue_type_upd_0924 <- dcis_clinical_0924[, c(1,2,3,5)]
tissue_type_upd_0924$PR <- NULL
tissue_type_upd_0924_merge <- merge(tissue_type_upd_0725, tissue_type_upd_0924, by = 'sample')

#plot figure legend
DCIS_hr_type_reord_24 <- tissue_type_upd_0924_merge[match(unique(dcis_tumor_lt_inte_freq_reord_24_sample_ord$Var1), tissue_type_upd_0924_merge$sample),]
DCIS_hr_type_reord1_24 <- DCIS_hr_type_reord_24 %>% dplyr::select(-c(ER.x, tissue_upd_2024_final, ER.y)) %>% pivot_longer(!sample, names_to = "status", values_to = "value")
ggplot(DCIS_hr_type_reord1_24, aes(x = sample, y = status, fill = value)) + 
  geom_tile(size = 0.3) +xlab("") + ylab("") +
  scale_fill_manual(values = c(tissue_type_col[2], tissue_type_col[4], '#ECC6A2FF', '#4B3A51FF', tissue_type_col[3], '#333333FF', tissue_type_col[2],  'grey')) + 
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")))+ 
  scale_x_discrete(limits = unique(dcis_tumor_lt_inte_freq_reord_24_sample_ord$Var1))


# metaprograms frequency plot-------------------------------------------------------------
meta_module_freq_mat <- matrix(NA, 14, 85)
n_row <- 1
for(i in unique(tumor_nmf_intersect_rmrg_group$grp)){
  meta_i <- rownames(tumor_nmf_intersect_rmrg_group)[tumor_nmf_intersect_rmrg_group$grp==i]
  tumor_nmf_programs_top50_dat <- tumor_nmf_programs_top50[, meta_i]
  meta_module_freq_mat[n_row, ] <- unlist(lapply(tumor_nmf_programs_top50_ori, check_existence_and_return, matrix1 = tumor_nmf_programs_top50_dat))
  n_row <- n_row + 1
}
meta_module_freq_mat_sum <- data.frame(m = names(module_gene_top50), freq = apply(meta_module_freq_mat, 1, sum))

ggplot(meta_module_freq_mat_sum, aes(x=m, y=freq)) + 
  geom_bar(stat = "identity", fill = "#FFB900FF") + scale_x_discrete(limits = meta_module_freq_mat_sum$m) +
  theme_classic() +theme(axis.text = element_text(size = 12)) + theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1,size = 10)) #coord_flip


# correlation between metaprograms identified in this work and pan-cancer metaporgrams-------------------------------------------------------------
pancancer_nmf <- read.csv("./pancancer_nmf.csv")
pancancer_jac_mat <- matrix(NA, 14, 41)

n_row <- 1
for (i in 1:14){
  for(j in 1:41){
    pancancer_jac_mat[n_row, j] <- jaccard_index(module_gene_top50[[i]], pancancer_nmf[, j])
  }
  n_row <- n_row + 1
}
colnames(pancancer_nmf) <- gsub('\\.+', '_', colnames(pancancer_nmf))
colnames(pancancer_jac_mat) <- colnames(pancancer_nmf)
#rownames(pancancer_jac_mat) <- paste("M", c(1:13), sep = "")
rownames(pancancer_jac_mat) <- names(module_gene_top50)
pancancer_jac_mat_w <- reshape2::melt(pancancer_jac_mat, varnames = c("Row", "Column"))
ggplot(data = pancancer_jac_mat_w, aes(x=Column, y=Row, fill=value)) + 
  geom_tile() + scale_fill_gradientn(colors = (paletteer::paletteer_c("viridis::rocket", n = 100)), name="Jaccard Index") + 
  theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1,size = 10),axis.text.y = element_text(size=10)) + ylab("Module") + xlab("Pan-cancer Module")

pheatmap(pancancer_jac_mat, color = (paletteer::paletteer_c("viridis::rocket", n = 100)), clustering_method = "ward.D2")


# correlation between metaprograms identified in this work-------------------------------------------------------------
er_tumor_obj_diet <- AddModuleScore(er_tumor_obj_diet, features = module_gene_top50, ctrl = 50)

colnames(er_tumor_obj_diet@meta.data)[grepl('^Cluster', colnames(er_tumor_obj_diet@meta.data))] <- names(module_gene_top50)
dcis_sc_all_module_score <- er_tumor_obj_diet@meta.data[, names(module_gene_top50)]
dcis_sc_all_module_score_cor <- cor(dcis_sc_all_module_score, method = "pearson")
#dcis_sc_all_module_score_cor <- cor(dcis_sc_all_module_score, method = "spearman")
diag(dcis_sc_all_module_score_cor) <- NA

pheatmap::pheatmap(
  dcis_sc_all_module_score_cor,
  color = viridis::viridis(50),
  angle_col = 45,
  clustering_method = 'ward.D2'
)


