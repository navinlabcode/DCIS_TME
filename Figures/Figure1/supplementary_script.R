library(Seurat)
library(dplyr)
library(infercna)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(cowplot)


###--------Supplementary figure 1a--------###
pdf('./all_tumor_cnv_heatmap.pdf', width=12, height=18.5)

sample_cell_df <- read.csv('./all_tumor_sample_cell.csv', header = T)
DCIS_copykat_cnv_sel <- read.csv('./DCIS_copykat_cnv_sel_all.csv', header = T, row.names = 1)
sample_record <- read.csv('./sample_record.csv', header = T, row.names = 1)
DCIS_copykat_cnv <- read.table("./BCMDCIS01T_copykat_CNA_results.txt", header = T)

DCIS_copykat_cnv_sel_all <- cbind(sample_record, DCIS_copykat_cnv_sel)
rownames(DCIS_copykat_cnv_sel_all) <- rownames(DCIS_copykat_cnv_sel)
colnames(DCIS_copykat_cnv_sel_all)[1] <- 'sample'

##clinical info
dcis_clinical_0924 <- read.csv('./dcis_project_tissuetype_0916.csv', header = T)
tissue_type_upd1 <- dcis_clinical_0924[, c(1,2,5)] 
tissue_type_upd1$sample <- sub('_pc1|_24hTis|_cryo|_24htis', '', tissue_type_upd1$sample)

DCIS_copykat_cnv_sel_all_sample_tissue <- right_join(tissue_type_upd1, DCIS_copykat_cnv_sel_all, by = 'sample')
rownames(DCIS_copykat_cnv_sel_all_sample_tissue) <- rownames(DCIS_copykat_cnv_sel_all)
DCIS_copykat_cnv_sel_all_sample_tissue <- DCIS_copykat_cnv_sel_all_sample_tissue %>% filter(!ER == 'neg') %>% filter(!tissue_upd_0924 == 'mucinous')
DCIS_copykat_cnv_sel_all_sample_tissue <- arrange(DCIS_copykat_cnv_sel_all_sample_tissue, match(DCIS_copykat_cnv_sel_all_sample_tissue$tissue_upd_0924, c("DCIS_yes", "synch_yes", "IDC_yes")))

sample_record <- DCIS_copykat_cnv_sel_all_sample_tissue$sample
tissue_record <- DCIS_copykat_cnv_sel_all_sample_tissue$tissue_upd_0924

DCIS_copykat_cnv_sel_all_sample_tissue$sample <- NULL
DCIS_copykat_cnv_sel_all_sample_tissue$ER <- NULL
DCIS_copykat_cnv_sel_all_sample_tissue$tissue_upd_0924 <- NULL
DCIS_copykat_cnv_sel_all_scale1 <- DCIS_copykat_cnv_sel_all_sample_tissue

corresponding_values <- paste0("p", seq_along(unique(sample_record)))
sample_record1 <- data.frame(x = corresponding_values[match(sample_record, unique(sample_record))])

chrom <- DCIS_copykat_cnv$chrom
hb <- HeatmapAnnotation(chr=factor(chrom,levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)),show_legend = F,annotation_name_side="left",which="column",
                        col=list(chr = c("1"="grey","2"="black","3"="grey","4"="black","5"="grey","6"="black","7"="grey","8"="black","9"="grey","10"="black","11"="grey",
                                         "12"="black","13"="grey","14"="black","15"="grey","16"="black","17"="grey","18"="black","19"="grey","20"="black","21"="grey","22"="black","23"="grey")))

ht <- ComplexHeatmap::Heatmap(DCIS_copykat_cnv_sel_all_scale1, cluster_rows = F ,cluster_columns=F, border="black", raster_device = "png", show_row_names = F, show_column_names = F, row_gap = unit(0,"mm"), 
                              top_annotation = hb, column_gap = unit(0,"mm"), column_split = DCIS_copykat_cnv$chrom, row_order = rownames(DCIS_copykat_cnv_sel_all_scale1), row_split = factor(sample_record1$x, 
                              levels = paste0('p', c(1:length(unique(sample_record1$x))))), col = circlize::colorRamp2(c(min(DCIS_copykat_cnv_sel_all_scale1), -0.4, 0, 0.55, max(DCIS_copykat_cnv_sel_all_scale1)), 
                              c("#26456EFF", "#1C5E9DFF", "white", "#B3101BFF", "#9C0824FF")), cluster_row_slices = F)

row_ht <- Heatmap(tissue_record, name = "type", width = 1)

p_com <- cowplot::plot_grid(grid::grid.grabExpr(ComplexHeatmap::draw(ht + row_ht)))

print(p_com)


###--------Supplementary figure 1b--------###
hbca_dcis_ibc.integrated.b1.b2.clean.cellanno <- readRDS('./hbca_dcis_ibc.integrated.b1.b2.clean.cellanno.rds')
hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$celltype_com <- hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$celltype_upd1
hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$celltype_com <- factor(hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$celltype_com, 
                                                                     levels = c("Basal", "lumhr", "lumsec",  "Bcell", "NK_T", "Myeloid", "Fibroblast", "Vascular", "Lymphatic", "PeriVas"))

hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$tumor_normal <- ifelse(rownames(hbca_dcis_ibc.integrated.b1.b2.clean.cellanno@meta.data) %in% rownames(tumor_obj_rm_meta_data), 
                                                                     'aneuploid', 'diploid')

hbca_dcis_ibc.integrated.all.meta <- data.frame(sample = hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$orig.ident, proj_ori = hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$project,
                                                celltype = hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$celltype_com, ploidy = hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$tumor_normal)


hbca_dcis_ibc.integrated.all.meta.tumor.freq <- data.frame(table(hbca_dcis_ibc.integrated.all.meta$sample, hbca_dcis_ibc.integrated.all.meta$ploidy))
hbca_dcis_ibc.integrated.all.meta.tumor.tis <- merge(hbca_dcis_ibc.integrated.all.meta.tumor.freq, dcis_clinical_0924, by.x = 'Var1', by.y = 'sample', all.x = T)
hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp <- hbca_dcis_ibc.integrated.all.meta.tumor.tis[, c(1,2,3,4,7)]

hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp <- hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp %>% filter(!ER == 'neg') %>% filter(!tissue_upd_0924 == 'mucinous')

hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp$Proportion <- hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp$Freq / ave(hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp$Freq, hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp$Var1, FUN = sum)
hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp$tissue <- factor(hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp$tissue_upd_0924, levels = c("Normal", "DCIS_no", "DCIS_yes", "synch_no", "synch_yes", "IDC_no", "IDC_yes"))

hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp <- hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp %>%
  arrange(tissue, desc(Var2), desc(Proportion))

ggplot(hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp, aes(fill=Var2, y=Proportion, x=Var1)) + 
  geom_bar(stat="identity") + scale_x_discrete(limits = unique(hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp$Var1)) + theme_bw()+
  theme(axis.text.x=element_text(angle=50,hjust=1,vjust=1),axis.text.y = element_text(size=12)) + scale_fill_manual(values = c('#DC863BFF', '#a4ceb7'))


###--------Supplementary figure 1c--------###
DCIS_copykat_cnv_bcmdcis01 <- read.table("./BCMDCIS01T_copykat_CNA_results.txt", header = T)

DCIS_copykat_cnv_sel_all <- read.csv('./DCIS_copykat_cnv_sel_all.csv', header = T)
sample_record <- read.csv('./sample_record.csv', header = T, row.names = 1)

colnames(DCIS_copykat_cnv_sel_all) <- c('cell', paste(DCIS_copykat_cnv_bcmdcis01$chrom, DCIS_copykat_cnv_bcmdcis01$chrompos, sep = '_'))

DCIS_copykat_cnv_sel_all_sample <- cbind(sample_record, DCIS_copykat_cnv_sel_all)
colnames(DCIS_copykat_cnv_sel_all_sample)[1] <- 'sample'

tissue_type_upd1 <- dcis_clinical_0924[, c(1,2,5)]
tissue_type_upd1$sample <- sub('_pc1|_24hTis|_cryo|_24htis', '', tissue_type_upd1$sample)

DCIS_copykat_cnv_sel_all_sample_tissue <- merge(tissue_type_upd1, DCIS_copykat_cnv_sel_all_sample, by = 'sample')

colnames(DCIS_copykat_cnv_sel_all_sample_tissue)[c(5:12171)] <- paste('chr', colnames(DCIS_copykat_cnv_sel_all_sample_tissue)[c(5:12171)], sep = '_')

DCIS_copykat_cnv_sel_all_sample_tissue <- DCIS_copykat_cnv_sel_all_sample_tissue %>% filter(!ER == 'neg') %>% filter(!tissue_upd_0924 == 'mucinous')

##this is for both gain and loss; so run twice but need to change parameters
{
  DCIS_copykat_cnv_sel_all_sample_tissue_bin <- DCIS_copykat_cnv_sel_all_sample_tissue %>%
    pivot_longer(cols = starts_with("chr"), names_to = "chr", values_to = "value") %>%
    #mutate(value = ifelse(value < -0.2, 1, 0)) %>%  ##loss
    mutate(value = ifelse(value > 0.2, 1, 0)) %>% #gain
    group_by(sample, tissue_upd_0924, chr) %>%
    summarize(count = sum(value > 0)) %>%
    ungroup()
  
  #only consider sample with cell # = 100
  DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel <- DCIS_copykat_cnv_sel_all_sample_tissue_bin %>% dplyr::filter(sample %in% names(table(DCIS_copykat_cnv_sel_all_sample_tissue$sample)[table(DCIS_copykat_cnv_sel_all_sample_tissue$sample) == 100]))
  #only consider events more than 5
  DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel$count <- ifelse(DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel$count < 5, 0, DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel$count)
  
  #frequency
  DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel_freq <- DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel %>%
    group_by(tissue_upd_0924, chr) %>%
    summarize(frequency = mean(count > 1))
  
  ######only for loss events
  DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel_freq$frequency <- -1 * DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel_freq$frequency
  
  tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq <- DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel_freq
  tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq <- tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq %>% mutate(split_strings = strsplit(chr, "_"),
                                                                                                                            chromosome = sapply(split_strings, function(x) x[[2]]))
  tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord <- tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq %>% 
    mutate(chr = factor(chr, levels = paste('chr', paste(DCIS_copykat_cnv_bcmdcis01$chrom, DCIS_copykat_cnv_bcmdcis01$chrompos, sep = '_'), sep = '_'))) %>%
    arrange(chr)
  tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord$chromosome <- factor(tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord$chromosome, levels = unique(DCIS_copykat_cnv_bcmdcis01$chrom))
}

#record gain
all_copykat_cnv_sel_all_sample_all_bin_sel_freq_amp <- tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord
colnames(all_copykat_cnv_sel_all_sample_all_bin_sel_freq_amp) <- c('all', 'pos', 'amp', 'str', 'chr')
#record loss
all_copykat_cnv_sel_all_sample_all_bin_sel_freq_loss <- tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord
colnames(all_copykat_cnv_sel_all_sample_all_bin_sel_freq_loss) <- c('all', 'pos', 'loss', 'str', 'chr')

all_copykat_gain_loss <- merge(all_copykat_cnv_sel_all_sample_all_bin_sel_freq_amp, all_copykat_cnv_sel_all_sample_all_bin_sel_freq_loss, by = c('pos', 'all'))
all_copykat_gain_loss_pos <- all_copykat_gain_loss[,c(1,2,3,5,6)]
colnames(all_copykat_gain_loss_pos) <- c('pos', 'tissue', 'gain', 'chr', 'loss')

all_copykat_gain_loss_pos1 <- all_copykat_gain_loss_pos
all_copykat_gain_loss_pos1 <- all_copykat_gain_loss_pos1 %>% separate(pos, into = c("chr", "chr_num", "chr_pos"), sep = "_")
all_copykat_gain_loss_pos1$chr_num <- factor(all_copykat_gain_loss_pos1$chr_num, unique(DCIS_copykat_cnv_bcmdcis01$chrom))
all_copykat_gain_loss_pos1$chr_pos <- as.numeric(all_copykat_gain_loss_pos1$chr_pos)

## Add breast cancer genes to freq plot ##
bc_gene_pos <- read.table('./BC_genes.sorted.bed', sep='\t') %>% 
  set_colnames(c('chr', 'start', 'end', 'gene')) %>% set_rownames(.$gene)

{
  bc_gene_select_df <- data.frame(bc_gene=c('NF2','AKT3','CCND1', 'ERBB2', 'MYC', 'PTEN', 'PPP2R2A', 'SHC1', 'TP53',  'RB1', 'GATA3', 'PLA2G10', 'WWOX', 'CCNE1', 'STK11'), 
                                  type=c('del','amp', 'amp', 'amp', 'amp', 'del', 'del', 'amp',  'del',  'del', 'amp',  'amp', 'amp', 'amp', 'amp'))
  bc_gene_select_df$chr <- bc_gene_pos[bc_gene_select_df$bc_gene, 'chr']
  bc_gene_select_df$start <- bc_gene_pos[bc_gene_select_df$bc_gene, 'start']
  bc_gene_select_df$end <- bc_gene_pos[bc_gene_select_df$bc_gene, 'end']
  bc_gene_select_df$mid <- (bc_gene_select_df$start + bc_gene_select_df$end)/2
  bc_gene_select_df$chr <- gsub("chr", "", bc_gene_select_df$chr)
  bc_gene_select_df$chr <- factor(bc_gene_select_df$chr, unique(DCIS_copykat_cnv_bcmdcis01$chrom))
  bc_gene_select_df$chr_num <- bc_gene_select_df$chr
}

all_copykat_gain_loss_pos1$tissue <- factor(all_copykat_gain_loss_pos1$tissue, levels = c("DCIS_yes",  "synch_yes", "IDC_yes"))

ggplot(all_copykat_gain_loss_pos1) + geom_line(aes(chr_pos, gain, group = tissue, color = tissue), linewidth = 0.75) + geom_line(aes(chr_pos, loss, group = tissue, color = tissue), linewidth = 0.75) + 
  ggrepel::geom_text_repel(data=bc_gene_select_df, aes(x = mid, y = 0.9, label = bc_gene), size = 3) +
  geom_vline(data = bc_gene_select_df, aes(xintercept = mid),  linetype = 'dotted', size = .5) +
  facet_grid(.~ chr_num, scales = "free_x", space = "free_x") + theme_bw() + theme(strip.text = element_text(size=10), panel.spacing = unit(-.5, 'pt'), plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  rremove('x.text') + rremove('x.ticks') + rremove('xlab') + rremove('grid') + scale_color_manual(values = c("#F9B90AFF", "#41B6E6FF", "#F24C3DFF"))+ ylim(-1, 1)


###--------Supplementary figure 1d--------###
hbca_dcis_ibc.integrated.all.meta <- data.frame(sample = hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$orig.ident, proj_ori = hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$project,
                                                celltype = hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$celltype_com, ploidy = hbca_dcis_ibc.integrated.b1.b2.clean.cellanno$tumor_normal)

hbca_dcis_ibc.integrated.all.celltype.freq <- data.frame(table(hbca_dcis_ibc.integrated.all.meta$sample, hbca_dcis_ibc.integrated.all.meta$celltype))

hbca_dcis_ibc.integrated.all.celltype.freq.sel <- hbca_dcis_ibc.integrated.all.celltype.freq %>% filter(Var1 %in% unique(hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp$Var1))

ggplot(hbca_dcis_ibc.integrated.all.celltype.freq.sel %>% filter(!Var2 %in% c('lumhr', 'lumsec')), aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position = 'fill', stat="identity") + scale_x_discrete(limits = unique(hbca_dcis_ibc.integrated.all.meta.tumor.tis.simp$Var1)) + theme_bw()+
  theme(axis.text.x=element_text(angle=50,hjust=1,vjust=1),axis.text.y = element_text(size=12)) + scale_fill_manual(values = c(celltype_col[1], celltype_col[4:10]))




