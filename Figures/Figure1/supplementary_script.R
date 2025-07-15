library(Seurat)


###--------Supplementary figure 1a--------###















###--------Supplementary figure 1b--------###

DCIS_copykat_cnv_bcmdcis01 <- read.table("/volumes/USR2/rumwei/Analysis/DCIS/CopyKat_res/CopyKat_DCIS_SC_NegCon/BCMDCIS01T_copykat_CNA_results.txt", header = T)

DCIS_copykat_cnv_sel_all <- read.csv('/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/tumor/CNV_heatmap/DCIS_copykat_cnv_sel_all.csv', header = T)
sample_record <- read.csv('/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/tumor/CNV_heatmap/sample_record.csv', header = T, row.names = 1)

colnames(DCIS_copykat_cnv_sel_all) <- c('cell', paste(DCIS_copykat_cnv_bcmdcis01$chrom, DCIS_copykat_cnv_bcmdcis01$chrompos, sep = '_'))

DCIS_copykat_cnv_sel_all_sample <- cbind(sample_record, DCIS_copykat_cnv_sel_all)
colnames(DCIS_copykat_cnv_sel_all_sample)[1] <- 'sample'

#tissue_type_upd1 <- tissue_type_upd[, c(1:3)]
tissue_type_upd1 <- dcis_clinical_0924[, c(1,2,5)] #using updated classfication
tissue_type_upd1$sample <- sub('_pc1|_24hTis|_cryo|_24htis', '', tissue_type_upd1$sample)

DCIS_copykat_cnv_sel_all_sample_tissue <- merge(tissue_type_upd1, DCIS_copykat_cnv_sel_all_sample, by = 'sample')

colnames(DCIS_copykat_cnv_sel_all_sample_tissue)[c(5:12171)] <- paste('chr', colnames(DCIS_copykat_cnv_sel_all_sample_tissue)[c(5:12171)], sep = '_')

#only select ER+ samples
#DCIS_copykat_cnv_sel_all_sample_tissue <- DCIS_copykat_cnv_sel_all_sample_tissue %>% filter(!ER == 'neg') %>% filter(!tissue_upd_2024_final == 'mucinous')
DCIS_copykat_cnv_sel_all_sample_tissue <- DCIS_copykat_cnv_sel_all_sample_tissue %>% filter(!ER == 'neg') %>% filter(!tissue_upd_0924 == 'mucinous')

#binarize the matrix and caculate the freq (thuis is for all and per tissue; cutoff for amp 0.2; loss -0.2)
##this is for both gain and loss; so run twice but need to change parameters!!!!!!!
{
  DCIS_copykat_cnv_sel_all_sample_tissue_bin <- DCIS_copykat_cnv_sel_all_sample_tissue %>%
    pivot_longer(cols = starts_with("chr"), names_to = "chr", values_to = "value") %>%
    #mutate(value = ifelse(value < -0.2, 1, 0)) %>%  ##loss
    mutate(value = ifelse(value > 0.2, 1, 0)) %>% #gain
    ##group_by(sample, tissue_upd_2024_final, chr) %>%
    group_by(sample, tissue_upd_0924, chr) %>%
    summarize(count = sum(value > 0)) %>%
    ungroup()
  
  #only consider sample with cell # = 100
  DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel <- DCIS_copykat_cnv_sel_all_sample_tissue_bin %>% dplyr::filter(sample %in% names(table(DCIS_copykat_cnv_sel_all_sample_tissue$sample)[table(DCIS_copykat_cnv_sel_all_sample_tissue$sample) == 100]))
  #only consider events more than 5
  DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel$count <- ifelse(DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel$count < 5, 0, DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel$count)
  
  #frequency
  DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel_freq <- DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel %>%
    #group_by(tissue_upd_2024_final, chr) %>%
    group_by(tissue_upd_0924, chr) %>%
    summarize(frequency = mean(count > 1))
  
  ######only for loss events!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel_freq$frequency <- -1 * DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel_freq$frequency
  
  ##choose tissues
  #tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq <- DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel_freq %>% filter(tissue_upd_2024_final == 'DCIS_yes')
  tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq <- DCIS_copykat_cnv_sel_all_sample_tissue_bin_sel_freq
  tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq <- tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq %>% mutate(split_strings = strsplit(chr, "_"),
                                                                                                                            chromosome = sapply(split_strings, function(x) x[[2]]))
  tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord <- tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq %>% 
    mutate(chr = factor(chr, levels = paste('chr', paste(DCIS_copykat_cnv_bcmdcis01$chrom, DCIS_copykat_cnv_bcmdcis01$chrompos, sep = '_'), sep = '_'))) %>%
    arrange(chr)
  tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord$chromosome <- factor(tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord$chromosome, levels = unique(DCIS_copykat_cnv_bcmdcis01$chrom))
  # ggplot(tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord, aes(x = chr, y = frequency, group = tissue_upd_2024_final)) + geom_line(color = 'blue') + geom_area(aes(y=frequency), fill='blue', alpha=1, position = "identity") + 
  #   facet_grid(.~ chromosome, scales = "free_x", space = "free_x") + theme_bw() + theme(strip.text = element_text(size=10), panel.spacing = unit(-.5, 'pt'), plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  #   scale_y_continuous(expand = c(0, 0)) + rremove('x.text') + rremove('x.ticks') + rremove('xlab') + rremove('grid')
  
}

####frequency plot for all
#record gain
all_copykat_cnv_sel_all_sample_all_bin_sel_freq_amp <- tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord
colnames(all_copykat_cnv_sel_all_sample_all_bin_sel_freq_amp) <- c('all', 'pos', 'amp', 'str', 'chr')
#record loss
all_copykat_cnv_sel_all_sample_all_bin_sel_freq_loss <- tissue_copykat_cnv_sel_all_sample_tissue_bin_sel_freq_ord
colnames(all_copykat_cnv_sel_all_sample_all_bin_sel_freq_loss) <- c('all', 'pos', 'loss', 'str', 'chr')

all_copykat_gain_loss <- merge(all_copykat_cnv_sel_all_sample_all_bin_sel_freq_amp, all_copykat_cnv_sel_all_sample_all_bin_sel_freq_loss, by = c('pos', 'all'))
all_copykat_gain_loss_pos <- all_copykat_gain_loss[,c(1,2,3,5,6)]
colnames(all_copykat_gain_loss_pos) <- c('pos', 'tissue', 'gain', 'chr', 'loss')

ggplot(all_copykat_gain_loss_pos) + geom_line(aes(pos, gain, group = tissue, color = tissue),linewidth = 0.7) + geom_line(aes(pos, loss, group = tissue, color = tissue),linewidth = 0.7) + 
  facet_grid(.~ chr, scales = "free_x", space = "free_x") + theme_bw() + theme(strip.text = element_text(size=10), panel.spacing = unit(-.5, 'pt'), plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  rremove('x.text') + rremove('x.ticks') + rremove('xlab') + rremove('grid') + scale_color_manual(values = c("#F9B90AFF", "#F24C3DFF"))

all_copykat_gain_loss_pos1 <- all_copykat_gain_loss_pos
all_copykat_gain_loss_pos1 <- all_copykat_gain_loss_pos1 %>% separate(pos, into = c("chr", "chr_num", "chr_pos"), sep = "_")
all_copykat_gain_loss_pos1$chr_num <- factor(all_copykat_gain_loss_pos1$chr_num, unique(DCIS_copykat_cnv_bcmdcis01$chrom))
all_copykat_gain_loss_pos1$chr_pos <- as.numeric(all_copykat_gain_loss_pos1$chr_pos)

ggplot(all_copykat_gain_loss_pos1) + geom_line(aes(chr_pos, gain, group = tissue, color = tissue),linewidth = 0.7) + geom_line(aes(chr_pos, loss, group = tissue, color = tissue),linewidth = 0.7) + 
  facet_grid(.~ chr_num, scales = "free_x", space = "free_x") + theme_bw() + theme(strip.text = element_text(size=10), panel.spacing = unit(-.5, 'pt'), plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  rremove('x.text') + rremove('x.ticks') + rremove('xlab') + rremove('grid') + scale_color_manual(values = c("#F9B90AFF", "#F24C3DFF", "blue"))


## Add breast cancer genes to freq plot ##
bc_gene_pos <- read.table('/volumes/USR1/siyuan/reference/Nick_BC_genes.sorted.bed', sep='\t') %>% 
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

#bc_gene_select_df$po <- paste(bc_gene_select_df$chr_s, bc_gene_select_df$mid, sep = '_')

# ggplot(all_copykat_gain_loss_pos1) + geom_line(aes(chr_pos, gain, group = tissue, color = tissue),linewidth = 0.85) + geom_line(aes(chr_pos, loss, group = tissue, color = tissue),linewidth = 0.85) + 
#   #ggrepel::geom_text_repel(data=bc_gene_select_df, aes(x=mid, y=0.8, label = bc_gene), size = 3) +
#   geom_vline(data=bc_gene_select_df, aes(xintercept=mid),  linetype='dotted', size=.5) +
#   facet_grid(.~ chr_num, scales = "free_x", space = "free_x") + theme_bw() + theme(strip.text = element_text(size=10), panel.spacing = unit(-.5, 'pt'), plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#   rremove('x.text') + rremove('x.ticks') + rremove('xlab') + rremove('grid') + scale_color_manual(values = c("#F9B90AFF", "#F24C3DFF"))

all_copykat_gain_loss_pos1$tissue <- factor(all_copykat_gain_loss_pos1$tissue, levels = c("DCIS_yes",  "synch_yes", "IDC_yes"))

ggplot(all_copykat_gain_loss_pos1) + geom_line(aes(chr_pos, gain, group = tissue, color = tissue),linewidth = 0.75) + geom_line(aes(chr_pos, loss, group = tissue, color = tissue), linewidth = 0.75) + 
  ggrepel::geom_text_repel(data=bc_gene_select_df, aes(x=mid, y=0.9, label = bc_gene), size = 3) +
  geom_vline(data=bc_gene_select_df, aes(xintercept=mid),  linetype='dotted', size=.5) +
  facet_grid(.~ chr_num, scales = "free_x", space = "free_x") + theme_bw() + theme(strip.text = element_text(size=10), panel.spacing = unit(-.5, 'pt'), plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  rremove('x.text') + rremove('x.ticks') + rremove('xlab') + rremove('grid') + scale_color_manual(values = c("#F9B90AFF", "#41B6E6FF", "#F24C3DFF"))+ ylim(-1, 1)