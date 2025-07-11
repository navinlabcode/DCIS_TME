library(dplyr)
library(ggpurb)

###--------Supplementary figure 4a--------###
basal_all_sum_count <- merge(basal_all_sum, dcis_clinical_0924[,c(1,2,5)], by.x = 'Var1', by.y = 'sample')
basal_all_sum_count_filter <- basal_all_sum_count %>% filter(tissue_upd_0924 %in% c("DCIS_yes", "synch_yes", "IDC_yes", "Normal"))
basal_all_sum_count_filter$tissue_upd_0924 <- factor(basal_all_sum_count_filter$tissue_upd_0924, levels = c("Normal", "DCIS_yes", "synch_yes", "IDC_yes"))
ggboxplot(basal_all_sum_count_filter, x = "Var2", y = "Freq", color = "tissue_upd_0924", add = "jitter",  ylab = "Count", palette = tissue_type_col)+
  stat_compare_means(aes(group = tissue_upd_0924), label = "p.format", method = "kruskal.test")+guides(color = guide_legend(nrow = 1))+
  rotate_x_text(angle = 45, hjust = 1, vjust = 1)


###--------Supplementary figure 4b--------###
cell_stat_plot <- 'basal-IFN'
supp_pwc <- cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot) %>% group_by(cells_State) %>% rstatix::dunn_test(Count ~ tissue_upd_0924, p.adjust.method = "BH") %>% rstatix::add_xy_position(x = "tissue_upd_0924")
ggplot(cells.pt.final_sum3_sel %>% filter(!cells_State %in% cell_stat_plot), aes(x = tissue_upd_0924, y = Count)) + geom_boxplot(outliers = F) + geom_jitter(aes(colour = tissue_upd_0924), size = 1) +
  facet_wrap(facets = vars(cells_State), nrow = 1) + theme_bw()+rremove("x.ticks") + rremove("x.text") + scale_color_manual(values = tissue_type_col) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + stat_pvalue_manual(supp_pwc, hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}")


###--------Supplementary figures 4c,d--------###



