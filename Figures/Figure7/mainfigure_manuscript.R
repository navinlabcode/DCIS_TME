# fig7a-cell state correlation-------------------------------------------------------------------
all_cell_count_24 <- read.csv('./filtered/all_cell_merged_anno.csv')
all_cell_count_24$unique_id <- paste(all_cell_count_24$sample, sub('-.*', '', all_cell_count_24$cell_id), sep = '_')

tumor_id <- data.frame(id = rownames(tumor_obj_rm_meta_data), sample = tumor_obj_rm_meta_data$orig.ident)
tumor_id$unique_id <- paste(tumor_id$sample, sub('-.*', '', tumor_id$id), sep = '_')

all_cell_count_24$celltype1 <- ifelse(all_cell_count_24$unique_id %in% tumor_id$unique_id, 'tumor', all_cell_count_24$celltype)

all_cell_count_24 <- read.csv('./filtered/all_cell_merged_anno_upd.csv', row.names = 1)

all_cellstate_count_24 <- all_cell_count_24[, c(2,1,6,4)]

#only extract tme cells and normal lumsec
all_tme_cell_count_24 <- all_cell_count_24 %>% filter(!celltype1 %in% c('lumhr', 'tumor', 'lumsec', 'Lymphatic'))
all_tme_cellstate_24_count1 <- all_tme_cell_count_24 %>% group_by(sample) %>% count(cellstate)
colnames(all_tme_cellstate_24_count1)[3] <- 'cellstate_n'
all_tme_celltype_24_count1 <- all_tme_cell_count_24 %>% group_by(sample) %>% count(celltype1)
colnames(all_tme_celltype_24_count1)[3] <- 'celltype_n'

all_cell_count_24_celltypestate <- all_tme_cell_count_24 %>% dplyr::select(sample, celltype1, cellstate) %>% distinct()

#keep celltype sum more than 20
all_tme_cellstate_24_count2_refor1_prop <- all_cell_count_24_celltypestate %>%
  left_join(all_tme_celltype_24_count1, by = c("sample", "celltype1")) %>%
  left_join(all_tme_cellstate_24_count1, by = c("sample", "cellstate")) %>%
  mutate(prop = ifelse(celltype_n > 20, cellstate_n / celltype_n, 0)) #mutate(prop = ifelse(sum_value > 20, count / sum_value, 0))

all_tme_cellstate_24_count2_refor1_prop_bac <- all_tme_cellstate_24_count2_refor1_prop
all_tme_cellstate_24_count2_refor1_prop1 <- all_tme_cellstate_24_count2_refor1_prop_bac %>% pivot_wider(names_from = cellstate, values_from = prop)
all_tme_cellstate_24_count2_cli_prop_sel1 <- merge(all_tme_cellstate_24_count2_refor1_prop1, tissue_type_upd[, c(1:3)], by = "sample")

all_tme_cellstate_24_count2_cli_prop_sel1[is.na(all_tme_cellstate_24_count2_cli_prop_sel1)] <- 0

#merge NMF and TME
tumor_obj_dat_type_trunc_sum_upd1 <- tumor_obj_dat_type_trunc_sum_upd
colnames(tumor_obj_dat_type_trunc_sum_upd1)[c(3:16)] <- paste('NMF', colnames(tumor_obj_dat_type_trunc_sum_upd1)[c(3:16)], sep = '_')
all_tme_nmf_prop_24 <- merge(tumor_obj_dat_type_trunc_sum_upd1, all_tme_cellstate_24_count2_cli_prop_sel1, by = "sample")

write.csv(all_tme_nmf_prop_24[, c(-1,-78, -79)], 'cellstate_nmf.csv', quote = F)

##generate the graph
#p values associated with cell state correlation
corr_fea_pval_24 <- matrix(NA, ncol = 75, nrow = 75)
colnames(corr_fea_pval_24) <- rownames(corr_fea_pval_24) <- names(all_tme_nmf_prop_24[, 3:77])
for (i in 1:ncol(all_tme_nmf_prop_24_cor)) {
  for (j in 1:ncol(all_tme_nmf_prop_24_cor)) {  # Skip diagonal elements (self-correlations)
    corr_fea_pval_24[i, j] <- cor.test(all_tme_nmf_prop_24[, 3:77][, i], all_tme_nmf_prop_24[, 3:77][, j], method = "spearman")$p.value
  }
}

all_tme_nmf_24_prop_cor_test <- all_tme_nmf_prop_24_cor
all_tme_nmf_24_prop_cor_test_network <- graph_from_adjacency_matrix(all_tme_nmf_24_prop_cor_test, weighted=T, mode="undirected", diag=F)
all_tme_nmf_24_prop_cor_test_network <- igraph::simplify(all_tme_nmf_24_prop_cor_test_network)
#all_tme_nmf_24_prop_cor_test_network <- igraph::delete.vertices(all_tme_nmf_24_prop_cor_test_network, which(igraph::degree(all_tme_nmf_24_prop_cor_test_network)==0))
edge_attr(all_tme_nmf_24_prop_cor_test_network, 'label') <- apply(as_edgelist(all_tme_nmf_24_prop_cor_test_network), 1, function(v) paste(v, collapse = '--'))

##remove edges with non-significant correlations
Gf_pval_24 <- graph_from_adjacency_matrix( 
  adjmatrix = corr_fea_pval_24, mode = 'undirected',
  weighted = TRUE, 
  diag = F)
edge_attr(Gf_pval_24, 'label') <- apply(as_edgelist(Gf_pval_24), 1, function(v) paste(v, collapse = '--'))
significant_edge_labels_24 <- E(Gf_pval_24)$label[E(Gf_pval_24)$weight < 0.01]; str(significant_edge_labels_24) #0.05
all_tme_nmf_24_prop_cor_test_network_rm1 <- delete_edges(all_tme_nmf_24_prop_cor_test_network, E(all_tme_nmf_24_prop_cor_test_network)[!label %in% significant_edge_labels_24])

##remove edges with negative correlation
all_tme_nmf_24_prop_cor_test_network_rm2 <- delete_edges(all_tme_nmf_24_prop_cor_test_network_rm1, E(all_tme_nmf_24_prop_cor_test_network_rm1)[weight <= 0])
all_tme_nmf_24_prop_cor_test_network_rm2 <- igraph::delete.vertices(all_tme_nmf_24_prop_cor_test_network_rm2, which(igraph::degree(all_tme_nmf_24_prop_cor_test_network_rm2)==0))

##determine # of clusters (res = 2)
num_clu_dat_res <- num_cluster(all_tme_nmf_24_prop_cor_test_network_rm2, res = seq(0.5, 3, by = 0.05), method = "mean")
num_clu_dat_res <- num_cluster(all_tme_nmf_24_prop_cor_test_network_rm2, res = seq(0.5, 3, by = 0.05), method = "median")

ggplot(num_clu_dat_res, aes(x = res)) +
  geom_line(aes(y = n_clu), size = 1, color = "blue") +
  geom_line(aes(y = n_mem), size = 1, color = "red") + theme_bw()

{
set.seed(13579)
##consensus clustering using louvain algorithm
cluster_memberships_24 <- list()
# Run Louvain clustering multiple times
num_runs <- 1000
num_nodes <- vcount(all_tme_nmf_24_prop_cor_test_network_rm2)
for (i in 1:num_runs) {
  cluster_obj <- igraph::cluster_louvain(all_tme_nmf_24_prop_cor_test_network_rm2, resolution = 2.1) #2.2; 2.1
  cluster_memberships_24[[i]] <- membership(cluster_obj)
}
# Create a consensus matrix
consensus_matrix_24 <- matrix(0, nrow = num_nodes, ncol = num_nodes)
for (i in 1:num_runs) {
  for (j in 1:num_nodes) {
    for (k in 1:num_nodes) {
      if (cluster_memberships_24[[i]][j] == cluster_memberships_24[[i]][k]) {
        consensus_matrix_24[j, k] <- consensus_matrix_24[j, k] + 1
      }
    }
  }
}
rownames(consensus_matrix_24) <- colnames(consensus_matrix_24) <- names(cluster_memberships_24[[1]])

# get final clustering results
graph_24 <- graph.adjacency(consensus_matrix_24, mode = "undirected", weighted = T, diag = F)
community_obj_24 <- igraph::cluster_louvain(graph_24, resolution = 2.1) #2.2
unique(igraph::membership(community_obj_24)); table(igraph::membership(community_obj_24))
}

all_tme_nmf_24_prop_cor_test_network1 <- all_tme_nmf_24_prop_cor_test_network_rm2
V(all_tme_nmf_24_prop_cor_test_network1)$clu <- as.character(membership(community_obj_24))

xyz_power <- 60 
l_wt_24 <- apply(get.edgelist(all_tme_nmf_24_prop_cor_test_network1), 1, 
              weight.community, 
              membership(community_obj_24), xyz_power, 1) #1

{
    l_24 <- layout_with_fr(all_tme_nmf_24_prop_cor_test_network1, weights=l_wt_24)
    ggraph(all_tme_nmf_24_prop_cor_test_network1, layout = l_24)+
    geom_edge_link0(aes(edge_width = weight), edge_colour = "#999BA0FF")+ #Draw edges as straight lines between nodes
    geom_node_point(aes(colour = clu), size = 5.5)+
    scale_colour_manual(values = cellstate_ecotype_col)+
    scale_edge_width(range = c(0.1, 1.5))+
    theme_graph()+
    theme(legend.position = "right")
}
ggraph(all_tme_nmf_24_prop_cor_test_network1, layout = l_24)+
  geom_edge_link0(aes(edge_width = weight),edge_colour = "#999BA0FF")+ #Draw edges as straight lines between nodes; alpha = weight
  geom_node_point(aes(fill = clu), size = 3, shape = 21)+
  scale_fill_manual(values = cellstate_ecotype_col)+
  scale_edge_width(range = c(0.1,1.5))+
  theme_graph()+
  theme(legend.position = "right")+
  geom_node_text(aes(label = name), repel = TRUE, point.padding = unit(0.2, "lines"), size = 3)


# fig7b-tissue distribution --------------------------------------------------------------
jaccard_similarity <- function(x, y) {
  intersection <- sum(x & y)  # Count of common 1s
  union <- sum(x | y)         # Count of 1s in either x or y
  return(intersection / union)
}
#construct a cluster list info
community_obj_24_info <- as.character(membership(community_obj_24))
names(community_obj_24_info) <- names(membership(community_obj_24))
community_obj_24_info_list <- list()
for (num in unique(community_obj_24_info)) {
  community_obj_24_info_list[[num]] <- names(community_obj_24_info[community_obj_24_info == num])
}

#construct a cluster list info (test)
community_obj_24_info_test <- as.character(membership(community_obj_24_test))
names(community_obj_24_info_test) <- names(membership(community_obj_24_test))
community_obj_24_info_list_test <- list()
for (num in unique(community_obj_24_info_test)) {
  community_obj_24_info_list_test[[num]] <- names(community_obj_24_info_test[community_obj_24_info_test == num])
}

names(cell_cooccur_col_test) <- c(1, 10, 11,2:9)
names(cell_cooccur_col) <- c(1, 10, 11, 2:9)

#for colors
cellstate_ecotype_col1 <- cellstate_ecotype_col
names(cellstate_ecotype_col1) <- c(1:8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_tme_nmf_prop_0924$tissue_upd_0924_upd <- all_tme_nmf_prop_0924$tissue_upd_0924
all_tme_nmf_prop_0924$tissue_upd_0924 <- NULL

mat_use_scaled_24_binarize <- apply(all_tme_nmf_prop_0924[,2:76], 2, function(x) ifelse(x > quantile(x, prob = .55), 1, 0)) #0.5

tissue_binarize <- mat_use_scaled_24_binarize[which(all_tme_nmf_prop_0924$tissue_upd_0924_upd == "Normal"),] #"DCIS_yes"  "synch_yes" "IDC_yes"   "Normal" 
#t_i <- 1

{
  n <- ncol(tissue_binarize)
  tissue_jacard_fea_run_24 <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      similarity <- jaccard_similarity(tissue_binarize[, i], tissue_binarize[, j])
      tissue_jacard_fea_run_24[i, j] <- similarity
      tissue_jacard_fea_run_24[j, i] <- similarity  # Symmetric matrix
    }
  }
  colnames(tissue_jacard_fea_run_24) <- rownames(tissue_jacard_fea_run_24) <- colnames(tissue_binarize)
  
  tissue_jacard_fea_run_24 <- tissue_jacard_fea_run_24[-7, -7] #remove hypoxia
  diag(tissue_jacard_fea_run_24) <- NA
  #construct the graph
  {
    tissue_jacard_fea_run_24[is.na(tissue_jacard_fea_run_24)] <- 0
    tissue_tme_nmf_24_prop_scale_cb_sel_jac_network <- graph_from_adjacency_matrix(tissue_jacard_fea_run_24, weighted=T, mode="undirected", diag=F)
    tissue_tme_nmf_24_prop_scale_cb_sel_jac_network <- delete_edges(tissue_tme_nmf_24_prop_scale_cb_sel_jac_network, E(tissue_tme_nmf_24_prop_scale_cb_sel_jac_network)[weight < 0.5]) #0.5
    #V(tissue_tme_nmf_24_prop_scale_cb_sel_jac_network)$clu <- as.character(membership(community_obj_24))
  }
  
  V(tissue_tme_nmf_24_prop_scale_cb_sel_jac_network)$clu <- as.character(membership(community_obj_24_test))
  jam_igraph(tissue_tme_nmf_24_prop_scale_cb_sel_jac_network,
             layout=l_24_test_opt, edge.color = adjustcolor('grey', alpha.f = .8), edge.width=E(tissue_tme_nmf_24_prop_scale_cb_sel_jac_network)$weight, edge_factor = 1, #edge.width=seq(0.1,5),
             edge_bundling="nodegroups",
             nodegroups=community_obj_24_info_list_test, vertex.size = 4, vertex.label=NA, vertex.color=cell_cooccur_col[V(tissue_tme_nmf_24_prop_scale_cb_sel_jac_network)$clu])
  #'#7B8C95FF'
}

legend("topright", legend = c("low", "high"), 
       lwd = c(0.5, 1),  # Line widths to represent the legend
       col = "grey",     # Color of the lines in the legend
       title = "Co-presence")


