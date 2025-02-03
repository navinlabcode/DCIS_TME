library(Seurat)
library(dplyr)
library(magrittr)
library(RcppML)
library(future) ##run parallelization

plan("multisession", workers = 5)
options(future.globals.maxSize = 50000 * 1024^2)

##tumor NMF analysis on each tumor sample

args <- commandArgs(T)
input_file <- as.character(args[1])
#out_dir <- as.character(args[2])
setwd("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/tumor/NMF_analysis_upd/")

##example
#nohup Rscript Run_NMF.R tumor_obj_for_nmf.txt > tumor_obj_for_nmf.log 2>&1 &

# function ----------------------------------------------------------------
nmf_programs <- function(cpm, is.log=F, rank, gene_use=NULL, seed=30) {
  if(is.log==F) cpm_normalized <- log2((cpm/10) + 1) else cpm_normalized <- cpm
  if (is.null(gene_use)) {
    cpm_normalized <- cpm_normalized[apply(cpm_normalized, 1, function(x) length(which(x > 0)) > nco
l(cpm_normalized)*0.02),]
  } else {
    cpm_normalized <- cpm_normalized[rownames(cpm_normalized) %in% gene_use, ]
  }

  cpm_normalized <- cpm_normalized - rowMeans(cpm_normalized)
  cpm_normalized[cpm_normalized < 0] <- 0

  nmf_programs <- RcppML::nmf(cpm_normalized, k=rank, seed=seed)
  w_mat <- nmf_programs$w %>% set_rownames(rownames(cpm_normalized)) %>% set_colnames(paste0('NMF', 1:rank))
  h_mat <- nmf_programs$h %>% set_colnames(colnames(cpm_normalized)) %>% set_rownames(paste0('NMF', 1:rank)) %>% t
  nmf_programs_scores <- list(w_dat = w_mat, h_dat = h_mat)
  return(nmf_programs_scores)
}

# read data ---------------------------------------------------------------
input_file_list <- read.table(input_file)

tumor_count <- list()
tumor_w_list <- list()
tumor_h_list <- list()

tumor_obj_list <- list()
num <- 1
for(sample in input_file_list$V1){
  sample_name <- sub(".rds", "", sample)
  tumor_rds_raw <- readRDS(paste0("/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/tumor/refine_tumor_anno_upd/tumor_obj/", sample))
  tumor_rds_raw$sample_id <- sample_name
  
  #if(dim(tumor_rds_raw)[2] < 15) {next}
  if(dim(tumor_rds_raw)[2] < 20) {next}  

  tumor_obj_list[[num]] <- tumor_rds_raw

  message("Normalize data to CPM.....")
  tumor_rds <- CreateSeuratObject(counts = tumor_rds_raw@assays$RNA@counts, meta.data = tumor_rds_raw@meta.data)
  DefaultAssay(tumor_rds) <- "RNA"
  tumor_rds <- NormalizeData(tumor_rds, normalization.method = "RC", scale.factor = 1e6)
  tumor_count[[sample_name]] <- as.matrix(tumor_rds@assays$RNA@data)
  
  message("Run NMF.....")
  w <- NULL
  h <- NULL
  for(j in 4:20) {
    print(j)
    n <- nmf_programs(tumor_count[[sample_name]], is.log=F, rank=j)
    colnames(n$w_dat) <- paste0(sample_name, "_", j, "_", 1:j)
    colnames(n$h_dat) <- paste0(sample_name, "_", j, "_", 1:j)
    w <- cbind(w, n$w_dat)
    h <- cbind(h, n$h_dat)
  }
  tumor_w_list[[sample_name]] <- w
  tumor_h_list[[sample_name]] <- h
 
  num <- num + 1

}

saveRDS(tumor_count, "tumor_count_normalized_list_all95.rds")
saveRDS(tumor_w_list, "tumor_w_list_all95.rds")
saveRDS(tumor_h_list, "tumor_h_list_all95.rds")




