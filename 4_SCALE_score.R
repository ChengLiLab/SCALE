.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/")

library(Matrix)
library(Seurat)
library(ggplot2)
library(ggsignif)
library(ggsci)
library(stringr)
library(minerva)
library(energy)
library(nlcor)
library(Hmisc)
library(ggpubr)
library(gghalves)


# ------------------------------------------------------------
# (1) Calculate SCALE score
# ------------------------------------------------------------

age_data <- readLines("../age_data_all_tissue.txt")
age_data_list <- list()
for (line in age_data)
{
  context <- strsplit(line, "\t")
  age_data_list[unlist(context)[1]] <- context
}


aging_gene_glmnet_all_tissue <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/aging_gene_glmnet_all_tissue_add_50.rds")
aging_score_list <- list()
aging_score_pvalue <- list()
aging_score_lm <- list()
aging_score_cor <- data.frame(row.names = c("pearson", "dcor", "mic", "nlcor"))
total <- length(dir(".", pattern = "h5ad$"))
done <- 1
for (h5ad_file in dir(".", pattern = "h5ad$"))
{
  # loading datasets
  tissue_name <- unlist(strsplit(unlist(strsplit(str_extract(h5ad_file, pattern = "ns-.*\\.h5ad$"), "-"))[2], "\\."))[1]
  if (stringr::str_detect(h5ad_file, "droplet")){
    method <- "droplet" 
  }else{
    method <- "facs"
  }
  tissue_method <- paste(method, tissue_name, sep = "-")
  
  tissue <- ReadH5AD(h5ad_file, assay = "RNA", verbose = TRUE)
  
  age_seq_length <- length(age_data_list[[h5ad_file]])
  tissue$age_value = as.numeric(as.character(factor(
    tissue$age, 
    levels = c(0:(age_seq_length-2)), 
    labels = age_data_list[[h5ad_file]][2:age_seq_length])))
  
  if(age_seq_length <= 2){
    print(paste(done, "/", total, sep = ""))
    done <- done + 1
    next
  }
  
  # loading aging gene set
  aging_genes <- c(as.character(aging_gene_glmnet_all_tissue[[tissue_method]][1:50, "up"]),
                   as.character(aging_gene_glmnet_all_tissue[[tissue_method]][1:50, "down"]))
  
  # calculate SCALE score
  weights = c(rep(1, 50), rep(-1, 50))*
    Matrix::rowSums(GetAssayData(tissue, "counts")[aging_genes,]>0)/dim(GetAssayData(tissue, "counts"))[2]

  z_matrix = t(
    scale(
      t(GetAssayData(tissue, "data")[aging_genes,]),
      center = T, scale = T))
  aging.score = t(z_matrix[aging_genes,]) %*% weights  
  
  tissue$aging_scores <- aging.score
  aging_score_list[[tissue_method]] <- tissue[[c("age_value", "aging_scores")]]
  aging_score_list[[tissue_method]]$age_value <- as.factor(aging_score_list[[tissue_method]]$age_value)
  
  pvalue <- c()
  for (age in 1:(age_seq_length-2))
  {
    pvalue[paste(0, "-", age)] <- 
      t.test(tissue$aging_scores[tissue$age==0], tissue$aging_scores[tissue$age==age], alternate= "l")$p.value 
  }
  aging_score_pvalue[[tissue_method]] <- pvalue
  
  regression <- summary(lm(age_value~aging_scores, data = tissue@meta.data))
  aging_score_lm[[tissue_method]] <- c(regression$coefficients[2,4], regression$adj.r.squared)
  
  # plot
  age_list <- age_data_list[[h5ad_file]][2:age_seq_length]
  if (length(age_list)>=3){
    signif_list <- list(age_list[c(1, length(age_list))],
      age_list[c(1, length(age_list)-1)], 
      age_list[c(2, length(age_list))])
  }else{
    signif_list <- list(age_list[c(1, length(age_list))])
  }
  ggplot(aging_score_list[[tissue_method]], aes(x = age_value, y = aging_scores)) +
    geom_violin(aes(fill = age_value)) +
    stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
    geom_signif(comparisons = signif_list, 
                step_increase = 0.1, map_signif_level = T, test = t.test) +
                #step_increase = 0.14, map_signif_level = T, test = t.test, textsize = 7) +
    theme_classic() +
    scale_fill_npg() + 
    theme_set(theme_classic(base_size = 20)) +
    guides(fill = F) +
    # coord_cartesian(ylim = c(-25, 45)) +
    labs(x="Chronological age (m)",y="SCALE score", title=tissue_method)
  ggsave(paste("../figure/aging_score/larger/", tissue_method, ".pdf", sep = ""), height = 5, width = 7)
  
  print(paste(done, "/", total))
  done <- done + 1
}

saveRDS(aging_score_pvalue, "aging_score_pvalue.rds")
saveRDS(aging_score_cor, "aging_score_cor.rds")
saveRDS(aging_score_lm, "aging_score_lm.rds")

# ------------------------------------------------------------
# (2) permutation-based testing for SCALE score
# ------------------------------------------------------------

age_data <- readLines("../age_data_all_tissue.txt")
age_data_list <- list()
for (line in age_data)
{
  context <- strsplit(line, "\t")
  age_data_list[unlist(context)[1]] <- context
}
aging_gene_glmnet_all_tissue <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/aging_gene_glmnet_all_tissue_add_50.rds")
aging_score_cor <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/aging_score_cor.rds")
aging_score_lm <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/aging_score_lm.rds")
aging_score_pvalue <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/aging_score_pvalue.rds")
aging_score_random <- data.frame(row.names = c("pvalue", "regression.pval", "regression.radj", "dcor.pval", "mic.pvalue", "nlcor.pval"))
total <- length(dir(".", pattern = "h5ad$"))
done <- 1

for (h5ad_file in dir(".", pattern = "h5ad$"))
{
  age_seq_length <- length(age_data_list[[h5ad_file]])
  if(age_seq_length <= 2){
    print(paste(done, "/", total, sep = ""))
    done <- done + 1
    next
  }
  
  tissue_name <- unlist(strsplit(unlist(strsplit(str_extract(h5ad_file, pattern = "ns-.*\\.h5ad$"), "-"))[2], "\\."))[1]
  if (stringr::str_detect(h5ad_file, "droplet")){
    method <- "droplet" 
  }else{
    method <- "facs"
  }
  tissue_method <- paste(method, tissue_name, sep = "-")
  
  tissue <- ReadH5AD(h5ad_file, assay = "RNA", verbose = TRUE)
  
  tissue$age_value = as.numeric(as.character(factor(
    tissue$age, 
    levels = c(0:(age_seq_length-2)), 
    labels = age_data_list[[h5ad_file]][2:age_seq_length])))
  
  df_aging_score_random <- data.frame(row.names = 1:1000)
  counts <- GetAssayData(tissue, "counts")
  for (i in 1:1000)
  {
    aging_genes <- sample(rownames(tissue), 100)
    weights = c(rep(-1, 50), rep(1, 50))*
      rowSums(counts[aging_genes,]>0)/dim(counts)[2]
    
    z_matrix = t(
      scale(
        t(GetAssayData(tissue, "counts")[aging_genes,]),
        center = T, scale = T))
    z_matrix[is.na(t(z_matrix))[1, ], ] <- 0
    
    tissue$aging_scores <- t(z_matrix[aging_genes,]) %*% weights
    
    df_aging_score_random[i, "pearson.pval"] <- cor(tissue$age_value, tissue$aging_scores)
    df_aging_score_random[i, "mic.pval"] <- cstats(as.matrix(as.numeric(as.character(tissue$age_value))), as.matrix(tissue$aging_scores), alpha = 0.6)[3]

    if(i %% 100==0){
      print(i)
      flush.console()
    }
  }
 
  aging_score_random["pearson.pval", tissue_method] <-
    (1001 - rank(c(df_aging_score_random$pearson.pval, aging_score_cor["pearson", tissue_method]))[1001])/1001
  aging_score_random["mic.pval", tissue_method] <-
    (1001 - rank(c(df_aging_score_random$mic.pval, aging_score_cor["mic", tissue_method]))[1001])/1001
 
  print(paste(done, "/", total))
  done <- done + 1
}


# ------------------------------------------------------------
# (3) Mutation analysis
# ------------------------------------------------------------

setwd("/lustre/user/liclab/liocean/maosl/tms/rds_data/")
fail <- c()

t_test_result <- list()
cor_result <- list()
df_mutation_rate <- read.csv("/lustre/user/liclab/liocean/maosl/tms/mutation/tms_mutation.csv",
                             header = T,
                             stringsAsFactors = F)
mutation_vs_residual_pval <- c()
done <- 1
for (file_name in dir(".", pattern = "facs-*")){
  tissue_method <- unlist(str_split(file_name, pattern = ".rds"))[1]
  if(tissue_method == "facs-Brain_Non")  tissue_method = "facs-Brain_Non-Myeloid"
  tissue_name <- unlist(str_split(tissue_method, pattern = "facs-"))[2]
  tissue <- readRDS(file_name)
  
  tissue_mutation_rate <- df_mutation_rate[df_mutation_rate$tissue == tissue_name, c("cell", "adata")]
  
  cell_code <- unlist(as.character(tissue_mutation_rate$cell))
  
  df_tissue_age_info <- tissue[[c("cell", "age_value", "aging_scores")]]
  colnames(df_tissue_age_info) <- c("cell", "age", "aging_scores")
  for (cell_index in 1:length(tissue$cell)){
    if (sum(str_detect(cell_code, pattern = paste("^", tissue$cell[cell_index], sep = "")))){
      df_tissue_age_info[cell_index, "mutation_rate"] <- 
        tissue_mutation_rate[
          which(str_detect(cell_code, pattern = paste("^", tissue$cell[cell_index], sep = ""))),
          "adata"]
    }else{
      df_tissue_age_info[cell_index, "mutation_rate"] <- NaN
    }
  }
  df_tissue_age_info <- na.omit(df_tissue_age_info)
  
  df_tissue_age_info$age <- as.factor(df_tissue_age_info$age)
  
  ggplot(df_tissue_age_info, aes(x = aging_scores, y = mutation_rate)) + 
    geom_point(aes(color = age)) +
    geom_smooth(method=lm, se=F, color = "black") +
    ggpubr::stat_cor(method = "pearson", size = 4.5) +
    theme_classic() +
    scale_color_npg() +
    theme_set(theme_classic(base_size = 20)) +
    labs(x="SCALE score",y="Mean number of somatic mutations in genes", col="Age (m)", title=tissue_method)
  ggsave(paste("/lustre/user/liclab/liocean/maosl/tms/figure/mutation_agingscore_dot/", tissue_method, ".pdf", sep = ""),
         height = 4, width = 6)
  
  regression <- lm(data = df_tissue_age_info, formula = aging_scores ~ age)
  df_tissue_age_info$aging_scores_residuals <- regression$residuals
  
  mutation_vs_residual_pval[tissue_method] <- 
    summary(lm(data = df_tissue_age_info, formula = mutation_rate ~ aging_scores_residuals))$coefficients[2, 4]
  
  ggplot(df_tissue_age_info, aes(x = aging_scores_residuals, y = mutation_rate)) + 
    geom_point(aes(color = age)) +
    geom_smooth(method=lm, se = F, color = "black") +
    ggpubr::stat_cor(method = "pearson", size = 4.5) +
    theme_classic() +
    scale_color_npg() +
    theme_set(theme_classic(base_size = 20)) +
    labs(x="SCALE score residuals",y="Mean number of somatic mutations in genes", col="Age (m)", title=tissue_method)
  ggsave(paste("/lustre/user/liclab/liocean/maosl/tms/figure/mutation_deltaagingscore_dot/", tissue_method, ".pdf", sep = ""),
         height = 4, width = 6)
  
  done <- done + 1
  print(done)
  
}

names(mutation_vs_residual_pval)[mutation_vs_residual_pval < 0.05]
sum(mutation_vs_residual_pval >= 0.05)
