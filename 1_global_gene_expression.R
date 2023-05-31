.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/")

library(ggplot2)
library(Matrix)
library(Seurat)
library(stringr)
library(ggsignif)
library(ggsci)
library(latex2exp)

# ------------------------------------------------------------
# global gene expression
# ------------------------------------------------------------

age_data <- readLines("../age_data_all_tissue.txt")
age_data_list <- list()
for (line in age_data)
{
  context <- strsplit(line, "\t")
  age_data_list[unlist(context)[1]] <- context
}

done <- 1
total_num <- length(dir(".", pattern = "h5ad$"))
lm_result <- data.frame()
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
  age_tissue = as.numeric(as.character(factor(
    tissue$age, 
    levels = c(0:(age_seq_length-2)), 
    labels = age_data_list[[h5ad_file]][2:age_seq_length])))
  
  if(age_seq_length <= 2){
    print(paste("Have done:", done, "/", total_num, sep = ""))
    done <- done + 1
    next
  }
  
  # select genes that expressed in samples
  keep.genes = Matrix::rowSums(GetAssayData(tissue)>0) >= dim(tissue)[2]*0.05
  
  # gene expression level
  gene_expression <- tissue[[c("age", "nFeatures_RNA", "nCount_RNA")]]
  gene_expression$row_nCounts_per_gene <- rowMeans(t(GetAssayData(tissue, "counts")[keep.genes,]))
  gene_expression$log_row_nCounts_per_gene <- log10(gene_expression$row_nCounts_per_gene + 1)
  gene_expression$nCounts_per_gene <- gene_expression$nCount_RNA / gene_expression$nFeatures_RNA
  gene_expression$log_nCounts_per_gene <- log10(gene_expression$nCounts_per_gene + 1)
  gene_expression$age <-  factor(
    gene_expression$age, 
    levels = c(0:(age_seq_length-2)), 
    labels = age_data_list[[h5ad_file]][2:age_seq_length])
  gene_expression$log_nFeatures_RNA <- log10(gene_expression$nFeatures_RNA + 1)
  gene_expression$log_nCounts_RNA <- log10(gene_expression$nCount_RNA + 1)
  
  # linear regression
  regression <- summary(lm(data = gene_expression, log_nCounts_per_gene ~ as.numeric(age)))
  
  gene_expression <- na.omit(gene_expression)
  lm_result[tissue_method, "counts_beta"] <- regression$coefficients[2,1]
  lm_result[tissue_method, "counts_radj"] <- regression$adj.r.squared
  lm_result[tissue_method, "counts_pvalue"] <- regression$coefficients[2,4]
  
  ggplot(gene_expression, aes(x = age, y = log_nCounts_per_gene)) +
    geom_violin(aes(fill = age)) +
    stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
    theme_classic() + 
    scale_fill_npg() +
    theme_set(theme_classic(base_size = 20)) +
    guides(fill = F) +
    labs(x="Chronological age",y=paste("Global gene expression", '\n', "(log(No. of counts per gene + 1))"),
         title=tissue_method, subtitle=paste("Expr", "=",
                                             round(regression$coefficients[1,1], 3), 
                                             "+", round(regression$coefficients[2,1], 3), 
                                             "*", "age", "(", "p",
                                             ifelse(regression$coefficients[2,4]<0.001, "< 0.001", 
                                                    paste("=", round(regression$coefficients[2,4], 3))), ")"))
  
  ggsave(filename = paste0("../figure/counts_age/", tissue_method, ".pdf"), height = 5, width = 7)
  done <- done + 1
  print(paste(done, "/", total_num))
}
saveRDS(lm_result, "../temp_result/counts_vs_age_lm_result.rds")

lm_result <- readRDS("/lustre/user/liclab/liocean/maosl/tms/temp_result/counts_vs_age_lm_result.rds")
#df_mutation_rate <- read.csv("/lustre/user/liclab/liocean/maosl/tms/mutation/tms_mutation.csv", header = T, stringsAsFactors = F)
sum(lm_result$counts_beta<0)
sum(lm_result$counts_beta<0 & lm_result$counts_pvalue>0.01)
sum(lm_result$counts_beta>0 & lm_result$counts_pvalue<0.01)
p3<-ggplot(lm_result, aes(x=counts_pvalue)) +
  geom_histogram(binwidth = 0.01) +
  theme_classic() + 
  scale_fill_npg() +
  labs(x = "p value", y = "No. of tissues")

p4<-ggplot(lm_result, aes(x=counts_beta)) +
  geom_histogram() +
  theme_classic() + 
  scale_fill_npg() +
  labs(x = TeX("$\\beta$"), y = "No. of tissues")
ggsave(filename = "../figure/counts_age/different_tissue_regression_result.pdf", 
       plot = p3+p4, height = 5, width = 12)
