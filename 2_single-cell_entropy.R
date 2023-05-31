.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/")

library(ggplot2)
library(gghalves)
library(cowplot)
library(Matrix)
library(Seurat)
library(ggsci)
library(ggpubr)
library(minerva)
library(energy)
library(stringr)

source("../aging.function.R")

# ------------------------------------------------------------
# Entropy
# ------------------------------------------------------------
age_data <- readLines("../age_data_all_tissue.txt")
age_data_list <- list()
for (line in age_data)
{
  context <- strsplit(line, "\t")
  age_data_list[unlist(context)[1]] <- context
}

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
  
  # calculate entropy
  df_aging_score = data.frame(age = tissue$age)
  counts.mat <- GetAssayData(tissue, "counts")
  df_aging_score$entropy = calcEntropy(counts.mat)

  # cor
  entropy_cor["all", tissue_method] <- dcor(df_aging_score$age, df_aging_score$entropy)
  
  # plot
  df_aging_score$age <- factor(tissue$age_value)
  
  ggplot(df_aging_score, aes(y=entropy, x=age, fill=age)) + geom_violin() + 
  theme_classic() + 
  scale_fill_npg() + 
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  labs(x = "Chronological age", y = "Shannon entropy", title = tissue_method) + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  stat_compare_means(method="kruskal.test", size = 5)

  ggsave(paste("../figure/shannon_entropy/", tissue_method, ".pdf", sep = ""), height = 5, width = 7)
  
  done <- done + 1
  print(done)
}