.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/tms/rds_data/")

library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggsci)
library(ggpubr)
library(stringr)
library(tidyr)
library(matrixStats)


# ------------------------------------------------------------
# (1) Calculate SCALE score for each cell type
# ------------------------------------------------------------

cell_type <- readRDS("/lustre/user/liclab/liocean/maosl/tms/cell_type.rds")
colnames(cell_type) <- c("tissue_type", "cell_type")

fail <- c()

t_test_result <- list()
cor_result <- list()

done <- 1
for (file_name in dir(".")){
  tissue_name <- unlist(str_split(file_name, pattern = ".rds"))[1]
  if(tissue_name == "facs-Brain_Non")
    tissue_name = "facs-Brain_Non-Myeloid"
  tissue <- readRDS(file_name)
  
  cell_ontology <- cell_type[cell_type$tissue_type==tissue_name, "cell_type"]
  if (length(cell_ontology) != length(levels(as.factor(tissue$cell.ontology.class)))){
    fail <- c(fail, tissue_name)
    print("failed")
    done <- done + 1
    next
  }
  
  tissue$cell_type <- factor(tissue$cell.ontology.class, labels = cell_ontology)
  
  df_aging_score <- tissue[[c("age_value", "cell_type", "aging_scores")]]
  
  cell_type_ttest_list <- list()
  cell_type_cor <- c()
  cell_type_cor["all"] <- cor(df_aging_score$age_value, df_aging_score$aging_scores) 
  
  for (cell in cell_ontology){
    df_aging_score_cell_type <- df_aging_score[df_aging_score$cell_type == cell, ]
    df_aging_score_cell_type$age_value <- factor(df_aging_score_cell_type$age_value)
    cell_type_cor[cell] <- cor(as.numeric(as.character(df_aging_score_cell_type$age_value)), 
                               df_aging_score_cell_type$aging_scores) 
    
    df_cell_type_ttest <- data.frame()
    for (age in levels(df_aging_score_cell_type$age_value)){
      if (sum(df_aging_score_cell_type$age_value==age)<2){
        next
      }
      t_test <- t.test(x = df_aging_score[df_aging_score$age_value==age, "aging_scores"],
                       y = df_aging_score_cell_type[df_aging_score_cell_type$age_value==age, "aging_scores"])
      df_cell_type_ttest[age, "t_value"] <- t_test$statistic
      df_cell_type_ttest[age, "p_value"] <- t_test$p.value  
    }
    cell_type_ttest_list[[cell]] <- df_cell_type_ttest 
    
    age_list <- levels(df_aging_score_cell_type$age_value)
    if (length(age_list)>=3){
      signif_list <- list(age_list[c(1, length(age_list))],
                          age_list[c(1, length(age_list)-1)], 
                          age_list[c(2, length(age_list))])
    }else{
      signif_list <- list(age_list[c(1, length(age_list))])
    }
    
    p <- ggplot(df_aging_score_cell_type, aes(x = age_value, y = aging_scores)) +
      geom_violin(aes(fill = age_value)) +
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
      geom_signif(comparisons = signif_list, step_increase = 0.1, map_signif_level = F, test = t.test) +
      scale_fill_npg() + 
      theme_classic() + 
      guides(fill = F) +
      labs(x="Chronological age",y="SCALE score", title = paste0(cell, " (", tissue_name, ")"))
    ggsave(p, filename = paste0("../figure/aging_score/cell_type/", tissue_name, "/", cell, ".pdf"), width = 5, height = 4)
  }
  
  t_test_result[[tissue_name]] <- cell_type_ttest_list
  cor_result[[tissue_name]] <- cell_type_cor
  print(done)
  done <- done + 1
}

saveRDS(t_test_result, "../temp_result/cell_type_t_test_result.rds")
saveRDS(cor_result, "../temp_result/cell_type_cor_result.rds")


# ------------------------------------------------------------
# (2) Correlation between chronological age and SCALE score
# ------------------------------------------------------------
cell_type <- readRDS("/lustre/user/liclab/liocean/maosl/tms/cell_type.rds")
colnames(cell_type) <- c("tissue_type", "cell_type")


done <- 1
df_cor_pvalue <- data.frame(row.names = c("tissue", "cell type", "cor_pvalue"))
fail <- c()
for (file_name in dir(".")){
  tissue_name <- unlist(str_split(file_name, pattern = ".rds"))[1]
  tissue <- readRDS(file_name)
  
  if(tissue_name == "facs-Brain_Non")
    tissue_name = "facs-Brain_Non-Myeloid"
  
  cell_ontology <- cell_type[cell_type$tissue_type==tissue_name, "cell_type"]
  if (length(cell_ontology) != length(levels(as.factor(tissue$cell.ontology.class)))){
    fail <- c(fail, tissue_name)
    print("failed")
    done <- done + 1
    next
  }
  
  tissue$cell_type <- factor(tissue$cell.ontology.class, labels = cell_ontology)
  df_aging_score <- tissue[[c("age_value", "cell_type", "aging_scores")]]
  
  cor_pvalue <- c()
  for (cell in cell_ontology){
    cor_actual <- cell_type_cor_result[[tissue_name]][cell]
    if (is.na(cor_actual)){
      cor_pvalue[cell] <- NA
      next
    }
    cell_num <- sum(df_aging_score$cell_type == cell)
    cor_record <- c()
    for (i in 1:1000){
      select <- sample(1:dim(df_aging_score)[1], cell_num, replace = T)
      sub_sample <- df_aging_score[select, ]
      cor_record <- c(cor_record, cor(sub_sample$age_value, sub_sample$aging_scores)) 
    }
    cor_record <- c(cor_record, cor_actual) 
    cor_pvalue <- c(cor_pvalue, rank(cor_record)[cell]/length(cor_record))
  }
  df_cor_pvalue <- cbind(df_cor_pvalue, 
                         rbind(rep(tissue_name, length(cor_pvalue)), cell_ontology, cor_pvalue))
  print(done)
  done <- done + 1
}

saveRDS(df_cor_pvalue, "../temp_result/cor_pvalue_cell_type.rds")

cor_pvalue_cell_type <- readRDS("/lustre/user/liclab/liocean/maosl/tms/temp_result/cor_pvalue_cell_type.rds")
cell_type_cor_result <- readRDS("/lustre/user/liclab/liocean/maosl/tms/temp_result/cell_type_cor_result.rds")
cor_pvalue_cell_type <- readRDS("/lustre/user/liclab/liocean/maosl/tms/temp_result/cor_pvalue_cell_type.rds")

cor_pvalue_cell_type["direction", ] <- 0
for(i in 1:dim(cor_pvalue_cell_type)[2]){
  cor_pvalue_cell_type[,i] <- as.character(cor_pvalue_cell_type[,i])
  cor_pvalue_cell_type["cor_pvalue",i] <- as.numeric(cor_pvalue_cell_type["cor_pvalue",i])
}
for(i in 1:dim(cor_pvalue_cell_type)[2]){
  if(!is.na(cor_pvalue_cell_type["cor_pvalue", i]) & as.numeric(cor_pvalue_cell_type["cor_pvalue", i])>0.5){
    cor_pvalue_cell_type["direction", i] <- 1
    cor_pvalue_cell_type["cor_pvalue", i] <- 1 - as.numeric(cor_pvalue_cell_type["cor_pvalue", i])
  }else{
    cor_pvalue_cell_type["direction", i] <- -1
  }
}


cor_cell_type <- data.frame()
for (tissue in names(cell_type_cor_result)){
  cor <- as.data.frame(cell_type_cor_result[[tissue]])
  cor[["cell type"]] <- row.names(cor)
  cor_pvalue <- cor_pvalue_cell_type[c("cell type", "cor_pvalue", "direction"), cor_pvalue_cell_type["tissue", ] == tissue]
  cor_pvalue <- as.data.frame(t(cor_pvalue))
  cor <- left_join(cor, cor_pvalue, by = "cell type")
  cor[["tissue"]] <- rep(tissue, dim(cor)[1])
  cor <- cor[, c(5, 2, 1, 3, 4)]
  colnames(cor) <- c("tissue", "cell type", "Pearson Correlation", "p value", "direction")
  cor$direction <- as.numeric(as.character(cor$direction))
  cor <- cor[!is.na(cor$`Pearson Correlation`), ]
  cor$`p value` <- as.numeric(as.character(cor$`p value`))
  cor$`Pearson Correlation` <- as.numeric(cor$`Pearson Correlation`)
  cor$`cell type` <- as.character(cor$`cell type`)
  cor <- cor[order(cor$`Pearson Correlation`), ]
  row.names(cor) <- 1:dim(cor)[1]
  cor$`cell type` <- ordered(cor$`cell type`, level = cor$`cell type`)
  cor[["significance"]] <- ifelse((!is.na(cor$`p value`) & (cor$`p value` < 0.05)),
                                  "+", "")
  cor[["significance"]] <- as.character(cor$significance)
  cor[(cor$direction < 0) & (cor$significance == "+"), "significance"] <- rep("-", sum((cor$direction < 0) & (cor$significance == "+")))

  ggplot(cor, aes(x = `cell type`, y = tissue)) + 
    geom_point(aes(size = `Pearson Correlation`, color = `p value`)) +
    geom_text(aes(label = cor$significance), color = "red", size = 5, vjust = 0.35) +
    #geom_text(aes(label = cor$`cell type`)) +
    scale_size(rang=c(0.1,10)) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
          axis.ticks = element_blank()) +
    labs(x = "Cell type", y = "")
  ggsave(paste0("../figure/aging_score/cell_type/", tissue, "/bubble.pdf"), width = 15, height = 7)
  
  cor_cell_type <- rbind(cor_cell_type, cor)
}
cor_cell_type[["tissue_method"]] <- cor_cell_type$tissue
cor_cell_type <- separate(cor_cell_type, "tissue", into = c("method", "tissue"), sep = "-")
cor_cell_type[cor_cell_type$tissue == "Brain_Non", "tissue"] <- "Brain_Non-Myeloid"
cor_cell_type <- cor_cell_type[order(cor_cell_type$tissue),]
cor_cell_type$`cell type` <- as.character(cor_cell_type$`cell type`)
cor_cell_type$`cell type` <- as.factor(cor_cell_type$`cell type`)
cor_cell_type$`cell type` <- relevel(cor_cell_type$`cell type`, "all")
cor_cell_type$tissue_method <- factor(cor_cell_type$tissue_method,
                                      level = c("facs-Aorta",
                                                "droplet-Fat",
                                                "facs-BAT",
                                                "facs-GAT",
                                                "facs-MAT",
                                                "facs-SCAT",
                                                "droplet-Bladder",
                                                "facs-Bladder",
                                                "facs-Brain_Myeloid",
                                                "facs-Brain_Non-Myeloid",
                                                "facs-Diaphragm",
                                                "facs-Heart",
                                                "droplet-Heart_and_Aorta",
                                                "droplet-Kidney",
                                                "facs-Kidney",
                                                "facs-Large_Intestine",
                                                "droplet-Limb_Muscle",
                                                "facs-Limb_Muscle",
                                                "droplet-Liver",
                                                "facs-Liver",
                                                "droplet-Lung",
                                                "facs-Lung",
                                                "droplet-Mammary_Gland",
                                                "facs-Mammary_Gland",
                                                "droplet-Marrow",
                                                "facs-Marrow",
                                                "droplet-Pancreas",
                                                "facs-Pancreas",
                                                "droplet-Skin",
                                                "facs-Skin",
                                                "droplet-Spleen",
                                                "facs-Spleen",
                                                "droplet-Thymus",
                                                "facs-Thymus",
                                                "droplet-Tongue",
                                                "facs-Tongue",
                                                "facs-Trachea"))


ggplot(cor_cell_type, aes(y = `cell type`, x = tissue_method)) + 
  geom_point(aes(size = `Pearson Correlation`, color = `p value`)) +
  geom_text(aes(label = cor_cell_type$significance), color = "red", size = 5, vjust = 0.5) +
  #geom_text(aes(label = cor$`cell type`)) +
  scale_size(rang=c(0.1,8)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(size="Pearson\ncorrelation", color='P value') +
  labs(x = "Cell type", y = "Tissue")
ggsave("../figure/aging_score/cell_type/all.pdf", width = 14, height = 34)

select <- names(table(cor_cell_type$`cell type`)[table(cor_cell_type$`cell type`) >= 7])
cor_cell_type <- cor_cell_type[cor_cell_type$`cell type` %in% select, ]
ggplot(cor_cell_type, aes(y = `cell type`, x = tissue_method)) + 
  geom_point(aes(size = `Pearson Correlation`, color = `p value`)) +
  geom_text(aes(label = cor_cell_type$significance), color = "red", size = 5, vjust = 0.5) +
  #geom_text(aes(label = cor$`cell type`)) +
  scale_size(rang=c(0.1,8)) +
  theme_classic() + 
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=12),
        axis.text.y = element_text(size=12)) +
  labs(size="Pearson\ncorrelation", color='P value') +
  labs(x = "Cell type", y = "Tissue")
ggsave("../figure/aging_score/cell_type/select.pdf", width = 12, height = 6)

rm(cor)
all_cor <- cor_cell_type[cor_cell_type$`cell type` == "all", c("Pearson Correlation", "tissue_method")]
cor_cell_type_cor_vs_all_cor <- data.frame()
cor <- data.frame()
for (i in 2:length(select)){
  cor_select <- cor_cell_type[cor_cell_type$`cell type` == select[i], c("Pearson Correlation", "tissue_method")]
  cor_select <- inner_join(cor_select, all_cor, "tissue_method")
  cor_cell_type_cor_vs_all_cor[i-1, "celltype"] <- select[i]
  cor_cell_type_cor_vs_all_cor[i-1, "cor"] <-  
    cor(cor_select$`Pearson Correlation.x`, cor_select$`Pearson Correlation.y`)
  colnames(cor_select) <- c("cor_cell_type", "tissue", "cor_tissue")
  cor_select[["cell_type"]] <- rep(select[i], dim(cor_select)[1])
  cor <- rbind(cor, cor_select)
}


ggplot(cor, aes(x = cor_tissue, y = cor_cell_type, color = cell_type)) +
  geom_point() +
  theme_classic() + 
  geom_smooth(method = "lm", se = F) +
  #stat_cor(method = "pearson", position = "identity") +
  scale_colour_aaas() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  labs(x="Tissue level correlation",y="Cell type level correlation", color = "Cell type")
ggsave("../figure/aging_score/cell_type/cell_type_cor_vs_all_cor.pdf", height = 4, width = 7)

saveRDS(cor_cell_type, "../temp_result/agingscore_cor_cell_type.rds")




