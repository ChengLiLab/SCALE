.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/")
source("../aging.function.R")

library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(Matrix)
library(Seurat)
library(glmnet)
library(caTools)
library(stringr)

# ------------------------------------------------------------
# (1) Summarize age-related genes
# ------------------------------------------------------------
capitalize <- function(x) {
  x = tolower(x)
  substr(x, 1, 1) = toupper(substr(x, 1, 1))
  return(x)
}

# manual_check_aging_gene
af_mc_genage <- read.table(
  "/lustre/user/liclab/liocean/maosl/tms/manual_select_gene/manual_check_mouse_csgene_result.txt",
  header = TRUE, sep = "\t")
af_mc_genage$Symbol = sapply(af_mc_genage$Symbol, capitalize)

# GenAge
af_genage = read.csv(
  "/lustre/user/liclab/licontinent/sujy/Thesis/AgingFeatures/GenAge/genage_mice.csv",
  header = T, stringsAsFactors = F)
af_genage$Symbol = sapply(af_genage$symbol, capitalize)

# CSGene
af_csgene = read.table(
  "/lustre/user/liclab/licontinent/sujy/Thesis/AgingFeatures/csgene_mouse.txt", 
  header = T, stringsAsFactors = F)
af_csgene$Symbol = sapply(af_csgene$MouseGeneSymbol, capitalize)

# AgeFactDB
af_agefactdb = read.csv(
  "/lustre/user/liclab/licontinent/sujy/Thesis/AgingFeatures/AgeFactDB/afdb_mice.csv", 
  header = T, stringsAsFactors = F)
af_agefactdb$Symbol = sapply(af_agefactdb$Gene.Symbol..Homology.Group.Gene., capitalize)

# agingatlas
af_agingatlas <- read.table("../AgingAtlas_mouse.txt",
                            header = T, stringsAsFactors = F, sep = ",")
af_agingatlas$Symbol = sapply(af_agingatlas$Symbol, capitalize)

# GO: Aging
af_go_aging = read.table(
  "/lustre/user/liclab/licontinent/sujy/Thesis/AgingFeatures/GO/GO_Aging_mice.txt",
  header = F, sep = "\t", stringsAsFactors = F)
colnames(af_go_aging) = 
  c("GeneID","Symbol","FullName","Organism","GOClass","Annotation",
    "Evidence","EvidenceWith","Type","Isoform","Refence")
af_go_aging$Symbol = sapply(af_go_aging$Symbol, capitalize)

# GO: CellAge
af_go_ca = read.table(
  "/lustre/user/liclab/licontinent/sujy/Thesis/AgingFeatures/GO/GO_CellAge_mice.txt",
  header = F, sep = "\t", stringsAsFactors = F)
colnames(af_go_ca) = 
  c("GeneID","Symbol","FullName","Organism","GOClass","Annotation",
    "Evidence","EvidenceWith","Type","Isoform","Refence")
af_go_ca$Symbol = sapply(af_go_ca$Symbol, capitalize)

# GO: CellularSenescence
af_go_cs = read.table(
  "/lustre/user/liclab/licontinent/sujy/Thesis/AgingFeatures/GO/GO_CellularSenescence_mice.txt",
  header = F, sep = "\t", stringsAsFactors = F)
colnames(af_go_cs) = 
  c("GeneID","Symbol","FullName","Organism","GOClass","Annotation",
    "Evidence","EvidenceWith","Type","Isoform","Refence")
af_go_cs$Symbol = sapply(af_go_cs$Symbol, capitalize)

# GO: RegulationOfCellAging
af_go_roca = read.table(
  "/lustre/user/liclab/licontinent/sujy/Thesis/AgingFeatures/GO/GO_RegulationOfCellAging_mice.txt",
  header = F, sep = "\t", stringsAsFactors = F)
colnames(af_go_roca) = 
  c("GeneID","Symbol","FullName","Organism","GOClass","Annotation",
    "Evidence","EvidenceWith","Type","Isoform","Refence")
af_go_roca$Symbol = sapply(af_go_roca$Symbol, capitalize)

# GO: ReplicativeSenescence
af_go_rs = read.table(
  "/lustre/user/liclab/licontinent/sujy/Thesis/AgingFeatures/GO/GO_ReplicativeAge_mice.txt",
  header = F, sep = "\t", stringsAsFactors = F)
colnames(af_go_rs) = 
  c("GeneID","Symbol","FullName","Organism","GOClass","Annotation",
    "Evidence","EvidenceWith","Type","Isoform","Refence")
af_go_rs$Symbol = sapply(af_go_rs$Symbol, capitalize)

db_aging_gene <- c(af_genage$Symbol,
                      af_csgene$Symbol,
                      af_agefactdb$Symbol,
                      af_agingatlas$Symbol,
                      af_go_aging$Symbol,
                      af_go_ca$Symbol,
                      af_go_cs$Symbol,
                      af_go_roca$Symbol,
                      af_go_rs$Symbol)
diff_literature_db <- setdiff(af_mc_genage$Symbol, db_aging_gene)
writeLines(diff_literature_db, "./different_genes_literature_db.txt")

# overlap
library(VennDiagram)
png("../figure/summarize aging genes/venn_fig.png", height = 600, width = 600)
v1 = venn.diagram(
  list("GenAge"=af_genage$Symbol,
       "CSGene"=af_csgene$Symbol,
       "Literature"=af_mc_genage$Symbol,
       "GO: Aging"=af_go_aging$Symbol),
  filename = NULL,
  cex = 2,
  cat.cex = 2,
  alpha=0.4,
  cat.fontfamily = "Arial",
  fontfamily = "Arial",
  fill = c("LightSalmon", "LightYellow", "PaleGreen", "Lavender"))
grid.draw(v1)
dev.off()

library(UpSetR)
pdf("../figure/summarize aging genes/upsetr_fig.pdf", width = 10, height = 6)
listInput = list("GenAge"=af_genage$Symbol,
                 "CSGene"=af_csgene$Symbol,
                 "Literature based aging genes"=af_mc_genage$Symbol,
                 "AgeFactDB"=af_agefactdb$Symbol,
                 "AgingAtlas"=af_agingatlas$Symbol,
                 "GO: Aging"=af_go_aging$Symbol,
                 "GO: Cell Age"=af_go_ca$Symbol,
                 "GO: Cellular Senescence"=af_go_cs$Symbol,
                 "GO: Regulation of Cell Aging"=af_go_roca$Symbol,
                 "GO: Replicative Senescence"=af_go_rs$Symbol)
upset(fromList(listInput), order.by = "freq", 
      nsets = 8, text.scale = 1.5,
      matrix.color = "dodgerblue4",
      main.bar.color = "Darkmagenta",
      sets.bar.color = "Darkorange2")
dev.off()

af_merge = unique(c(
  af_agefactdb$Symbol,
  af_csgene$Symbol,
  af_mc_genage$Symbol,
  af_genage$Symbol,
  af_go_aging$Symbol,
  af_go_ca$Symbol,
  af_go_cs$Symbol,
  af_go_roca$Symbol,
  af_go_rs$Symbol
))

af_merge_no_go = unique(c(
  af_agefactdb$Symbol,
  af_csgene$Symbol,
  af_mc_genage$Symbol,
  af_genage$Symbol
))

write.table(
  af_merge, 
  "/lustre/user/liclab/liocean/maosl/tms/af_merge_mice.txt", 
  quote = F, row.names = F, col.names = F)

# GO analysis
library(clusterProfiler)
library(topGO) 
library(Rgraphviz)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GOplot)

# for all genes collected from literature and databases
# gene name conversion
sample.genename.converted <- mapIds(
  org.Mm.eg.db, af_merge, 
  keytype ="SYMBOL", column="ENTREZID",
  multiVals = "first")
df_convert <- data.frame(unlist(sample.genename.converted))
colnames(df_convert) <- "entrez_id"
convert_fail <- row.names(df_convert)[is.na(df_convert$entrez_id)]
sample.gene.entrezid <- df_convert[!is.na(df_convert$entrez_id),]

# GO biological process enrichment
sample.enrichgo.bp = enrichGO(gene = sample.gene.entrezid,
                              OrgDb = org.Mm.eg.db,
                              keyType = "ENTREZID",
                              ont = "BP",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)

pdf("../figure/GO_analysis/database_literatures_GO_dotplot.pdf", width = 7.5, height = 5)
dotplot(sample.enrichgo.bp, showCategory=10, 
        title="Aging related genes (Databases & Literatures) \nEnrichment GO Biological Processes (Mouse)")
dev.off()

# for genes collected from literature
sample.genename.converted <- mapIds(
  org.Mm.eg.db, af_mc_genage$Symbol, 
  keytype ="SYMBOL", column="ENTREZID",
  multiVals = "first")
df_convert <- data.frame(unlist(sample.genename.converted))
colnames(df_convert) <- "entrez_id"
convert_fail <- row.names(df_convert)[is.na(df_convert$entrez_id)]
sample.gene.entrezid <- df_convert[!is.na(df_convert$entrez_id),]

# GO biological process enrichment
sample.enrichgo.bp = enrichGO(gene = sample.gene.entrezid,
                              OrgDb = org.Mm.eg.db,
                              keyType = "ENTREZID",
                              ont = "BP",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)

pdf("../figure/GO_analysis/literatures_GO_dotplot.pdf", width = 6, height = 3.5)
#dotplot(sample.enrichgo.bp, showCategory=10, font.size = 25, title="Aging-related genes (Literature)")
emapplot(sample.enrichgo.bp, color = "p.adjust", showCategory = 10, repel = TRUE,
         layout = "kk", title="Aging-related genes (Literature)")
dev.off()

# for genes collected from databases
sample.genename.converted <- mapIds(
  org.Mm.eg.db, af_merge_no_go, 
  keytype ="SYMBOL", column="ENTREZID",
  multiVals = "first")
df_convert <- data.frame(unlist(sample.genename.converted))
colnames(df_convert) <- "entrez_id"
convert_fail <- row.names(df_convert)[is.na(df_convert$entrez_id)]
sample.gene.entrezid <- df_convert[!is.na(df_convert$entrez_id),]

# GO biological process enrichment
sample.enrichgo.bp = enrichGO(gene = sample.gene.entrezid,
                              OrgDb = org.Mm.eg.db,
                              keyType = "ENTREZID",
                              ont = "BP",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)

pdf("../figure/GO_analysis/database_no_GO_literatures_GO_dotplot.pdf", width = 7.5, height = 5)
dotplot(sample.enrichgo.bp, showCategory=10, 
        title="Aging related genes (Databases & Literatures) \nEnrichment GO Biological Processes (Mouse)")
dev.off()

# ------------------------------------------------------------
# (2) Check expression in TMS datasets
# ------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(ggsci)

.findProp <- function(A, B){
  A = unique(A)
  B = unique(B)
  return(as.vector(table(A %in% B)/length(A)))
}

findExpProp <- function(obj, threshold=0.05){
  require(Seurat)
  # filter genes
  counts = GetAssayData(obj, "data") # use unscaled raw data
  keep.genes = Matrix::rowSums(counts>0) >= dim(obj)[2]*threshold
  counts = counts[keep.genes,]
  
  # calculate expression proportion
  df_af_expr = data.frame(
    "GenAge"=.findProp(af_genage$Symbol, row.names(counts)),
    "CSGene"=.findProp(af_csgene$Symbol, row.names(counts)),
    "Selected_aging_gene"=.findProp(af_mc_genage$Symbol, row.names(counts)),
    "AgeFactDB"=.findProp(af_agefactdb$Symbol, row.names(counts)),
    "GO_Aging"=.findProp(af_go_aging$Symbol, row.names(counts)),
    "GO_CellAge"=.findProp(af_go_ca$Symbol, row.names(counts)),
    "GO_Cellular_Senescence"=.findProp(af_go_cs$Symbol, row.names(counts)),
    "GO_Regulation_of_Cell_Aging"=.findProp(af_go_roca$Symbol, row.names(counts)),
    "GO_Replicative_Senescence"=.findProp(af_go_rs$Symbol, row.names(counts)),
    "Expr_Gene_Num"=as.vector(table(keep.genes))
  )
  rownames(df_af_expr) = c("Unexpressed", "Expressed")

  
  return (df_af_expr)
}

record_gene_expr <- function(obj, gene_list, threshold=0.05)
{
  require(Seurat)
  # filter genes
  counts = GetAssayData(obj, "data") # use unscaled raw data
  keep.genes = Matrix::rowSums(counts>0) >= dim(obj)[2]*threshold
  counts = counts[keep.genes,]
  return(gene_list %in% unique(row.names(counts)))
}

# ------------------------------------------------------------
# (2.1) Check expression in all tissues
# ------------------------------------------------------------
library(stringr)
# Sex: 0 female, 1 male
# Age: 0 3m, 1 24m
# Annotation: 
plot_af_expr_all_tissue <- list()
merge_af_expr_all_tissue <- 
  data.frame(row.names = c("GenAge","CSGene",
                           "Selected_aging_gene",
                           "AgeFactDB",
                           "GO_Aging",
                           "GO_CellAge",
                           "GO_Cellular_Senescence",
                           "GO_Regulation_of_Cell_Aging",
                           "GO_Replicative_Senescence",
                           "Expr_gene_num"))
df_mc_aging_gene_expr <- data.frame(row.names = af_mc_genage$Symbol)
df_merge_aging_gene_expr <- data.frame(row.names = af_merge)
total <- length(dir(".", pattern = "h5ad$"))
done <- 1
for (h5ad_file in dir(".", pattern = "h5ad$"))
{
  tissue <- ReadH5AD(h5ad_file,
                     assay = "RNA", 
                     verbose = TRUE)
  tissue_name <- unlist(strsplit(unlist(strsplit(str_extract(h5ad_file, pattern = "ns-.*\\.h5ad$"), "-"))[2], "\\."))[1]
  if (stringr::str_detect(h5ad_file, "droplet")){
    method <- "droplet" 
  }else{
    method <- "facs"
  }
  tissue_name <- paste(method, tissue_name, sep = '-')
  df_af_expr_tissue = findExpProp(tissue)
  merge_af_expr_all_tissue[tissue_name] <- t(df_af_expr_tissue["Expressed",])
  df_mc_aging_gene_expr[tissue_name] <- record_gene_expr(obj=tissue, gene_list=af_mc_genage$Symbol) 
  df_merge_aging_gene_expr[tissue_name] <- record_gene_expr(obj=tissue, gene_list=af_merge) 
  df_vs = reshape2::melt(t(df_af_expr_tissue[, 1:9]))
  colnames(df_vs) = c("source","expression", "value")
  plot_af_expr_all_tissue[[tissue_name]] = 
    ggplot(df_vs, aes(fill=expression, y=value, x=source)) + 
    geom_bar(position="fill", stat="identity") + 
    labs(x = "", y="Proportion", title=tissue_name) + 
    theme_classic() +
    theme(
      axis.text.x=element_text(angle=45, size=10, hjust = 1, color="black"),
      axis.text.y=element_text(size=10, color = "black"),
      legend.title = element_blank(),
      legend.position = "top") +
    scale_fill_manual(values = c("Expressed"="red","Unexpressed"="light gray"))
  print(paste(tissue_name, done, "/", total))
  done <- done + 1
}

# save data
saveRDS(merge_af_expr_all_tissue, "af_expr_prop_all_tissue.rds")
saveRDS(df_mc_aging_gene_expr, "select_aging_gene_expr_all_tissue.rds")
saveRDS(plot_af_expr_all_tissue, "figure_af_expr_all_tissue.rds")
saveRDS(df_merge_aging_gene_expr, "merge_aging_gene_expr.rds")

# ------------------------------------------------------------
# (3) Select aging-related genes
# ------------------------------------------------------------

# collected aging genes
af_merge = readLines("../af_merge_mice.txt")

age_data <- readLines("../age_data_all_tissue.txt")
age_data_list <- list()
for (line in age_data)
{
  context <- strsplit(line, "\t")
  age_data_list[unlist(context)[1]] <- context
}

# ------------------------------------------------------------
# (3.1) Perform Elastic-Net regression on all genes and CA
# ------------------------------------------------------------

glm_cor_list <- list()
glm_cor_df <- data.frame()
aging_gene_glmnet_all_tissue <- list()
total_num <- length(dir(".", pattern = "h5ad$"))
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
  tissue_method <- paste0(method, "-", tissue_name)
  tissue_glm <- ReadH5AD(h5ad_file, assay = "RNA", verbose = TRUE)
  
  # Filter genes that do not express
  keep.genes = rowSums(GetAssayData(tissue_glm)>0) >= dim(tissue_glm)[2]*0.05
  
  table(keep.genes)
  
  tissue_glm_counts.hq = GetAssayData(tissue_glm, "counts")[keep.genes,]
  age_seq_length <- length(age_data_list[[h5ad_file]])
  age_tissue = as.numeric(as.character(factor(
    tissue_glm$age, 
    levels = c(0:(age_seq_length-2)), 
    labels = age_data_list[[h5ad_file]][2:age_seq_length])))
  
  if(age_seq_length <= 2){
    print(paste("Have done:", done, "/", total_num, sep = ""))
    done <- done + 1
    next
  }
  
  # Prepare training and test data
  set.seed(0)
  split = sample.split(tissue_glm$age, SplitRatio = 0.67)
  x.train = tissue_glm_counts.hq[,split == T]
  y.train = age_tissue[split == T]
  x.test = tissue_glm_counts.hq[,split == F]
  y.test = age_tissue[split == F]
  
  df_glm_cor <- data.frame(row.names = c("alpha", "train", "test", "mse_train", "mse_test", "lambda.1se"))
  lambda.train.record <- c()
  for (alpha in 1:10)
  {
    df_glm_cor["alpha", alpha] = alpha/10
    # Run glm
    glmnet.train.CV = cv.glmnet(t(x.train), y.train, nfolds=10, 
                                family="gaussian",
                                alpha=alpha/10)
    lambda.train = glmnet.train.CV$lambda.1se
    lambda.train.record <- c(lambda.train.record, lambda.train)
    glmnet.train = glmnet(
      t(x.train), y.train, family="gaussian", alpha=alpha/10, lambda=glmnet.train.CV$lambda.1se)
    df_glm_cor["lambda.1se", alpha] = lambda.train
    
    # Check peformance on training set
    y.train.est = predict(
      glmnet.train, t(x.train), type="response", s=lambda.train)
    df_glm_cor["train", alpha] = cor(y.train, y.train.est)
    df_glm_cor["mse_train", alpha] = (mean((y.train - y.train.est)**2))**(1/2)
    
    # Check performance on test set
    y.test.est <- predict(
      glmnet.train, t(x.test), type="response", s=lambda.train)
    df_glm_cor["test", alpha] = cor(y.test, y.test.est)
    df_glm_cor["mse_test", alpha] = (mean((y.test - y.test.est)**2))**(1/2)
  }
  glm_cor_list[[paste(method, tissue_name, sep = '-')]] <- df_glm_cor
  
  # min mse model
  min_mse_model <- which(df_glm_cor["mse_test",]==min(df_glm_cor["mse_test",]))
  glmnet.train.result = glmnet(
    t(x.train), y.train, family="gaussian", alpha=min_mse_model/10, lambda=lambda.train.record[min_mse_model])
 
   # Calculate weighted average expression
  coef = sort(coef(glmnet.train.result)[,which(glmnet.train.result$lambda == lambda.train.record[min_mse_model])])
  y.train.est = predict(
    glmnet.train.result, t(x.train), type="response", s=lambda.train.record[min_mse_model])
  glm_cor_df[tissue_method, "train"] <- cor(y.train, y.train.est)
  
  df.train = data.frame(
    v1=y.train, v2=y.train.est, re=y.train.est-y.train, 
    sex = factor(tissue_glm$sex[split == T], levels = 0:1, labels = c("female", "male")))
  colnames(df.train) <- c("v1", "v2", "re", "sex")
  
  p1 <- ggplot(df.train, aes(x=v1, y=v2, color=sex)) +
    geom_point(size=0.3, alpha=0.5, position = "jitter") +
    geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1) +
    theme_classic() +
    labs(x = "Chronological age (m)",
         y = "Estimated age (m)",
         title = "Training set") +
    guides(colour=guide_legend(override.aes = list(size=5)))+
    scale_color_aaas()
  
  y.test.est <- predict(
    glmnet.train.result, t(x.test), type="response", s=lambda.train.record[min_mse_model])
  glm_cor_df[tissue_method, "test"] <- cor(y.test, y.test.est) # 0.708
  
  df.test = data.frame(
    v1=y.test, v2=y.test.est, re=y.test.est-y.test, 
    sex = factor(tissue_glm$sex[split == F], levels = 0:1, labels = c("female", "male")))
  colnames(df.test) <- c("v1", "v2", "re", "sex")
  
  p2 <- ggplot(df.test, aes(x=v1, y=v2, color=sex)) +
    geom_point(size=0.3, alpha=0.5, position = "jitter") +
    geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1) +
    theme_classic() +
    labs(x = "Chronological age (m)",
         y = "Estimated age (m)",
         title = "Test set") +
    guides(colour=guide_legend(override.aes = list(size=5)))+
    scale_color_aaas()
  
  cowplot::plot_grid(p1, p2, ncol = 2)
  ggsave(paste0("../figure/glmnet/regression_result_selected_genes/", tissue_method, ".pdf"), height = 3, width = 8)
  
  # record aging gene selected by glmnet
  aging_genes = data.frame(
    names(head(coef, 250)),
    names(rev(tail(coef, 251))[2:251]) # remove intercept
  )
  colnames(aging_genes) <- c("up", "down")
  aging_gene_glmnet_all_tissue[[paste(method, tissue_name, sep = '-')]] <- aging_genes
  print(paste("Have done:", done, "/", total_num, sep = ""))
  done <- done + 1
}
saveRDS(glm_cor_df, "glm_cor_manual_selected_gene.rds")
saveRDS(glm_cor_list, "glm_cor_list_all_gene.rds")
saveRDS(aging_gene_glmnet_all_tissue, "aging_gene_glmnet_all_tissue_500.rds")

# ------------------------------------------------------------
# (3.2) Add genes identified in (3.1)
# ------------------------------------------------------------

total_num <- length(dir(".", pattern = "h5ad$"))

for (add in seq(0, 300, 10))
{
  print(paste("add:", add))
  done <- 1
  aging_gene_glmnet_all_tissue_add <- list()
  df_glm_cor <- data.frame(row.names = c("train", "test", "mse_train", "mse_test", "lambda.1se"))
  for (h5ad_file in dir(".", pattern = "h5ad$"))
  {
    tissue_name <- unlist(strsplit(unlist(strsplit(str_extract(h5ad_file, pattern = "ns-.*\\.h5ad$"), "-"))[2], "\\."))[1]
    if (stringr::str_detect(h5ad_file, "droplet")){
      method <- "droplet" 
    }else{
      method <- "facs"
    }
    tissue_method <- paste(method, tissue_name, sep = '-')
    
    tissue_glm <- ReadH5AD(h5ad_file, assay = "RNA", verbose = TRUE)
    
    # Filter counts
    aging_gene_set <- unique(c(af_merge, 
                               as.character(aging_gene_glmnet_all_tissue[[tissue_method]][1:add, 1]), 
                               as.character(aging_gene_glmnet_all_tissue[[tissue_method]][1:add, 2]))) 
    keep.genes = Matrix::rowSums(GetAssayData(tissue_glm)>0) >= dim(tissue_glm)[2]*0.05 & rownames(tissue_glm) %in% aging_gene_set
    
    tissue_glm_counts.hq = GetAssayData(tissue_glm, "counts")[keep.genes,]
    age_seq_length <- length(age_data_list[[h5ad_file]])
    age_tissue = as.numeric(as.character(factor(
      tissue_glm$age, 
      levels = c(0:(age_seq_length-2)), 
      labels = age_data_list[[h5ad_file]][2:age_seq_length])))
    
    if(age_seq_length <= 2){
      print(paste("Have done:", done, "/", total_num, sep = ""))
      done <- done + 1
      next
    }
    
    # Prepare training and test data
    set.seed(0)
    split = sample.split(tissue_glm$age, SplitRatio = 0.67)
    x.train = tissue_glm_counts.hq[,split == T]
    y.train = age_tissue[split == T]
    x.test = tissue_glm_counts.hq[,split == F]
    y.test = age_tissue[split == F]
    
    glmnet.train.CV = cv.glmnet(t(x.train), y.train, nfolds=10, 
                                family="gaussian",
                                alpha=0.5)
    lambda.train = glmnet.train.CV$lambda.1se
    glmnet.train = glmnet(
      t(x.train), y.train, family="gaussian", alpha=0.5, lambda=lambda.train)
    df_glm_cor["lambda.1se", tissue_method] = lambda.train
    
    # Check peformance on training set
    y.train.est = predict(
      glmnet.train, t(x.train), type="response", s=lambda.train)
    df_glm_cor["train", tissue_method] = cor(y.train, y.train.est)
    df_glm_cor["mse_train", tissue_method] = (mean((y.train - y.train.est)**2))**(1/2)
    
    # Check performance on test set
    y.test.est <- predict(
      glmnet.train, t(x.test), type="response", s=lambda.train)
    df_glm_cor["test", tissue_method] = cor(y.test, y.test.est)
    df_glm_cor["mse_test", tissue_method] = (mean((y.test - y.test.est)**2))**(1/2)
    
    # Calculate weighted average expression
    coef = sort(coef(glmnet.train)[,which(glmnet.train$lambda == lambda.train)])
    
    # record aging gene selected by glmnet
    aging_genes = data.frame(
      names(head(coef, 150)),
      names(rev(tail(coef, 151))[2:151])) # remove intercept
    
    colnames(aging_genes) <- c("down", "up")
    aging_gene_glmnet_all_tissue_add[[tissue_method]] <- aging_genes
    print(table(c(as.character(aging_genes[,1]), as.character(aging_genes[,2])) %in% af_merge))
    print(paste("Have done:", done, "/", total_num, sep = ""))
    flush.console()
    done <- done + 1
  }
  
  saveRDS(aging_gene_glmnet_all_tissue_add, paste("aging_gene_glmnet_all_tissue_add_", add, ".rds", sep = ""))
  saveRDS(df_glm_cor, paste("glm_res_add", add, ".rds", sep = ""))
  
}


# ------------------------------------------------------------
# (3.3) Mse change comparison
# ------------------------------------------------------------

# loading data produced in (3.2)
glm_res_add10 <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/glm_res_add10.rds")
glm_res_add20 <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/glm_res_add20.rds")
glm_res_add30 <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/glm_res_add30.rds")
glm_res_add40 <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/glm_res_add40.rds")
glm_res_add50 <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/glm_res_add50.rds")
glm_res_add100 <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/glm_res_add100.rds")
glm_res_add150 <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/glm_res_add150.rds")
glm_res_add200 <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/glm_res_add200.rds")
glm_res_add250 <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/glm_res_add250.rds")

glm_result_msetest["add10", ] <- glm_res_add10["mse_test", ]
glm_result_msetest["add20", ] <- glm_res_add20["mse_test", ]
glm_result_msetest["add30", ] <- glm_res_add30["mse_test", ]
glm_result_msetest["add40", ] <- glm_res_add40["mse_test", ]
glm_result_msetest["add50", ] <- glm_res_add50["mse_test", ]
glm_result_msetest["add100", ] <- glm_res_add100["mse_test", ]
glm_result_msetest["add150", ] <- glm_res_add150["mse_test", ]
glm_result_msetest["add200", ] <- glm_res_add200["mse_test", ]
glm_result_msetest["add250", ] <- glm_res_add250["mse_test", ]

glm_result_msetest <- glm_result_msetest[c(2:11, 1),]
plot_glm_mse_data <- reshape2::melt(t(glm_result_msetest))
colnames(plot_glm_mse_data) <- c("tissue", "add", "mse")
pdf("../figure/glmnet/mse_add_genes.pdf", width = 10, height = 5)
ggplot(data = plot_glm_mse_data, aes(x = tissue, y = mse, fill = add)) +
  geom_bar(position="dodge", stat="identity") + 
  labs(x = "Tissue", y="Mean square error", title="") + 
  theme_classic() + 
  theme(
    axis.text.x=element_text(angle=45, size=10, hjust = 1, color="black"),
    axis.text.y=element_text(size=10, color = "black"),
    legend.title = element_blank())
dev.off()
