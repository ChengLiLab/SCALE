library(Seurat)
library(ggplot2)
library(ggpubr)


# ------------------------------------------------------------------------------
# read in guided forward-selection-based Elastic Net selected aging genes 
selected_aging_genes <- readRDS("./Selected_aging_gene_sets.rds")

# read single-cell RNA-seq data
tissue_name <- "facs-Aorta"
tissue <- readRDS("./data/facs-Aorta.rds")


# ------------------------------------------------------------------------------
# calculate SCALE scores
aging_gene_down <- as.character(selected_aging_genes[[tissue_name]][1:50, "down"])
aging_gene_up <- as.character(selected_aging_genes[[tissue_name]][1:50, "up"])
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(tissue@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(tissue@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(tissue@assays$RNA@counts[aging_gene,])>0)/dim(tissue@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(tissue@assays$RNA@data[aging_gene,])),
    center = T, scale = T))

SCALE_score = t(z_matrix[aging_gene,]) %*% weights
tissue$aging_scores = SCALE_score


# ------------------------------------------------------------------------------
# visualization
df_SCALE_score <- tissue[[c("age_value", "aging_scores")]]
df_SCALE_score$aging_scores <- as.numeric(df_SCALE_score$aging_scores)
df_SCALE_score$age_value <- factor(df_SCALE_score$age_value)

ggplot(df_SCALE_score, aes(x = age_value, y = aging_scores)) +
  geom_violin(aes(fill = age_value)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("3", "18"), c("18", "24"), c("3", "24")), 
              step_increase = 0.1, map_signif_level = T, test = t.test, vjust = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  guides(fill = "none") +
  labs(x = "Chronological age (m)", y = "SCALE score")
ggsave("./figure/facs-Aorta_SCALE_score.pdf", width = 4, height = 4)
