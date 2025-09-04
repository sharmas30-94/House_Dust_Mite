# Rnbeads Methylation analysis
BiocManager::install("RnBeads")
library("RnBeads")
BiocManager::install("RnBeads.mm10")
library("RnBeads.mm10")
BiocManager::install("wateRmelon")
BiocManager::install("org.Mm.eg.db")
library("wateRmelon")
library("org.Mm.eg.db")
library("wordcloud")
# Tumor Tissue
data.dir <- "/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Tumor_Tissue"
idat.dir <- file.path(data.dir, "idat")
sample.annotation <- file.path(data.dir, "sample_annotation.csv")
analysis.dir <- "/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Tumor_Tissue/analysis"
report.dir <- file.path(analysis.dir, "reports")
rnb.options(filtering.sex.chromosomes.removal=TRUE, identifiers.column="Sample_ID",normalization.method = "wm.dasen", normalization.background.method = "methylumi.noob", differential.enrichment.go = FALSE, assembly="mm10", differential.comparison.columns = "Treatment", differential.report.sites = TRUE, export.to.csv = TRUE, export.types=c("tiling", "genes", "promoters", "cpgislands"))
tumor <- rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation, data.dir=idat.dir, data.type="infinium.idat.dir")

save.rnb.set(tumor, , archive = TRUE)

## normal Tissue
data.dir <- "/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Normal_Tissue"
idat.dir <- file.path(data.dir, "idat")
sample.annotation <- file.path(data.dir, "sample_annotation.csv")
analysis.dir <- "/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Normal_Tissue/analysis"
report.dir <- file.path(analysis.dir, "reports")
rnb.options(filtering.sex.chromosomes.removal=TRUE, identifiers.column="Sample_ID",normalization.method = "wm.dasen", normalization.background.method = "methylumi.noob", assembly="mm10", differential.comparison.columns = "Treatment", differential.enrichment.go = FALSE,  differential.report.sites = TRUE, export.to.csv = TRUE, export.types=c("tiling", "genes", "promoters", "cpgislands"))
normal <- rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation, data.dir=idat.dir, data.type="infinium.idat.dir")
save.rnb.set(normal, path= "/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Normal_Tissue/rnbeads_normal_object/", archive = TRUE)




### Get Methylation Analysis ####

library(ggplot2)
library(reshape2)
df <- read.csv("/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Differential_sites.csv", header=TRUE, sep=",")
df_filtered <- df[df$Region %in% c("Genes", "CpG islands", "Promoters", "Tiling"), ]

# Melt for ggplot2
df_long <- melt(df_filtered, id.vars = c("Region", "Tissue"),
                variable.name = "Methylation_Status", value.name = "Count")

# Set facet labels
region_labels <- c("Genes" = "Gene Body", "CpG islands" = "CpG Islands",
                   "Promoters" = "Promoters", "Tiling" = "Tiling Regions")

# Plot
ggplot(df_long, aes(x = Tissue, y = Count, fill = Methylation_Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Region, scales = "free_y", labeller = as_labeller(region_labels)) +
  labs(title = "Differential Methylation by Region and Tissue",
       x = "Tissue Type", y = "Number of Sites") +
  theme_bw() +
  scale_fill_manual(values = c("Hypermethylated" = "red", "Hypomethylated" = "blue"))


#### Add heatmap of genes
setwd("/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Normal_Tissue/analysis/reports/tracks_and_tables_data/csv")
data1<-read.csv("gene_list.csv", header=TRUE, sep=",")
data2<-read.csv("betas_3.csv", header=TRUE, sep=",", check.names=FALSE)
data3<-read.csv("annotation_promoters.csv", header=TRUE, sep=",", row.names=1)
data2$Gene<-data3$symbol[match(data2[,1], rownames(data3))]
data2 <- data2 %>% select(last_col(), everything())
data1$R02C01_Veh_Normal<-data2[['207723860040_R02C01']][match(data1[,1], data2[,1])]
data1$R04C01_Veh_Normal<-data2[['207723860040_R04C01']][match(data1[,1], data2[,1])]
data1$R06C01_Veh_Normal<-data2[['207723860040_R06C01']][match(data1[,1], data2[,1])]
data1$R02C02_HDM_Normal<-data2[['207723860040_R02C02']][match(data1[,1], data2[,1])]
data1$R04C02_HDM_Normal<-data2[['207723860040_R04C02']][match(data1[,1], data2[,1])]
data1$R06C02_HDM_Normal<-data2[['207723860040_R06C02']][match(data1[,1], data2[,1])]

data4<-read.csv("/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Tumor_Tissue/analysis/reports/tracks_and_tables_data/csv/betas_3.csv", header=TRUE, sep=",", check.names=FALSE)
data3<-read.csv("annotation_promoters.csv", header=TRUE, sep=",", row.names=1)
data4$Gene<-data3$symbol[match(data4[,1], rownames(data3))]
data4 <- data4 %>% select(last_col(), everything())
data1$R01C01_Veh_Tumor<-data4[['207723860040_R01C01']][match(data1[,1], data4[,1])]
data1$R03C01_Veh_Tumor<-data4[['207723860040_R03C01']][match(data1[,1], data4[,1])]
data1$R05C01_Veh_Tumor<-data4[['207723860040_R05C01']][match(data1[,1], data4[,1])]
data1$R01C02_HDM_Tumor<-data4[['207723860040_R01C02']][match(data1[,1], data4[,1])]
data1$R03C02_HDM_Tumor<-data4[['207723860040_R03C02']][match(data1[,1], data4[,1])]
data1$R05C02_HDM_Tumor<-data4[['207723860040_R05C02']][match(data1[,1], data4[,1])]



#### Methylation Integration using the EpiMethEx ####
### Tumor Tissue ####
library(dplyr)
setwd("/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Intergration/Tumor_tissue/Promoter")
data1<-read.csv("EpiMethEx_Methylation_Tumor_Tissue_HDM.csv", header=TRUE, sep=",")
data2<-read.csv("annotation_promoters.csv", header=TRUE, sep=",")
data1$Gene<-data2$symbol[match(data1[,1], data2[,1])]
data1 <- data1 %>% select(last_col(), everything())
data1 <- data1[, c(1,7,8,9)]
data3<-read.csv("EpiMethEx_RNA-Seq_Tumor_Tissue_HDM.csv", header=TRUE, sep=",")
common_genes <- inner_join(data1, data3, by = "Gene")
common_genes_unique <- common_genes[!duplicated(common_genes$Gene), ]
Expression<-common_genes_unique[,c(1,5,6,7)]
Methylation<-common_genes_unique[,c(1,2,3,4)]
colnames(Expression)<-c("sample", "M_39_1_T3","M_39_3_T1", "M_40_2_T2")
colnames(Methylation)<-c("sample", "M_39_1_T3","M_39_3_T1", "M_40_2_T2")
write.csv(Methylation, file="Methylation_tumor.csv", row.names=FALSE)
write.csv(Expression, file="Expression_tumor.csv", row.names=FALSE)
Methylation<-read.csv("Methylation_tumor.csv", header=TRUE, sep=",", row.names=1)
Expression<-read.csv("Expression_tumor.csv", header=TRUE, sep=",", row.names=1)
# Ensure column names match
stopifnot(all(colnames(Methylation) == colnames(Expression)))
gene_cor <- data.frame(
  Gene = rownames(Methylation),
  Pearson_r = NA,
  Pearson_p = NA,
  Spearman_r = NA,
  Spearman_p = NA
)

# Loop through each gene
for (i in seq_len(nrow(Methylation))) {
  meth <- as.numeric(Methylation[i, ])
  expr <- as.numeric(Expression[i, ])
  
  if (all(!is.na(meth)) && all(!is.na(expr))) {
    # Pearson
    pear <- cor.test(meth, expr, method = "pearson")
    gene_cor$Pearson_r[i] <- pear$estimate
    gene_cor$Pearson_p[i] <- pear$p.value
    
    # Spearman
    spear <- cor.test(meth, expr, method = "spearman")
    gene_cor$Spearman_r[i] <- spear$estimate
    gene_cor$Spearman_p[i] <- spear$p.value
  }
}

# Remove NA rows
gene_cor <- gene_cor[complete.cases(gene_cor), ]

# Add FDR correction (Benjamini-Hochberg)
gene_cor$Pearson_FDR <- p.adjust(gene_cor$Pearson_p, method = "fdr")
gene_cor$Spearman_FDR <- p.adjust(gene_cor$Spearman_p, method = "fdr")

# View results
head(gene_cor[order(gene_cor$Pearson_FDR), ])
head(gene_cor[order(gene_cor$Spearman_FDR), ])

# Save to file
write.csv(gene_cor, "methylation_expression_correlation_all_stats.csv", row.names = FALSE)






#### Normal Tissue ####
library(dplyr)
setwd("/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Intergration/Normal_tissue/Promoter")
data1<-read.csv("EpiMethEx_Methylation_Normal_Tissue_HDM.csv", header=TRUE, sep=",")
data2<-read.csv("annotation_promoters.csv", header=TRUE, sep=",")
data1$Gene<-data2$symbol[match(data1[,1], data2[,1])]
data1 <- data1 %>% select(last_col(), everything())
data1 <- data1[, c(1,6,7,8)]
data3<-read.csv("EpiMethEx_HDM_Normal_Tissue_RNA_Seq.csv", header=TRUE, sep=",")
common_genes <- inner_join(data1, data3, by = "Gene")
common_genes_unique <- common_genes[!duplicated(common_genes$Gene), ]
Expression<-common_genes_unique[,c(1,5,6,7)]
Methylation<-common_genes_unique[,c(1,2,3,4)]
colnames(Expression)<-c("sample", "M_39_1_H2","M_39_3_H2", "M_40_2_H2")
colnames(Methylation)<-c("sample", "M_39_1_H2","M_39_3_H2", "M_40_2_H2")
write.csv(Methylation, file="Methylation_tumor.csv", row.names=FALSE)
write.csv(Expression, file="Expression_tumor.csv", row.names=FALSE)
Methylation<-read.csv("Methylation_tumor.csv", header=TRUE, sep=",", row.names=1)
Expression<-read.csv("Expression_tumor.csv", header=TRUE, sep=",", row.names=1)
# Ensure column names match
stopifnot(all(colnames(Methylation) == colnames(Expression)))
gene_cor <- data.frame(
  Gene = rownames(Methylation),
  Pearson_r = NA,
  Pearson_p = NA,
  Spearman_r = NA,
  Spearman_p = NA
)

# Loop through each gene
for (i in seq_len(nrow(Methylation))) {
  meth <- as.numeric(Methylation[i, ])
  expr <- as.numeric(Expression[i, ])
  
  if (all(!is.na(meth)) && all(!is.na(expr))) {
    # Pearson
    pear <- cor.test(meth, expr, method = "pearson")
    gene_cor$Pearson_r[i] <- pear$estimate
    gene_cor$Pearson_p[i] <- pear$p.value
    
    # Spearman
    spear <- cor.test(meth, expr, method = "spearman")
    gene_cor$Spearman_r[i] <- spear$estimate
    gene_cor$Spearman_p[i] <- spear$p.value
  }
}

# Remove NA rows
gene_cor <- gene_cor[complete.cases(gene_cor), ]

# Add FDR correction (Benjamini-Hochberg)
gene_cor$Pearson_FDR <- p.adjust(gene_cor$Pearson_p, method = "fdr")
gene_cor$Spearman_FDR <- p.adjust(gene_cor$Spearman_p, method = "fdr")

# View results
head(gene_cor[order(gene_cor$Pearson_FDR), ])
head(gene_cor[order(gene_cor$Spearman_FDR), ])

# Save to file
write.csv(gene_cor, "methylation_expression_correlation_all_stats.csv", row.names = FALSE)



### Scatter Plot Normal

setwd("/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Normal_Tissue/analysis/reports/differential_methylation_data")
data1<-read.csv("Genes_from_RNA-Seq.csv", header=TRUE, sep=",")
data2<-read.csv("Ensembl_ID_genes.csv", header=TRUE, sep=",")
data1$Ensembl_id<-data2$Gene.stable.ID[match(data1[,1], data2[,1])]
data1<-data1[,c(3,2,1)]
data3<-read.csv("diffMethTable_region_cmp1_promoters.csv", header=TRUE,sep=",")
data1$comb.p.val<-data3$comb.p.val[match(data1[,1], data3[,1])]
data1$mean.mean.quot.log2<-data3$mean.mean.quot.log2[match(data1[,1], data3[,1])]
data1<-subset(data1, data1$comb.p.val < 0.05)
data1<-subset(data1, data1$comb.p.val != "NA" & data1$mean.mean.quot.log2 != "NA")


# Define colors
custom_colors <- c(
  "Low_effect_size" = "grey80",
  "Hyper-Up"        = "cyan",  # red
  "Hypo-Up"         = "gold",  # green
  "Hyper-Down"      = "pink",  # blue
  "Hypo-Down"       = "#984EA3"   # purple
)

# Make categories

data1$Category <- "Low_effect_size"  # initialize all rows
data1$Category[data1$log2FoldChange >= 1 & data1$mean.mean.quot.log2 >= 0.1]  <- "Hyper-Up"
data1$Category[data1$log2FoldChange <= -1 & data1$mean.mean.quot.log2 >= 0.1] <- "Hyper-Down"
data1$Category[data1$log2FoldChange >= 1 & data1$mean.mean.quot.log2 <= -0.1] <- "Hypo-Up"
data1$Category[data1$log2FoldChange <= -1 & data1$mean.mean.quot.log2 <= -0.1] <- "Hypo-Down"
write.csv(data1, file="Normal.csv", row.names=FALSE)
normal<-ggplot(data1, aes(x = mean.mean.quot.log2, y = log2FoldChange)) + 
  geom_point(aes(fill = Category), shape = 21, size = 2, color = "black") +
  scale_fill_manual(values = custom_colors) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-1, 1)) +  # <-- zooms in without removing points
  labs(
    x = "log2 Fold Change (Methylation)", 
    y = "log2 Fold Change (RNA-seq)",
    fill = "Category"
  ) +
  theme_classic() + ggtitle("Normal Tissue")


## Scatter Plot Tumor 
setwd("/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Tumor_Tissue/analysis/reports/differential_methylation_data")
data1<-read.csv("genes_from_RNA-Seq.csv", header=TRUE, sep=",")
data2<-read.csv("biomart.csv", header=TRUE, sep=",")
data1$ Gene.stable.ID<-data2$ Gene.stable.ID[match(data1[,1], data2[,1])]
data1<-data1[,c(3,1,2)]
data3<-read.csv("diffMethTable_region_cmp1_promoters.csv", header=TRUE, sep=",")
data1$mean.mean.quot.log2<-data3$mean.mean.quot.log2[match(data1[,1], data3[,1])]
data1$comb.p.val<-data3$comb.p.val[match(data1[,1], data3[,1])]
data1<-subset(data1, data1$comb.p.val < 0.05)
data1<-subset(data1, data1$comb.p.val != "NA" & data1$mean.mean.quot.log2 != "NA")
# Define colors
custom_colors <- c(
  "Low_effect_size" = "grey80",
  "Hyper-Up"        = "cyan",  # red
  "Hypo-Up"         = "gold",  # green
  "Hyper-Down"      = "pink",  # blue
  "Hypo-Down"       = "#984EA3"   # purple
)

# Make categories

data1$Category <- "Low_effect_size"  # initialize all rows
data1$Category[data1$log2FoldChange >= 1 & data1$mean.mean.quot.log2 >= 0.1]  <- "Hyper-Up"
data1$Category[data1$log2FoldChange <= -1 & data1$mean.mean.quot.log2 >= 0.1] <- "Hyper-Down"
data1$Category[data1$log2FoldChange >= 1 & data1$mean.mean.quot.log2 <= -0.1] <- "Hypo-Up"
data1$Category[data1$log2FoldChange <= -1 & data1$mean.mean.quot.log2 <= -0.1] <- "Hypo-Down"
write.csv(data1, file="Tumor.csv", row.names=FALSE)
tumor<-ggplot(data1, aes(x = mean.mean.quot.log2, y = log2FoldChange)) + 
  geom_point(aes(fill = Category), shape = 21, size = 2, color = "black") +
  scale_fill_manual(values = custom_colors) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-1, 1)) +  # <-- zooms in without removing points
  labs(
    x = "log2 Fold Change (Methylation)", 
    y = "log2 Fold Change (RNA-seq)",
    fill = "Category"
  ) +
  theme_classic() + ggtitle("Tumor Tissue")


library(gridExtra)
grid.arrange(normal, tumor, nrow = 2)




### Dotplot
library(ggpubr)

# Read and format data
data <- read.csv("/Users/sharmas30/Documents/WGS3/rna-seq/Methylation/Normal_Tissue/analysis/reports/differential_methylation_data/Dotplot.csv", 
                 header = TRUE, sep = ",")
data <- data[, c(1, 3)]  # Keep only needed columns
rownames(data) <- data$Gene
data$Gene <- NULL 
# Create red-white-blue gradient
my_cols <- colorRampPalette(c("purple", "white", "yellow"))(50)

# Plot
b<-ggballoonplot(data, fill = "value") +
  scale_fill_gradientn(colors = my_cols) + theme_bw()
library(gridExtra)
grid.arrange(a, b, ncol = 2)



#### Immunedeconv ###
library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)


install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")
library("immunedeconv")
install.packages(
  "/Users/sharmas30/Downloads/Matrix_1.6-5.tar",
  repos = NULL,
  type  = "source",
  dependencies = TRUE
)

data<-read.csv("/Users/sharmas30/Documents/WGS3/rna-seq/Immunedeconv_tmp_counts.csv", header=TRUE, sep=",", row.names=1)
res_mMCPcounter <- deconvolute_mouse(data, "DCQ")
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# 1) generate 14 distinct colours
cols14 <- viridis(14)

# 2) capture the original sample‐order from your wide data
samples_order <- colnames(res_mMCPcounter)[-1]

res_mMCPcounter %>%
  gather(sample, score, -cell_type) %>%
  # force sample into a factor with the original levels
  mutate(sample = factor(sample, levels = samples_order)) %>%
  ggplot(aes(x = sample, y = score, color = cell_type)) +
  geom_point(size = 4) +
  facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
  scale_color_manual(values = cols14, guide = FALSE) +
  scale_x_discrete(limits = samples_order) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
data <- read.csv("/Users/sharmas30/Documents/WGS3/rna-seq/mmcp_counter_stats.csv", stringsAsFactors = FALSE)
# 2. Pivot every score column into long format
long <- data %>%
  pivot_longer(
    cols      = -c(sample_id, Type, Exposure),
    names_to  = "cell_type",
    values_to = "score"
  ) %>%
  mutate(
    # create a combined group factor in the order you want
    group = factor(
      paste(Type, Exposure, sep = "_"),
      levels = c("Normal_Veh", "Normal_HDM", "Tumor_Veh", "Tumor_HDM")
    )
  )

# 3. Specify the four cell types to plot (exact names from your data)
selected <- c(
  "B.cell",
  "Macrophage.Monocyte",
  "Endothelial.cell",
  "Cancer.associated.fibroblast"
)

# 4. Filter and plot
ggplot(
  long %>% filter(cell_type %in% selected),
  aes(x = group, y = score)
) +
  geom_boxplot(outlier.shape = NA) +                         # hide boxplot outliers
  geom_jitter(aes(color = cell_type),                        # overlay points
              position = position_jitter(width = 0.2),
              size = 2, alpha = 0.7) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 2) +     # one panel per cell type
  scale_color_manual(values = viridis::viridis(4)) +         # four distinct colors
  labs(
    x     = "Group",
    y     = "mMCPcounter Score",
    title = "Score Distribution by Group (with Individual Points)"
  ) +
  theme_bw() +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


data<-read.csv("/Users/sharmas30/Documents/WGS3/rna-seq/Immunedeconv_tmp_counts.csv", header=TRUE, sep=",",row.names=1)
expr_mat <- as.matrix(data)
mode(expr_mat) <- "numeric"  
res_seqImmuCC <- deconvolute(expr_mat, "seqImmuCC")

res_mMCPcounter <- deconvolute_mouse(data, "seqImmuCC")
dataset_petitprez_humanGenes <- convert_human_mouse_genes(data, convert_to = 'human')
immunedeconv::deconvolute(dataset_petitprez_humanGenes, "seqImmuCC")





