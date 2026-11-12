## package setup

library(DESeq2) 
library(ggplot2) 
library(dplyr) 
library(pheatmap)

# TASK 1

## 1a
counts <- read.csv("airway_counts.csv", row.names = 1) 
samples <- read.csv("airway_samples.csv")
dim(counts) 
##63677 genes & 8 samples
str(samples) 
colnames(counts) 
samples$run

## 1b
## Updating the library
samples$lib_size <- colSums(counts)
library_summary <- samples %>%
  group_by(group) %>%
  summarize(
    mean_library_size = mean(lib_size),
    number_of_samples = n())

##overall library summary
library_summary

##1c
## plot the bar chart
ggplot(samples, aes(x = run,
                    y = lib_size,
                    fill = group)) + geom_col() +
  labs(
    title = "Library Size per Sample",
    x = "Sample",
    y = "Total Read Counts") +
  theme_minimal()

## save the plot
ggsave("plots/library_sizes.png")


#TASK 2
##2a
## Match sample names with count matrix columns
rownames(samples) <- samples$run

## Set untreated as reference
samples$group <- factor(samples$group, levels = c("untrt", "trt"))

## Creating DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                              colData = samples,
                              design = ~ group)

##Check group levels
levels(dds$group)

##2b
## removing the low count genes with (at least 10 reads in at least 2 samples)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
## Number of genes remaining
nrow(dds)
##17199 genes are remaining

##2c
## Variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

## PCA data
pca_data <- plotPCA(vsd, intgroup = "group")
pca_data
# save the data file
ggsave("plots/pca_data.png", pca_data, width = 7, height = 5, dpi = 300)  

## PC1 variance = 41% and PC2 variance = 26%
## Yes the treated and untreated samples gets separated, along the direction of x-axis.   


#TASK 3
## 3a
dds <- DESeq(dds)
## Normalised counts
write.csv(counts(dds, normalized = TRUE), "results/normalized_counts.csv", quote = FALSE)

## 3b
numerator <- "trt"
denominator <- "untrt"
res <- results(dds, contrast = c("group", numerator, denominator))
## summary of overall DEGs
summary(res)

## create the data frame and add gene id
res <- res[order(res$padj), ]
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

## sorting the up and down genes with respect to Log2fold change 
# centre value genes that are not significant
res_df$status <- "Not significant"
##genes with log2fold greater or equal to 1 
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange >=  1] <- "Up"
##genes with log2fold lesser or equal to 1 
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Down"
##overall table stating the up and downregulated gene
table(res_df$status)
## 388 down regulated and 465 upregulated genes
## total up and downregulated gene
degs <- subset(res_df, status %in% c("Up", "Down"))
nrow(degs)
## 853 genes were up and down-regulated

##3c
##saving the DEG file
write.csv(res_df, paste0("results/", "_all_genes.csv"),  row.names = FALSE)
write.csv(degs,   paste0("results/", "_DEGs.csv"),       row.names = FALSE)
## overall DEG table
deg_table <- read.csv(paste0("results/", "_DEGs.csv"))
head(deg_table, 5)

##3d
# I
## DESeq2 uses raw count data to estimate size factors, dispersion and differential expression. Giving normalized counts would interfere with its statistical modelling assumptions.

## II 
## Some genes have insufficient information for reliable testing, like low counts. DESeq2 assign NA values during filtering or when a gene has outlier issues.

## III
## The gene has an estimated fold of 8 times (2^3) higher expression in treated samples compared with untreated samples.

# TASK 4
## 4a
## remove genes with padj = NA
volcano_data <- res_df[!is.na(res_df$padj), ]    
p_volcano <- ggplot(volcano_data,
                    aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "firebrick",
                                "Down" = "steelblue",
                                "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Comparison of Up and Downregulated gene",
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  theme_bw()
## visualize and save the volcano plot
p_volcano    
ggsave(paste0("plots/volcano", ".png"), p_volcano, width = 6, height = 5, dpi = 300)

## Inference
## In the plot, towards the right (red) are up-regulated, while towards the left are down-regulated (blue). While the middle area represents the not significant genes (statistically). While very few genes has higher p-value indicating the possible outliers present in the data. Also, it is clearly visible that, there are more up-regulated gene (as counted, 465) compared to down-regulated (as counted, 388)

## 4b
## filter the dusp1 gene using Ensembl ID
dusp1 <- res_df %>%
  filter(gene == "ENSG00000120129")
## retrive its log2fold change data and p-value
dusp1[,c("gene","log2FoldChange","padj","status")]

## The gene DUSP1 has an adjusted p-value of 3.63 x 10^45 and log2fold of 2.94, which indicates that it is upregulated by dexamethasone

