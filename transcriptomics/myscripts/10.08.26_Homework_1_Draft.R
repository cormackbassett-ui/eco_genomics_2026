# Basic plan for the code is to reuse the WGCNA code provided to us in day 6, 
# fiddle around with the soft threshold powers to create 3 separate sharpness
#  contrasts in the data, pull out the differences in the number of modules, 
# genes per module, and strength of correlation with traits, AND Go functional 
# module for targeted traits - traits investigated by the 2023 DeMayo paper.


setwd("~/projects/eco_genomics_2026/transcriptomics")

# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

library(DESeq2)
library(ggplot2)

library(tidyverse)

library(CorLevelPlot) 
library(gridExtra)

library(Rmisc) 


# 1. Import the counts matrix and metadata and filter using DESeq2

countsTable <- read.table("mydata/salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
# conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
# head(conds)

sample_metadata = read.table(file = "mydata/Ahud_trait_data.txt",header=T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=sample_metadata, 
                              design= ~ 1)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)


# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree) 



# pca - method 2

pca <- prcomp(t(countsTable))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# All of our samples look good; the genes that don't look good shouldn't cluster well and 
# should end up in the "grey", catch all module for uncorrelated genes.

# 3. Normalization ----------------------------------------------------------------------

colData <- row.names(sample_metadata)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(countsTableRound)) # to see if all samples are present in both
all(rownames(colData) == colnames(countsTableRound))  # to see if all samples are in the same order



# perform variance stabilization
dds_norm <- vst(dds)
# dds_norm <- vst(normalized_counts)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function; this step takes a couple minutes
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices 

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)
# based on this plot, choose a soft power to maximize R^2 (above 0.8) and minimize connectivity
# for these ahud data: 6-8; Higher R2 should yield more modules.


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

#Labeling soft powers
soft_power <- 4
soft_power <- 6
soft_power <- 8

temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package

# Labelling the .rds files
bwnet_4 <- readRDS("bwnet.4.rds")
bwnet_6 <- readRDS("bwnet.6.rds")
bwnet_8 <- readRDS("bwnet.8.rds")

cor <- temp_cor

###Bwent with threshold of 4 
module_eigengenes_4 <- bwnet_4$MEs

head(module_eigengenes_4)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_4$dendrograms[[1]], cbind(bwnet_4$unmergedColors, bwnet_4$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

##All Together
nets <- list("Power 4" = bwnet_4,
             "Power 6" = bwnet_6,
             "Power 8" = bwnet_8)

# Function: dendrogram + unmerged/merged module colors
plot_dendro <- function(net, title, block = 1) {
  # genes belonging to this block (handles multi-block runs safely)
  idx <- net$blockGenes[[block]]
  
  plotDendroAndColors(
    net$dendrograms[[block]],
    cbind(net$unmergedColors[idx], net$colors[idx]),
    groupLabels = c("Unmerged", "Merged"),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste("Gene Dendrogram and Module Colors -", title)
  )
}

###Bwent with threshold of 6 ---------------------------------------------------------
module_eigengenes_6 <- bwnet_6$MEs

head(module_eigengenes_6)

# get number of genes for each module
table(bwnet_6$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_6$dendrograms[[1]], cbind(bwnet_6$unmergedColors, bwnet_6$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

###Bwent with threshold of 8 ---------------------------------------------------------
module_eigengenes_8 <- bwnet_8$MEs

head(module_eigengenes_8)

# get number of genes for each module
table(bwnet_8$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_8$dendrograms[[1]], cbind(bwnet_8$unmergedColors, bwnet_8$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

#Stitching them together
for (nm in names(nets)) {
  png(paste0("dendrogram_", gsub(" ", "_", nm), ".png"),
      width = 2400, height = 1200, res = 200)
  plot_dendro(nets[[nm]], nm)
  dev.off()
}

###Bwent with threshold of 4 ---------------------------------------------------------
module_eigengenes_4 <- bwnet_4$MEs

head(module_eigengenes_4)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_4$dendrograms[[1]], cbind(bwnet_4$unmergedColors, bwnet_4$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

##All Together
nets <- list("Power 4" = bwnet_4,
             "Power 6" = bwnet_6,
             "Power 8" = bwnet_8)

# Function: dendrogram + unmerged/merged module colors
plot_dendro <- function(net, title, block = 1) {
  # genes belonging to this block (handles multi-block runs safely)
  idx <- net$blockGenes[[block]]
  
  plotDendroAndColors(
    net$dendrograms[[block]],
    cbind(net$unmergedColors[idx], net$colors[idx]),
    groupLabels = c("Unmerged", "Merged"),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste("Gene Dendrogram and Module Colors -", title)
  )
}

###Bwent with threshold of 6 ---------------------------------------------------------
module_eigengenes_6 <- bwnet_6$MEs

head(module_eigengenes_6)

# get number of genes for each module
table(bwnet_6$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_6$dendrograms[[1]], cbind(bwnet_6$unmergedColors, bwnet_6$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

###Bwent with threshold of 8 ---------------------------------------------------------
module_eigengenes_8 <- bwnet_8$MEs

head(module_eigengenes_8)

# get number of genes for each module
table(bwnet_8$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_8$dendrograms[[1]], cbind(bwnet_8$unmergedColors, bwnet_8$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

#Stitching them together
for (nm in names(nets)) {
  png(paste0("dendrogram_", gsub(" ", "_", nm), ".png"),
      width = 2400, height = 1200, res = 200)
  plot_dendro(nets[[nm]], nm)
  dev.off()
}


install.packages("magick")
library(magick)
setwd("~/projects/eco_genomics_2026/transcriptomics/myimages")

files <- c("Cluster_Dendrogram_4.png",
           "Cluster_Dendrogram_6.png",
           "Cluster_Dendrogram_8.png")

imgs <- image_read(files)

# Stack vertically (one above the other), which suits wide dendrograms
stacked <- image_append(imgs, stack = TRUE)
image_write(stacked, "dendrograms_stacked.png", format = "png")

# View it
print(stacked)


setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")

###Stats for different thresholds
#Number of modules
length(unique(bwnet_4$colors)) 
length(unique(bwnet_6$colors)) 
length(unique(bwnet_8$colors))

#Percentage of genes changing module assignment
mean(bwnet_4$colors != bwnet_6$colors) 
mean(bwnet_4$colors != bwnet_8$colors) 
mean(bwnet_6$colors != bwnet_8$colors)

#Module eigengene similarity
cor(bwnet_4$MEs, bwnet_6$MEs)
cor(bwnet_4$MEs, bwnet_8$MEs)

#Graph comparing eigengenes

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_8$dendrograms[[1]], cbind(bwnet_8$unmergedColors, bwnet_8$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

##All Together
nets <- list("Power 4" = bwnet_4,
             "Power 6" = bwnet_6,
             "Power 8" = bwnet_8)

# 1. Per-sample eigengene values
p_me <- ggplot(me_long, aes(x = Sample, y = Eigengene, fill = Power)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  facet_wrap(~ Module, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(title = "Module eigengenes from WGCNA Soft Thresholds",
       x = "Sample", y = "Eigengenes", fill = "Soft power") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

p_me

ggsave("~/projects/eco_genomics_2026/transcriptomics/myimages/eigengene_values_by_power.png",
       p_me, width = 14, height = 3 * length(shared), dpi = 200, limitsize = FALSE)

# 2. How similar is each module's eigengene to power 4
cor_df <- rbind(
  data.frame(Module = shared,
             Comparison = "Power 4 vs 6",
             Correlation = sapply(shared, function(m) WGCNA::cor(ME6[[m]], ME4[[m]]))),
  data.frame(Module = shared,
             Comparison = "Power 4 vs 8",
             Correlation = sapply(shared, function(m) WGCNA::cor(ME8[[m]], ME4[[m]])))
)
cor_df$Module <- sub("^ME", "", cor_df$Module)

p_cor <- ggplot(cor_df, aes(x = Module, y = abs(Correlation), fill = Comparison)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  labs(title = "Eigengene similarity to power 4",
       x = "Modules", y = "Correlation of eigengenes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_cor

ggsave("~/projects/eco_genomics_2026/transcriptomics/myimages/eigengene_similarity_barplot.png",
       p_cor, width = 10, height = 5, dpi = 250)


# # of Genes across thresholds

get_module_sizes <- function(colors, threshold) {
  as.data.frame(table(colors)) |>
    dplyr::rename(Module = colors, Count = Freq) |>
    dplyr::mutate(Threshold = threshold)
}

df <- rbind(
  get_module_sizes(bwnet_4$colors, "4"),
  get_module_sizes(bwnet_6$colors, "6"),
  get_module_sizes(bwnet_8$colors, "8")
)

# Order modules by size at power 4 so the x-axis is easier to read
mod_order <- df$Module[df$Threshold == "4"][order(-df$Count[df$Threshold == "4"])]
df$Module <- factor(df$Module, levels = unique(c(as.character(mod_order),
                                                 as.character(df$Module))))

p <- ggplot(df, aes(x = Module, y = Count, fill = Threshold)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.75) +
  theme_bw() +
  labs(
    title = "Module Size Across WGCNA Soft Thresholds",
    x = "Module Color",
    y = "Gene Count",
    fill = "Soft power"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

ggsave("~/projects/eco_genomics_2026/transcriptomics/myimages/module_sizes_barplot.png",
       p, width = 12, height = 6, dpi = 250)


###Visualizing module correlation
library(pheatmap)

#4 v 6
corME_46 <- cor(bwnet_4$MEs,
                bwnet_6$MEs)

pheatmap(corME_46,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE)

#4 v 8
corME_48 <- cor(bwnet_4$MEs,
                bwnet_8$MEs)

pheatmap(corME_48,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE)

##All together
ME4 <- bwnet_4$MEs
ME6 <- bwnet_6$MEs
ME8 <- bwnet_8$MEs

# Add prefixes so modules are identifiable
colnames(ME4) <- paste0("P4_", colnames(ME4))
colnames(ME6) <- paste0("P6_", colnames(ME6))
colnames(ME8) <- paste0("P8_", colnames(ME8))

# Combine eigengenes
allMEs <- cbind(ME4, ME6, ME8)

# Correlation matrix
MEcor <- cor(allMEs)


pheatmap(
  MEcor,
  color = colorRampPalette(c("#E5F5F9", "#2CA25F", "#1A1A1A"))(100),
  border_color = NA,
  main = "Correlation Across Soft Threshold Powers"
)

#Clustering it
pheatmap(
  MEcor,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  color = colorRampPalette(c("#E5F5F9", "#2CA25F", "#1A1A1A"))(100)
)

#Gene Reassignment
table46 <- table(
  Power4 = bwnet_4$colors,
  Power6 = bwnet_6$colors
)

table48 <- table(bwnet_4$colors,
                 bwnet_8$colors)

table68 <- table(bwnet_6$colors,
                 bwnet_8$colors)

pheatmap(log10(table46 + 5))
pheatmap(log10(table48 + 5))
pheatmap(log10(table68 + 5))


ggsave("Correlation Across Soft Threshold Power.png",
       p_cor, width = 12, height = 6, dpi = 250)
