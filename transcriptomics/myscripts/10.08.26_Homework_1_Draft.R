# Basic plan for the code is to reuse the WGCNA code provided to us in day 6, 
# fiddle around with the soft threshold powers to create 3 separate sharpness
#  contrasts in the data, pull out the differences in the number of modules, 
# genes per module, and strength of correlation with traits
# - traits investigated by the 2023 DeMayo paper.


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
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

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

module_eigengenes <- bwnet$MEs

head(module_eigengenes)


# getting number of genes per module
table(bwnet_4$colors)
table(bwnet_6$colors)
table(bwnet_6$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")

#Generating the comparison to power 4 from 6 and 8

library(WGCNA)
library(ggplot2)

# Make sure to debug this chunk of code, I think it's configured to set row count
# as 0, which is annoying, but it still works, so interesting! 
cor_df$Module <- sub("^ME", "", cor_df$Module)

# Plotting the two separate power comparisons
correlation_df <- rbind(
  data.frame(
    Module     = shared_valid,
    Comparison = "Power 6 vs 4", "Power 6 vs 4",
    Correlation = sapply(shared_valid, function(m) {
      cor(ME4[[m]], ME6[[m]])
    })
  ),
  data.frame(
    Module     = shared_valid,
    Comparison = "Power 8 vs 4","Power 8 vs 4",
    Correlation = sapply(shared_valid, function(m) {
      cor(ME4[[m]], ME8[[m]])
    })
  )
)

# Cleaning up module genes


cor_df$Module <- sub("^ME", "", cor_df$Module)

p_cor <- ggplot(cor_df, aes(x = Module, y = abs(Correlation), fill = Comparison)) +
  geom_col(position = position_dodge(width = 1.0), width = 0.75) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  labs(
    title = "Eigengene Module Similarity",
    x = "Modules",
    y = "Correlation"
  ) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

p_cor


ggsave("~/projects/eco_genomics_2026/transcriptomics/myimages/Eigengene_Module_Similarity.png",
       p_cor, width = 10, height = 5, dpi = 250)




# Creating the module size across WGCNA thresholds graph - make sure to debug
df <- rbind(
  as.data.frame(table(bwnet_colors)) %>%
    
    rename(Module = Var1, Count = Freq)%>%
    mutate(Threshold = "4"),
  
  as.data.frame(table(bwnet_6$colors)) %>%
    rename(Module = bwnet_6.colors, Count = Freq) %>%
    mutate(Threshold = "6"),
  
  as.data.frame(table(bwnet_8$colors)) %>%
    rename(Module = bwnet_8.colors, Count = Freq) %>%
    mutate(Threshold = "8")
)

#Actually buildign module size graph
mod_order <- df$Module[df$Threshold == "4"][order(-df$Count[df$Threshold == "4"])]
df$Module <- factor(df$Module, levels = mod_order)

p <- ggplot(df, aes(x = Module, y = Count, fill = Threshold)) +
  geom_col(width = 0.75) +
  facet_wrap(~ Threshold, nrow = 1) +
  theme_bw() +
  labs(
    title = "Different Module Size Across WGCNA Soft Thresholds",
    x = "Module Color",
    y = "Gene Count"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

ggsave("~/projects/eco_genomics_2026/transcriptomics/myimages/Different_Module_Size_Across_WGCNA_Soft_Thresholds.png",
       p_cor, width = 10, height = 5, dpi = 250)



# SETTING UP HEATMAP FOR 6 YAYAAY (REMEMBER TO CHANGE THE DIRECTORY OR THIS WONT WORK AND UR COMPUTER EXPLODES)

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

traits <- sample_metadata[, c(5,8,11,14,17)]


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data)

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[12:16],
             y = names(heatmap.data)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"))

# SETTING UP 8 SOFT THRESHOLD HEATMAP YIPPEE HOORAYYYY


soft_power <- 8 # should try running with 7 based on R^2
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package



saveRDS(bwnet, file = "bwnet.8.rds")

bwnet <- readRDS("bwnet.8.rds")

cor <- temp_cor
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                   c("unmerged", "merged"),
                   dendroLabels = FALSE,
                   addGuide = TRUE,
                   hang= 0.03,
                   guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module


# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

traits <- sample_metadata[, c(5,8,11,14,17)]


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)


# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data)

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[12:16],
             y = names(heatmap.data)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"))

# CODE FOR 4 HEATMAP WAHOOO


bwnet <- readRDS("bwnet.4.rds")

cor <- temp_cor

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module


# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

traits <- sample_metadata[, c(5,8,11,14,17)]


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data)

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[12:16],
             y = names(heatmap.data)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"))
