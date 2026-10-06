# Transcriptomics Notebook

**Course:** Ecological Genomics 2026

**Name:** Cormac Bassett

------------------------------------------------------------------------

## 9.16.2026 - Setting up lab notebook and learning markdown

-   Setting up transcriptomics notebook

-   Learn how to take notes in markdown

-   Push notes to github

**Working Directory:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

**Scripts**

`none`

**Code:**

\`\`\` r

print("Hello World")

\`\`\`

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

![](images/markdown-syntax-cheatsheet.webp)

**Notes/Observations:**

# jajaja

**Next Steps:**

# Excited to learn more about how the VACC works!

# Transcriptomics Notebook

**Course:** Ecological Genomics 2026

**Name:** Cormac Bassett

------------------------------------------------------------------------

## 9.17.2026 - Learning basic commands in the VACC shell

-   Figuring out cool commands to muck about in the VACC computing shell

**Working Directory:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

**Scripts**

`none`

**Code:**

\`\`\`

```         
# /gfps/1/cl/biol3990 #how to get into class data set
# ll #long list
# ls short list
# zcat #look into a zip file. DON'T RUN it without 'piping' it to a head command!
# head -n #number of to display after a zcat
# -wc #number of lines in file
# history #shows all recent commands
```

\`\`\`

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

![This was the work done in the VACC shell today!](file:///users/c/k/ckbasset/projects/eco_genomics_2026/images%20for%20notebook/9.17.2026%20VACC%20shell.png)

**Notes/Observations:** The VACC shell interface looks really scary

**Next Steps:**

## 9.22.2026 - Continuing the Gene expression analysis tutorial

-   Figuring out how to set up our R working environment and copied the data to import into DESeq2

**Working Directory:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics/`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1/tidyverse`

-   `R-Studio`

**Scripts**

`none`

**Code:**

\`\`\`

```         
# mkdir("something") lets us make a new directory at whatever level we are in
# countsTable will format a data packet for us into a chart
# countsTable round will round the figures to the closest integer
# hist(apply(countsTableRound,1,mean),xlim=c(0,1000), ylim=c(0,10000),breaks=10000) turns the data into a readable histogram to understand
# dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ generation + treatment)
# The above command let us create a DeSeq object and labe lit via the tilda
# dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
# nrow(dds) 
# Above command filtered out the genes for only those btween 15 and 28 strand longs!
# dds <- DESeq(dds) This command runs a differential gene exrpession anaysis on the function
# resultsNames(dds) This command spits out the names of what groups we just made
# vsd <- vst(dds, blind=FALSE)
# meanSdPlot(assay(vsd))
# sampleDists <- dist(t(assay(vsd))) This whole section of commands bascially just tries to vet the data for variance and stabilizes it
# 
library("RColorBrewer")
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$line, vsd$generation, sep="-")
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
# This WHOOOLE section just sets up the histogram that allows us to analyze for outliers.
# 
sampleTree <- hclust(dist(sampleDists), method="average")
# plot
# plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="",cex.lab=1.5, cex.axis=1.5, cex.main=2)
# The above section lets us actually check for outliers!
# 
```

\`\`\`

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

**Notes/Observations:** \# This is really complicated so far but I enjoy this!

**Next Steps:** \# Work more in DESeq2 in the future!

-   

    # Transcriptomics Notebook

**Course:** Ecological Genomics 2026

**Name:** Cormac Bassett

## 9.24.2026 - Setting up lab notebook and learning markdown

-   Figuring out how to set up our R working environment, as well as noting common bash commands

------------------------------------------------------------------------

**Working Directory:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1/tidyverse`

-   `R-Studio`

**Scripts**

`none`

**Code:**

\`\`\` r

print("Hello World")

\`\`\`

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

**Notes/Observations:**

**Next Steps:**

------------------------------------------------------------------------

## 9.29.2026 - Day 4 of Differential gene expression analysis

-   Analyzing the counts matrix using a simplified data set
-   Understand what a contrast is and up vs. down regulated
-   Learning how to make common types fo differential gene expression evaluations **Working Directory:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics/mydata`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

**Scripts**

```         
# BASH CODE
pwd - print working directory - tells you what directory you are on
# zcat - prints a file (don't use in large files, as you'll open the whole thing!)
#head - just scans the top of a file
# ll and ls - lists possible directories to move into
# cd - change directory
# history - shows all past commands
# cp - copy
# rm - remove
# .. - moves back a directory
# . - from the directory I am in
# /gpfs1
  |
  |->/cl
  |
  |-/biol3390
  |     |
  |     |->/Transcriptomics
  |            |
  |            |->CountsMatrix
  |
  |->/ecogen
      |
      |->/sw
        |->setup.sh
#
#
```

**Code:**

```         
# # RStudio Coding
# x -> 5 lets you define an 'object'or variable
# data.frame( ) lets you group multiple data pieces together into a graph
# Using "####" at the end of a # piece of text creates a section that you can jump back to later!
```

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

**Notes/Observations:**

**Next Steps:**

## 10.01.2026 - Day 5 of Differential gene expression analysis - wrapup, and GO Expression Tutorial

-   Going over how to create a scatter plot to compare relative differences between populations, understand how GO expression and gene enrichment works, and understand how WGCNA works

<!-- -->

-   `/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics/mydata`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1/tidyverse`

-   `R-Studio`

**Scripts**

`~/projects/eco_genomics_2026/transcriptomics/myscripts`

**Code:**

```         
# filter() to remove rows
# mutate() to add a new variable
# case_when() to classify genes into categories
# arrange() to sort the rows
# 
# 
```

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

**Notes/Observations:**

**Next Steps:**

# Transcriptomics Notebook

**Course:** Ecological Genomics 2026

**Name:** Cormac Bassett

------------------------------------------------------------------------

## 10.06.2026 - Setting up lab notebook and learning markdown

-   Working on WGCNA analysis and also GO analysis

**Working Directory:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/c/k/ckbasset/Projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1/tidyverse`

-   `R-Studio`

**Scripts**

```         
 res_OWAvsAM.df <- as.data.frame(res_OWAvsAM)
res_OWAvsAM.df$fullID <- rownames(res_OWAvsAM.df)

parts <- strsplit(res_OWAvsAM.df$fullID, "::")

res_OWAvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM.df,
  "myresults/F0_OWAvsAM_results.csv",
  row.names = FALSE
)
# This is an example piece of code that was used to save the created results 
# file with the abbreviated Trinity ID's - This was used AFTER importing the 
# counts matrix and setting up + processing with DESeq2

library(topGO)

mappingFile <- "mydata/trinotate_annotation_GOblastx_forTopGO.txt"

geneID2GO <- readMappings(
  file = mappingFile
)

cat("Genes in GO mapping:",
    length(geneID2GO),
    "\n")
# This chunk of code prepares the data folders we pulled for Top GO analysis
# by pulling out the relative gene ID's and lining them up to be run

allGO <- table(unlist(geneID2GO))

keepTerms <- names(allGO)[
  allGO >= 5 &
  allGO <= 500
]

geneID2GO.filtered <- lapply(
  geneID2GO,
  function(x) intersect(x, keepTerms)
)

geneID2GO.filtered <- geneID2GO.filtered[
  lengths(geneID2GO.filtered) > 0
]
# This chunk filters our GO terms by size to kick out the little ones that 
# would muck up the data and result in no conclusive answers. Without doing
# this, the GO results would be just saying that our protein does 'biology'!


run_topGO_contrast <- function(
    infile,
    outfile,
    ontology = "BP",
    padj.cutoff = 0.05){

  deseq <- read.csv(
    infile,
    stringsAsFactors = FALSE
  )

  deseq <- subset(
    deseq,
    !is.na(padj)
  )

  geneList <- factor(
    as.integer(deseq$padj < padj.cutoff)
  )

  names(geneList) <- deseq$shortID

  cat("\nGenes tested:",
      length(geneList))

  cat("\nSignificant genes:",
      sum(geneList == 1),
      "\n")

  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )

  resultWeight <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )

  GOresults <- GenTable(
    GOdata,
    weightFisher = resultWeight,
    orderBy = "weightFisher",
    topNodes = 100
  )

  write.csv(
    GOresults,
    outfile,
    row.names = FALSE
  )

  return(GOresults)
}
# THIS WHOLE chunk sets up the GO function by defining it, setting our
# variables, and setting up the DESeq2 connection.

GO_OA <- run_topGO_contrast(
  "myresults/F0_OAvsAM_results.csv",
  "myresults/GO_F0_OAvsAM_BP.csv"
)

GO_OW <- run_topGO_contrast(
  "myresults/F0_OWvsAM_results.csv",
  "myresults/GO_F0_OWvsAM_BP.csv"
)

GO_OWA <- run_topGO_contrast(
  "myresults/F0_OWAvsAM_results.csv",
  "myresults/GO_F0_OWAvsAM_BP.csv"
)
# These three above set the 3 separate contrasts for the GO function to munch on
# Read TopGO results
go <- read.csv(
  "myresults/GO_F0_OWAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: OWA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Significant\nGenes"
  )
#This whole cluster of code sets up a bubble plot for the GO function to 
# express itself, like converting p-values adn counts to numerics, filter
# to the top 10 GO terms, and plot it in GGPlot!
# REMEMBER: The GGplot significant gene numberdoesn't compare to the whole, so 
# this may be misleading! Weighting it as a proportion helps to clear this up

## Set your working directory
setwd("~/YOURWORKINGDIRCTORY")

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

countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
# conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
# head(conds)

sample_metadata = read.table(file = "Ahud_trait_data.txt",header=T, row.names = 1)

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

soft_power <- 6 # should try running with 7 based on R^2
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package
```

**Code:**

```         

print("Hello World")
```

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

Here are the three GO plots created today!

![](~/projects/eco_genomics_2026/transcriptomics/mydata/myresults/Top%20GO%20Terms:%20OA%20vs.%20AM.png)

![](~/projects/eco_genomics_2026/transcriptomics/mydata/myresults/Top%20GO%20Terms:%20OWA%20vs.%20AM.png)

![](~/projects/eco_genomics_2026/transcriptomics/mydata/myresults/Top%20GO%20Terms:%20OW%20vs.%20AM.png)

**Notes/Observations:**

**Next Steps:**

# Excited to learn more about how the VACC works!
