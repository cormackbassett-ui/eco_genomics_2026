# Transcriptomics Notebook

**Course:** Ecological Genomics 2026

**Name:** Cormac Bassett

------------------------------------------------------------------------

## 9.22.2026 - Continuing the Gene expression analysis tutorial

-   Figuring out how to set up our R working environment and copied the data to import itno DESeq2


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



Principle Component Analysis collapses a matrix of different attributes to data into two separate axes to see how close data clusters are in relation to one another!
# library ("something") will pull up that module for us to use in the R Studio time (like DeSeq2!)
# 
```

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

![This was the work done in the VACC shell today!](file:///users/c/k/ckbasset/projects/eco_genomics_2026/images%20for%20notebook/9.17.2026%20VACC%20shell.png)

**Notes/Observations:**

The VACC shell interface looks really scary

**Next Steps:**
