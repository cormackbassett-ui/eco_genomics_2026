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
