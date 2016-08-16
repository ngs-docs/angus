---
title: "Differential Expression analysis using RNA-Seq and DESeq2"
author: "Ian Dworkin"
date: "August 16th, 2016"
output:
  html_document:
    keep_md: yes
---
# NGS2016 Tutorial on Differential Expression analysis using RNA-Seq data and DESeq2

Counts from Salmon are found [here](https://github.com/ngs-docs/angus/blob/2016/_static/quantification.tgz)

## Background
In this tutorial I will provide a basic overview of differential expression analysis for transcriptional profiling using RNA-Seq data. We will be using the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) library in R. This approach utilizes a variant on the assumption of a negative binomially set of counts. This approach assumes that all you have going in are counts, that have not been normalized either for library size (or number of mapped reads), not for transcript length.

Instead of running these analyses on an Amazon EC2 instance, we'll run this locally on our own computers. Before you begin, you will need to download all of the count files we generated using Salmon. Rob and the author of DESeq2 (Mike Love) have developed some nice import tools to get everything into `R` relatively efficiently. However, you will need to use a very recent version of `R` to use these functions.

## Data provenance
This is a small part of the data for a project examining how both variation in genotype and wild type genetic background influence wing morphology in *Drosophila melanogaster*. The sdE3 is a mutant allele in the *scalloped* gene of Drosophila that causes partial loss of wing tissue compared with wild type (wt). The two wild type strains (SAM and ORE) both have normal looking wings as wild type, but differ in the extent of wing loss when the sdE3 allele is introduced into each wild type genetic background.

This image should give you some idea of what it looks like:
![Effects of sdE3 mutation on Drosophila wing morphology in two wild type backgrounds](https://haldanessieve.files.wordpress.com/2013/09/figure1.jpg)

## How counts were generated

See Rob's tutorial on using Salmon. These are directly from that. Please unzip the file in an easy to find place in your directory.


## Get `R` loaded, and let's get started.

First we load in libraries in `R`. In addition to `DESeq2`, there are a few other libraries we will need.

```{r}
library(DESeq2)
library(tximport)
library(readr)
library("RColorBrewer")
library("gplots")
```

It is possible or even likely that you will get an error for some of these, as you have not yet installed the appropriate library. Some are from CRAN (where most R libraries are available), while others are part of bioconductor.

To install for base R (like `gplots`) you can use the:

```{r}
install.packages("gplots")
```

While for `tximport` you will need to do this via bioconductor:

```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("tximport")
```

Depending on the implementation of R you are using (R on mac, R on windows, R studio), there may be some slight differences, so grab an instructor.

Set the working directory for the raw count data (you need to know where you put it). I will go over how I organize my projects to keep this simple.  Just like Unix, R has a current working directory. You can set the working directory using the `setwd()` function. Once you have unzipped the file Rob has provided, you need to navigate to that directory. You will want to be inside the quantification folder. For me this will look like this.

```{r}
#setwd("../data")
#setwd("/Users/ian/NGS_Summer_course/NGS2016/quantification") # This will differ for you!!!

# Setting it up for the import (this is not the import itself)
quant_files <- file.path("quants", list.files("quants"), "quant.sf")

# Let's take a look at this variable we have created
quant_files
```

## Loading the count data into R

DESeq2 and other libraries often have helper functions for getting your count data in. In particular if you are using objects created from other tools that the same authors generated. However, if you are going to make your own pipeline, it is important to know how to write some simple R code to be able to get your data in, and to parse it so that it is in the format you need. I will (if we have time) go through a more typical example where there is no helper functions (so you write it yourself). However, we will use the ones available from `tximport`.

```{r}
# Names of samples.
samples <- c("ORE_sdE3_rep1", "ORE_sdE3_rep2", "ORE_wt_rep1","ORE_wt_rep2", "HYB_sdE3_rep1", "HYB_wt_rep1", "SAM_sdE3_rep1","SAM_sdE3_rep2", "SAM_wt_rep1","SAM_wt_rep2", "HYB_sdE3_rep2", "HYB_wt_rep2")

names(quant_files) <- samples
```

We also need to load in the table of gene and transcript names (with correspondence).  We can first take a look at this file (`txp_to_gene.tsv`) in your terminal. How do you do this?

To import it into `R` we are going to use the most basic and general import tool in base `R`. `read.table()`. There are many options (and I often use read.csv). The data.table library has many useful approaches to this problem (FYI.)


```{r}
tx2gene <- read.table("txp_to_gene.tsv", col.names=c("TXNAME", "GENEID"))

head(tx2gene)
```


Now we can go ahead and read in the input — this will automatically sum results to the gene level.  You can check out the tximport documentation for some other, potentially useful options.

```{r}
txi <- tximport(quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  reader = read_tsv)

# Let's also take a quick look at txi

summary(txi) # why 166884
str(txi)
head(txi$counts) # note these are not integers!


# We can look at patterns of correlations among samples
cor(txi$counts)


pairs(txi$counts[,1:2])
```

## Setting up our experimental design.

DESeq2 needs the information about your experiment set up, so it knows the various predictors in the model (in this case genotype and background). The easiest way to do this is by setting it up as a data frame in R (which is a specialized version of a list). I will (time permitting) show you a more general way of doing this with another example, but for now we are explictly writing this out.

```{r}
background <- c(rep("ORE", 4), rep("HYB", 2), rep("SAM", 4), rep("HYB", 2))

genotype <- c("sdE3", "sdE3", "wt", "wt", "sdE3", "wt", "sdE3", "sdE3", "wt", "wt", "sdE3", "wt")

lane <- c(4,5,1,2,2,6,6,7,3,4,1,5)

lane <- factor(lane) # we will want to treat this as a factor

rna.design <- data.frame(sample=samples,
  file=quant_files,
  background=background,
  genotype=genotype,
  lane = lane)

rna.design

# and we can start with a simple model (back to this later)
load.model <- formula(~ genotype)
```

Below is the crucial function, `DESeqDataSetFromTximport()`, that gives you a DESeq2 count matrix from the txt object. interestingly, this is the step that converts the estimated counts to integers (again we can take a look at the quant_files).

```{r}
all.data <- DESeqDataSetFromTximport(txi, rna.design, design=load.model)
```

## Data is in, now what?

This is normally a good opportunity to do some simple visulizations to look at the distributions of the estimates and the correlations among samples (what should we be looking for). Since Rob already went over this in his tutorial, and you have some sense of the properties of the data, we are going to skip that for the moment.


## Preliminary Quality Control analysis
Before we begin any real analysis. It pays to take some looks at the data. I am not going to go through a full exploratory data analysis session here. But some obvious plots

It is well known that there can be substantial lane to lane variation. For this experiment, it was designed so that 8 samples were run in each lane (barcoded), in a complete block randomized design. This enables us to control for lane effects if necessary. As it turns out I picked a somewhat useless sub-sample of the full data set, so we can not look at the lane effects (as we don't have enough samples in each lane for this subset of data we provide). But normally do something like this (and include a lane effect at a covariate)

First we create a DESeq data object using our counts, experimental design and a simple statistical model (more on this later)

```{r}
load.model <- formula(~ lane)

test_lane_effects <- DESeqDataSetFromTximport(txi,
  rna.design, design=load.model)

test_lane_effects2 <- DESeq(test_lane_effects)
# We now fit the simple model
```

This generates a fairly complex object
```{r, echo=TRUE}
str(test_lane_effects2)
```

For the moment we can ask whether any genes show evidence of different expression based solely on lane to lane variation. We use the results() to summarize some of the results.

```{r}
test_lane_effects2_results <- results(test_lane_effects2, alpha = 0.05)
# alpha = 0.05 is the  "cut-off" for significance (not really - I will discuss).

summary(test_lane_effects2_results)
# 2 genes which may show  evidence of lane effects, but this is a bit incomplete for the full data set.

head(test_lane_effects2_results)

# let's re-order the data to look at the two genes.
test_lane_effects2_results <- test_lane_effects2_results[order(test_lane_effects2_results$padj),]

```

We can also plot the mean-dispersion relationship for this data.

```{r plot Dispersion, echo=TRUE }
plotDispEsts(test_lane_effects2)
```
Let's talk about what this means.

### Principal Components analysis and hierarchical clustering are useful tools to visualize patterns (and to identify potential confounds)

 We can also use some multivariate approaches to look at variation. For PCA (checking it with a "blind" dispersion estimate to look for any funky effects. Not for biological inference).

```{r}
for_pca <- rlog(test_lane_effects2, blind=TRUE)
```
`rlog` is one approach to adjusting for both library size and dispersion among samples. `blind=TRUE`, has it ignore information from the model (in this case lane).

```{r lane effects, echo=TRUE}
plotPCA(for_pca, intgroup=c("lane")) # no obvious lane effects.
```

and now for sex and size

```{r, echo=TRUE}
plotPCA(for_pca, intgroup=c("genotype", "background"))
```

### We can also use some hierarchical clustering

For distance matrix for clustering QC
```{r}
rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion
distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix
```

We need to rename our new matrix of distances based on the samples.
```{r}
rownames(mat) <- colnames(mat) <-   with(colData(test_lane_effects2), paste(genotype, background, sep=" : "))

hc <- hclust(distsRL)  # performs hierarchical clustering
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)  # picking our colours
```

Now we generate the plot
```{r heatmap, echo=TRUE}
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13))
```

 While I checked that there was no evidence of lane effects (see the Likelihood Ratio Test below), I am keeping it in the model as it seems to have little effects. However, it is using up "df" (in the parameter estimation sense), so it may be worth ultimately getting rid of it.



## Proceeding with the real analysis we care about!
Given the results from above, I am removing lane entirely.

```{r}

load.model <- formula(~ genotype)

test_genotype_effects <- DESeqDataSetFromTximport(txi,
  rna.design, design=load.model)

test_genotype_effects2 <- DESeq(test_genotype_effects)

plotDispEsts(test_genotype_effects2)

plotMA(test_genotype_effects2, ylim =c(-1, 1))

genotype_results <- results(test_genotype_effects2, alpha = 0.05)
summary(genotype_results)

# reorder
genotype_results <- genotype_results[order(genotype_results$padj),]
```

While everything is stored, by default DESeq2 is evaluating the final term in the model. In this case evidence of interactions between sex and size. We can look at these in the model we actually want to fit

```{r, echo=TRUE}
load.model <- formula(~ genotype + background + genotype:background) # Let me go over the matrix rank issue here

test_GxB_effects <- DESeqDataSetFromTximport(txi,
  rna.design, design=load.model)

test_GxB_effects2 <- DESeq(test_GxB_effects)

plotDispEsts(test_GxB_effects2)

plotMA(test_GxB_effects2, ylim =c(-2, 2))

GxB_results <- results(test_GxB_effects2, alpha = 0.05, pAdjustMethod="BH")
summary(GxB_results)

# reorder
GxB_results <- GxB_results[order(GxB_results$padj),]
head(GxB_results)

```



Relatively few "significant" genes. Just to keep in mind. We expect a priori, with no true "significant" hits an approximately uniform distribution that looks like (for 12627 genes)

```{r null dist, echo=TRUE}
p_fake <- rbeta(12627, 1,1) # you could also use runif(12627,1,1)
hist(p_fake)
```

But we actually observe
```{r, echo=TRUE}
hist(GxB_results$pvalue)
```

FDR methods exploit this.


### Contrasts with DESeq2

We can now look at specific planned contrasts (in this case male VS female).

```{r}
res_contrast_genotype <- results(test_GxB_effects2,
    contrast=c("genotype", "wt", "sdE3"),
    pAdjustMethod="BH")
summary(res_contrast_genotype, alpha= 0.05)

attr(test_GxB_effects2, "modelMatrixType")   # how is it setting up the design matrix.
#Type of model matrix, expanded or treatment contrast. For more complex models this is very important to understand.

```
DESeq2 can handle a pretty wide array of contrasts. One thing to note is that the notation for some contrasts is quite different from standard `R` contrasts.


## Comparing a complex and reduced model

```{r, echo=TRUE}
full.model <- formula(~ genotype + background + genotype:background) # Let me go over the matrix rank issue here
reduced.model <- formula(~ genotype + background) # Let me go over the matrix rank issue here

test_GxB_effects <- DESeqDataSetFromTximport(txi,
  rna.design, design=load.model)

test_GxB_effects2 <- DESeq(test_GxB_effects, full = full.model, reduced = reduced.model, test= "LRT")

GxB_results2 <- results(test_GxB_effects2 , alpha = 0.05, pAdjustMethod="BH")
summary(GxB_results2)

```
