---
title: "Setting up"
date: "08/15/2015"
output: html_document
layout: topic
---

# Setting up
Script used to generate the bottomly_eset.RData

```{r setup}

# Load data from eset from ReCount
load('bottomly_eset.RData')

```

We want to evaluate our tools in the context of replicates and within group variability. Therefore we need to create datasets with varying number of replicates.

It is usually worthwhile to remove genes that appear not be expressed in any of the experimental conditions (although filtering is done by default with some of the DE tools). We will do that by removing any rows in which three or more samples have less than logCPM of 1. Note that filtering methods involving variances should not be used. The limma algorithm analyses the spread of the genewise variances. 

```{r create-data}

#Filtering out non-expressors
isexpr <- rowSums(cpm(exprs(bottomly.eset))>1) >= 3
sum(isexpr)
bottomly.eset <- bottomly.eset[isexpr,]

# Convert metadata to factors
pData(bottomly.eset)$experiment.number <- factor(pData(bottomly.eset)$experiment.number)
pData(bottomly.eset)$lane.number <- factor(pData(bottomly.eset)$lane.number)

# Get indices for each sample group 
C57 <- which(pData(bottomly.eset)$strain == 'C57BL/6J')
DBA <- which(pData(bottomly.eset)$strain == 'DBA/2J')

# Randomly sample 5 reps from each
reps <- c(sample(C57, 5), sample(DBA, 5))
bottomly.5reps <- bottomly.eset[ ,reps]

# Randomly sample 2 reps from each
reps <- c(sample(C57, 2), sample(DBA, 2))
bottomly.2reps <- bottomly.eset[ ,reps]

# Clean up
rm(reps)
rm(C57)
rm(DBA)
```
