# Exploratory RNAseq data anlysis using RMarkdown

## Introduction to RMarkdown 

1. Why RMarkdown?  
	- YAML Header --> renders pretty docs (we will use HTML)   
2. Markdown  
3. Code Chunks  
    - Knitr 
    - Figures 
    - Tables 
4. Inline Code  
5. Breifly  
    - Can do citations & bibliography 
    - Share easily on Rpubs/github  

## Exploratory data analysis with Yeast RNAseq data  

1. Wild Type vs Mutant plots  
	- Facet the plots  
	- Run linear regressions on correlations between the replicates  
2. More Plots  
3. Heatmap


```
# Install all of the necessary packages
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("edgeR")
install.packages("gplots")
install.packages("GGally")
```

I changed the gather_counts.py script to gather the TPM data (transcripts per million) in order to make a plot of this data. This is a simple normalization that can be used to correct for differences in library size. 
```
# load package for edgeR to be able to make a DGE object
library("edgeR")
# Create an object with the tpm file names
files_tpm<-c("ERR458493.fastq.gz.quant.tpm",
             "ERR458494.fastq.gz.quant.tpm",
             "ERR458495.fastq.gz.quant.tpm",
             "ERR458500.fastq.gz.quant.tpm",
             "ERR458501.fastq.gz.quant.tpm",
             "ERR458502.fastq.gz.quant.tpm")

# Create a DGE object from the TPM data. The DGEList contains many data frames that are contained within an object. After reading in the data, we have $samples and we have $counts. We will begin working with the counts data. 
tpm<-readDGE(files_tpm)
```

We will use the package ggplot2 to perform "sniff tests" like Adrienne spoke about. As a first pass, we will make pairwise comparison plots to see how well the replicates actually replicate each other.
Unfortunately, ggplot does not know how to deal with objects of class matrix, so we will convert the tpm data to a dataframe first.
```
tpm_df<-as.data.frame(tpm$counts)
head(tmp_df)
```

We will use a ggplot-esque package to visualize pairwise comparisons between our triplicate samples. 
```
  # Load the package GGally
  library(GGally)
  
  # create a plot using the first three columns of tpm data frame (the WT samples). 
      pairwise_wt <- ggpairs(tpm_df[,1:3], aes(alpha = .05), axisLabels = "internal", title = "WT pairwise comparisons")
      # add a layer to get rid of the grey background that is default in ggplot
      pairwise_wt + theme_minimal()

  # now create a plot using the last three columns of the tpm data frame (the MUT samples)
      pairwise_mut <- ggpairs(tpm_df[,4:6], aes(alpha = .05), axisLabels = "internal", title = "MUT pairwise comparisons")
      # add a layer to get rid of the grey background that is default in ggplot
      pairwise_mut + theme_minimal()
```
These plots show us that our data look good -- it is valid to continue on and do differential expression with these data, because the replicates are true biological replicates

Questions:
1. What species do you work with? Do you expect your biological replicates to be this similar?
2. Would you expect laboratory or environmental samples to be more variable?
3. What would alarming data look like?
4. How is TPM data different than raw counts? How is it different from edgeR or DeSEQ2 corrected results?



Since our data appears to be of high quality, let's perform differential expression with edgeR like we did before. We need to read in the proper count data, because edgeR expects raw counts, not counts that have been corrected in anyway. After we read in our data as a DGEList, we will get:
```
      files <- c(
      "ERR458493.fastq.gz.quant.counts",
      "ERR458494.fastq.gz.quant.counts",
      "ERR458495.fastq.gz.quant.counts",
      "ERR458500.fastq.gz.quant.counts",
      "ERR458501.fastq.gz.quant.counts",
      "ERR458502.fastq.gz.quant.counts"
      )
      
    # Make a character vector named labels that will label the 6 samples above. 
    # Note these are in the same order as the files -- WT are the first three, and MUT are the last three.
    # They will be used later for the MDS plot
    labels=c("WT_1", "WT_2", "WT_3", "MUT_1", "MUT_2", "MUT_3")

    # Read in the data
    data <- readDGE(files)

    # check that the data are correct
    print(data)
```
We then need to describe the group that each sample is a part of, and estimate dispersion
```
    # Create a vector to assign the "group" to the DGE object. This tells edgeR what the treatments here. 
    # We have 6 samples and two treatments (3 WT, and 3 MUT).
    group <- c(rep("wt", 3), rep("mut", 3))
    # create a DGEList with the group variable
    dge <- DGEList(counts=data, group=group)
    # Make dispersion estimates. 
    dge <- estimateCommonDisp(dge)
    dge <- estimateTagwiseDisp(dge)
    # Look at the group variable. 
    head(dge$samples) 
```
Lastly, we perform differential expression
```
    # perform the exact test one the DGE object for differential expression
    et <- exactTest(dge, pair=c("wt", "mut"))
    # grab the top 100,000 tags (n sets the number of tags that you want, the idea for using 100,000 is that all will be returned)
    etp <- topTags(et, n=100000)
    # let's take a look at the differentially expressed genes
    etp # The first thing that is printed show "Comparison of groups:  mut-wt"
```
Differential expression outputs a large file, which may leave you thinking...what next?? We can make an MDS plot to understand how our samples relate to one another. Our WT samples cluster more closely to each other than to the MUT samples, which is what we would expect.
```
  # Make an MDS plot
  plotMDS(dge, labels=labels)
```
We can also make an MA plot
```
    # Make an MA plot
    # because our comparison is mut vs. wt, let's flip the tags for our plot
    etp$table$logFC = -etp$table$logFC
    plot(etp$table$logCPM, etp$table$logFC, xlim=c(-3, 20), ylim=c(-12, 12), pch=20, cex=.3,
        col = ifelse( etp$table$FDR < .2, "red", "black" ) )
```

So now we have our differentially expressed genes. Let's make a volcano plot with the top differentially expressed genes labeled.
```
```

We can also make a histogram of the corrected p value (FDR) can also be informative. 
```
library(ggplot2)
ggplot(etp$table, aes(FDR)) + geom_histogram() + ggtitle("Histogram of FDR for all Differentially Expressed Genes")
```

Heatmaps are also a very popular way to visualize data. Let's make a heat map of the 75 most variable genes!
```
    # modified from http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
    # make a dataframe of the counts that were normalized by edgeR
    norm_counts<-dge$counts

    # We estimate the variance for each row in the logcounts matrix
    var_genes <- apply(norm_counts, 1, var)
    head(var_genes)

    # Get the gene names for the top 75 most variable genes
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:75]
    head(select_var)
    
    # Subset the normalized counts dataframe, grabbing the select_var genes
    highly_variable_norm_counts <- norm_counts[select_var,]
    
    # check the dimensions and the first few rows
    dim(highly_variable_norm_counts)
    head(highly_variable_norm_counts)
 
    # Prepare for the heatmap
          library(RColorBrewer) 
      # Get some nicer colours
          mypalette <- brewer.pal(11,"RdYlBu")
          morecols <- colorRampPalette(mypalette)
    
    # Plot the heatmap
    library(gplots)
    heatmap.2(highly_variable_norm_counts, col=rev(morecols(50)), trace="none", 
              main="Top 75 most variable genes across samples", scale="row",
              labCol = c("WT", "WT", "WT", "MUT", "MUT", "MUT"))
```
Let's do the same thing again, but this time we will plot the top 75 least variable genes
```
    # First, let's get rid of the entries that are not expressed acrossed the board
      keep_counts<-rowSums(cpm(dge)>1) >=2
      dge_expressed <- dge[keep_counts, , keep.lib.sizes=FALSE]
      norm_counts_expressed<-dge_expressed$counts
      
    # We estimate the variance for each row in the logcounts matrix
      var_genes_expressed <- apply(norm_counts_expressed, 1, var)
      
    # Get the gene names for the top 75 most variable genes
      select_var_least<- names(sort(var_genes_expressed, decreasing=FALSE))[1:75]
      head(select_var_least)
    
    # Subset the normalized counts dataframe, grabbing the select_var genes
      not_variable_norm_counts <- norm_counts_expressed[select_var_least,]
    
    # check the dimensions and the first few rows
      dim(not_variable_norm_counts)
      head(not_variable_norm_counts)
   
    # Plot the heatmap
    # install.packages("gplots")
    library(gplots)
    heatmap.2(not_variable_norm_counts, col=rev(morecols(50)), trace="none", 
              main="Top 75 least variable genes across samples", scale="row",
              labCol = c("WT", "WT", "WT", "MUT", "MUT", "MUT"))
```    
Even though we did cpm normalization, a lot of these genes are still very lowly expressed. This is important to keep in mind when interpretting what this heatmap means.
