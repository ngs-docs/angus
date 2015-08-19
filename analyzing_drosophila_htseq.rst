===================================
Analyzing RNA-seq counts with DESeq
===================================

In this tutorial, we will continue analyzing the Drosophila RNA-seq data from earlier to look for differentially expressed genes.

Instead of running these analyses on an Amazon EC2 instance, we'll run this locally on our own computers. Before you begin, you will need to download all of the count files we generated using HTSeq. (You can use scp at the command line, or WinSCP to download them.) Place them all in a single folder on your computer. Then, start R.

First, we need to make sure that R can find the files when we try to load them. Just like Unix, R has a current working directory. You can set the working directory using the setwd() function::

    setwd("...")

Be sure to replace the "..." with the path to the folder where you have placed your files.

Right now, we have our data spread out across multiple files, but we need to combine them all into a single dataset for DESeq to be able to use them. There are some convenience functions (see below) to perform this. But for many programs you need to write your own, so we will go through getting it all in writing the function ourselves.

 Let's start by creating a vector, to store the names of the samples and then write a function that will read in a single dataset given a sample name::

    samples <- c("ORE_wt_rep1","ORE_wt_rep2","ORE_sdE3_rep1","ORE_sdE3_rep2","SAM_wt_rep1","SAM_wt_rep2","SAM_sdE3_rep1","SAM_sdE3_rep2","HYB_wt_rep1","HYB_wt_rep2","HYB_sdE3_rep1","HYB_sdE3_rep2")

    #A function to read one of the count files produced by HTSeq
    read.sample <- function(sample.name) {
	    file.name <- paste(sample.name, "_htseq_counts.txt", sep="")
	    result <- read.delim(file.name, col.names=c("gene", "count"), sep="\t", colClasses=c("character", "numeric"))
    }

Now let's try it out::

    #Read the first sample
    sample.1 <- read.sample(samples[1])

Let's double check that we've loaded the first sample properly by looking at the first few rows and seeing how many rows there are::

    head(sample.1)
    nrow(sample.1)
    
Do you get the output you expected? Remember you can always check by using :code:`wc -l yourfile` at the shell

Now let's read the second sample, and double check that it has the same number of rows as the first sample. Before we merge the two datasets  by :code:`cbind()`ing the count columns together, we should make sure that the same gene is represented in each row in the two datasets::

    #Read the second sample
    sample.2 <- read.sample(samples[2])

    #Let's make sure the first and second samples have the same number of rows and the same genes in each row
    nrow(sample.1) == nrow(sample.2)
    all(sample.1$gene == sample.2$gene)

Looks good. Now let's do all the merging using a simple for loop::

    #Now let's combine them all into one dataset
    all.data <- sample.1
    all.data <- cbind(sample.1, sample.2$count)
    for (c in 3:length(samples)) {
	    temp.data <- read.sample(samples[c])
	    all.data <- cbind(all.data, temp.data$count)
    }

    #We now have a data frame with all the data in it:
    head(all.data)

You'll notice that the column names are not very informative. We can replace the column names manually to something more useful like this::

    colnames(all.data)[2:ncol(all.data)] <- samples

    #Now look:
    head(all.data)

    #Let's look at the bottom of the data table
    tail(all.data)


When you look at the bottom of the dataset, you'll notice that there are some rows we don't want to include in our analysis. We can remove them easily by taking a subset of the data that includes everything except the last 5 rows::

    all.data <- all.data[1:(nrow(all.data)-5),]

    tail(all.data)

Now we're ready to start working with DESeq. If you don't already have it installed on your computer, you will want to install it like this::

    source("http://bioconductor.org/biocLite.R")
    biocLite("DESeq")

Once it's installed, we need to load the library like this::

    library("DESeq")


Before we're quite ready to work with the data in DESeq, we need to re-format it a little bit more. DESeq wants every column in the data frame to be counts, but we have a gene name column, so we need to remove it. We can still keep the gene names, though, as the row names (just like each column has a name in a data frame in R, each row also has a name).
::

    #Remove the first column
    raw.deseq.data <- all.data[,2:ncol(all.data)]
    #Set row names to the gene names
    rownames(raw.deseq.data) <- all.data$gene

    head(raw.deseq.data)

Now we have our data, but we need to tell DESeq what our experimental design was. In other words, say we want to look for genes that are differentially expressed between the two genetic backgrounds, or between the mutant and wild-type genotypes--DESeq needs to know which samples belong to each treatment group. We do this by creating a second data table that has all the sample information. We could do this by creating another data table in a text editor (or if you absolutely must, something like Excel, but be careful because sometimes R has trouble reading files generated by Excel, even if you've saved them as "flat" tab-delimited or comma-delimited files).

But instead, we'll generate this sample information table manually in R since it's not very complicated in this case::

    #Create metadata
    wing.design <- data.frame(
	    row.names=samples,
	    background=c(rep("ORE", 4), rep("SAM", 4), rep("HYB", 4)),
	    genotype=rep( c("wt", "wt", "sdE3", "sdE3"), 3 ),
	    libType=rep("paired-end", 12)
    )
    #Double check it...
    wing.design 

By default R picks the "reference" level for each treatment alphanumerically. So in this case the reference level for background would be "HYB" and for genotype it would be "sdE3". However it will be a bit easier for us to intepret the data if we used the wild type (genotype="wt") as a reference. Also for the background, using one of the pure strains, and not their F1 hybrid may help. We accomplish this as follows::

     wing.design$genotype <- relevel(wing.design$genotype, ref="wt")
     wing.design$background <- relevel(wing.design$background, ref="ORE")
     
Now we can create a DESeq data object from our raw count table and our experimental design table::

    deseq.data <- newCountDataSet(raw.deseq.data, wing.design)

Note that DESeq also has a function newCountDataSetFromHTSeqCount that can automatically handle merging all the raw data files together (from HTSeq), but it's useful to be able to do this manually because not every statistical package you work with will have this functionality.

Now we have to do the "normalization" step and estimate the dispersion for each gene::

    deseq.data <- estimateSizeFactors(deseq.data)
    deseq.data <- estimateDispersions(deseq.data)

We have used the "vanilla settings" for both, but it highly recommedable to look at the defaults for the size factors and dispersions. As we discussed, how the gene-wise estimates for dispersions are estimated can have pretty substantial effects.

Let's make sure the dispersion estimates look reasonable::

    plotDispEsts(deseq.data)

In this case it looks ok, although the number of points on the graph is relatively modest compared to most RNA-seq studies, since we have intentionally included only genes on the X chromosome. Note that if the dispersion estimates don't look good you may need to tweak the parameters for the estimateDispersions() function (e.g., maybe try using fitType="local").

Now let's fit some models. These steps are a little bit slower, though not terribly so::

    fit.full <- fitNbinomGLMs(deseq.data, count ~ background + genotype + background:genotype)
    fit.nointeraction <- fitNbinomGLMs(deseq.data, count ~ background + genotype)
    fit.background <- fitNbinomGLMs(deseq.data, count ~ background)
    fit.genotype <- fitNbinomGLMs(deseq.data, count ~ genotype)
    fit.null <- fitNbinomGLMs(deseq.data, count ~ 1)

Now that we've fit a bunch of models, we can do pairwise comparisons between them to see which one best explains the data [for each gene]. For example, we can ask, for which genes does an interaction term between wt/mutant genotype and genetic background help explain variation in expression?::

    #Generate raw p-values for the first comparison: full model vs. reduced model without an interaction term; significant p-values will tell you that the full/more complex model does a "significantly" better job at explaining variation in gene expression than the reduced/less complex model (in this case, the one without the interaction term)
    pvals.interaction <- nbinomGLMTest(fit.full, fit.nointeraction)
    #Generate p-values adjusted for multiple comparisons using the Benjamini-Hochberg approach
    padj.interaction <- p.adjust(pvals.interaction, method="BH")
    #Look at the genes that have a significant adjusted p-value
    fit.full[(padj.interaction <= 0.05) & !is.na(padj.interaction),]

Note that the estimates are already log2 transformed counts, and that by default DESeq (and the glm() family of functions in R) use a "treatment contrast by default. So the genotype column represents a fold change relative to the intercept (in this case for the "ORE" & "wt" treatment combinations). If you had an intercept estimate of 5.2 for a given gene you could just use::

    2^5.2

Which would estimate the number of counts for ORE wild type flies.   
    
We could do a more extreme comparison between the full model and the null model::
    pvals.fullnull <- nbinomGLMTest(fit.full, fit.null)
    padj.fullnull <- p.adjust(pvals.fullnull, method="BH")
    fit.full[(padj.fullnull <= 0.05) & !is.na(padj.fullnull),]
    #And so on...

It is also very useful to do some plotting both for QC and for examining the totatily of your data. There are many different graphical approaches to examining your data, which we do not have time to get into here. For now we will just do an MA-plot and a volcano plot as a quick starting point. We will use these to compare the mutant (sdE3) from the wild type (i.e. the mutational treatment)::
    
    pvals.mutant <- nbinomGLMTest(fit.genotype, fit.null)
    padj.interaction <- p.adjust(pvals.interaction, method="BH")
    fit.full[(padj.interaction <= 0.05) & !is.na(padj.interaction),]

For the volcano plot we have fold change on the X-axis and -log10 of the p-value on the y-axis. First we create a dataframe containing the two relevant columns (the fold change from the treatment contrast, and the unadjusted p values)::

    volcano_mutant <- cbind(pvals.mutant, fit.genotype[,2])

Now we can plot this (I am going to produce the same plot but zoom in the second one)::

    plot(x=volcano_mutant[,2], y=-log10(volcano_mutant[,1]), xlab = "FOLD CHANGE", ylab="-log10(p)")

You probably notice the weird points that are extreme for fold change, but with low p-values. These are exactly the reason you always need to check your data graphically. What might be going on with these points? We could subset the data based on fold changes to pull out those genes. However, R has a reasonably useful function that might help for these cases :code:`identify()`. Once it is called scroll your mouse over the plot and click to highlight points. The numbers that appear are the index. Once you have finished press escape::

    identify(x=volcano_mutant[,2], y=-log10(volcano_mutant[,1]))

Remember to press :code:`esc` (or right click)!!! I have highlighted one point and it returned a value of 3713. So I will go back to the original count data and take a look at that row::
 
    raw.deseq.data[3713,]

Aha! This is a gene with very few counts (mostly zeroes), and we forgot to exclude such genes! In general if the mean number of counts (for a given gene) are below some threshold (say 5 or 10) you should probably exclude that gene since you have sampled so poorly from it that it would not be meaningful.

For now though, we will just zoom in to take a look::
 
    plot(x=volcano_mutant[,2], y=-log10(volcano_mutant[,1]), xlab = "FOLD CHANGE", ylab="-log10(p)", xlim=c(-5,5))
    abline

We can also fit an MA-plot, comparing the mean expression level of a gene with its fold change. As I am lazy, I did not actually compute the mean itself, but used the "intercept" (which represents the mean for the ORE wild type flies). However, this is likely sufficient for this plot::

    plot(x=fit.genotype[,1], y=fit.genotype[,2], xlim=c(0,20), ylim=c(-5,5),
    xlab= " ~ mean log2(counts)", ylab=" fold change", cex.lab=3)

We can also do some exploratory plotting to see if there's anything weird going on in our data. For example, we want to make sure that, if we cluster our samples based on overall gene expression patterns, we see clusters based on biological attributes rather than, say, the lane they were sequenced in or the day that the library preps were performed (the latter would indicate that there might be some unaccounted variable in our experimental design that is influencing gene expression)::

    #First, get dispersion estimates "blindly", i.e., without taking into account the sample treatments
    cdsFullBlind = estimateDispersions( deseq.data, method = "blind" )
    vsdFull = varianceStabilizingTransformation( cdsFullBlind )
    
    library("RColorBrewer")
    library("gplots")
    # Note: if you get an error message when you try to run the previous two lines,
    # you may need to install the libraries, like this:
    install.packages("RColorBrewer")
    install.packages("gplots")
    # After the libraries installed, don't forget to load them by running the library() calls again


Now let's make some heat maps::

    select = order(rowMeans(counts(deseq.data)), decreasing=TRUE)[1:30]
    hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
    
    # Heatmap of count table -- transformed counts
    heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
    
    # Heatmap of count table -- untransformed counts; you can see this looks pretty different
    # from the first heat map
    heatmap.2(counts(deseq.data)[select,], col = hmcol, trace="none", margin=c(10,6))
    
    # Heatmap of sample-to-sample distances
    dists = dist( t( exprs(vsdFull) ) )
    mat = as.matrix( dists )
    rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), paste(background, genotype, sep=" : "))
    heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))


We can also do a principal components analysis (PCA)::
    
    print(plotPCA(vsdFull, intgroup=c("background", "genotype")))
    
In this case, we see not only that the samples cluster by genotype and genetic background but also that PC1 represents genetic background, and PC2 seems to represent wt vs. mutant genotype.