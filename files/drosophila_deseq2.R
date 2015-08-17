#Install DESeq2
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

#Load the DESeq2 library
library(DESeq2)

#Generate a "data frame" to store the information about each of the samples in your
# RNA-seq experiment
#We could also have done this by creating a CSV or tab-delimited text file with
# all of this information and loading it into R with read.csv() or read.table()
# [That would probably be easier, for example, if we had a complicated design with
#  many samples and multiple variables]

samples <- c("ORE_wt_rep1","ORE_wt_rep2","ORE_sdE3_rep1","ORE_sdE3_rep2","SAM_wt_rep1","SAM_wt_rep2","SAM_sdE3_rep1","SAM_sdE3_rep2","HYB_wt_rep1","HYB_wt_rep2","HYB_sdE3_rep1","HYB_sdE3_rep2")

files <- paste(samples, "_htseq_counts.txt", sep="")

backgrounds <- c(rep("ORE", 4), rep("SAM", 4), rep("HYB", 4))

genotypes <- c(rep(c("wt", "wt", "sdE3", "sdE3"), 3))

rna.design <- data.frame(sample=samples, file=files, background=backgrounds, genotype=genotypes)

#DESeq2 needs a model to load your data; let's create a simple formula to start with, to look for
# genes that are differentially expressed between two genotypes [wild-type and scalloped-E3 mutants],
# without accounting for genetic background
load.model <- formula(~ genotype)

#Now load the data into R
all.data <- DESeqDataSetFromHTSeqCount(sampleTable=rna.design, directory="./counts", design=load.model)

#This next line will handle several steps for us--it will first estimate size factors, to account
# for differences in library size (total numbers of reads) across samples
#Then it will generate dispersion estimates
all.data <- DESeq(all.data)

#Always a good idea to plot dispersion estimates to make sure they look ok
plotDispEsts(all.data)

#Now let's look at differentially expressed genes
genotype.results <- as.data.frame(results(all.data, alpha=0.05))

#Viewing the top of the results data frame reveals that there's a lot of genes
# with missing data [no mapped reads]
#That makes sense for this sample dataset, because we kept only reads derived from
# X-linked genes [to speed up run times for the workshop]
head(genotype.results)

#Let's pull out the genes that are significant at a false discovery rate of 0.05
sig.genotype.results <- genotype.results[(genotype.results$padj <= 0.05) & !is.na(genotype.results$padj),]

#Sort them
sig.genotype.results <- sig.genotype.results[order(sig.genotype.results$padj, decreasing=F),]

#View them
sig.genotype.results

#You should see something like this:
               # baseMean log2FoldChange      lfcSE      stat       pvalue         padj
# FBgn0027287  363.507236      0.5801968 0.06718179  8.636222 5.810267e-18 7.919394e-15
# FBgn0030666  536.171769     -0.5022900 0.08275070 -6.069919 1.279752e-09 8.721507e-07
# FBgn0003345 4477.689176      0.4262786 0.07469411  5.706992 1.149905e-08 5.224402e-06
# FBgn0030608 7004.628284      0.4269717 0.07738605  5.517425 3.440022e-08 1.172187e-05
# FBgn0000022  381.986764      0.3922786 0.08078507  4.855830 1.198835e-06 3.268025e-04
# FBgn0030679  315.426730     -0.3511653 0.07424409 -4.729876 2.246575e-06 5.103469e-04
# FBgn0004170  481.871080      0.3638745 0.07765080  4.686037 2.785458e-06 5.423685e-04
# FBgn0027093 1576.040720     -0.2301959 0.04957300 -4.643573 3.424350e-06 5.834237e-04
# FBgn0003996 1234.285962     -0.2082758 0.04712886 -4.419284 9.902822e-06 1.499727e-03
# FBgn0030641   10.758415     -0.3002516 0.06895131 -4.354545 1.333437e-05 1.817474e-03
# FBgn0030696    8.133345     -0.2611802 0.06078021 -4.297126 1.730269e-05 2.143961e-03
# FBgn0028491   64.157746      0.3309560 0.08074447  4.098807 4.152845e-05 4.716940e-03
# FBgn0030628  301.998911      0.3490040 0.08995057  3.879953 1.044765e-04 1.017153e-02
# FBgn0030665   19.326091      0.2549655 0.06543734  3.896330 9.766126e-05 1.017153e-02
# FBgn0030611  467.278641     -0.3400626 0.09184314 -3.702645 2.133632e-04 1.938760e-02
# FBgn0029939  209.192726     -0.2998760 0.08210641 -3.652285 2.599174e-04 2.083926e-02
# FBgn0259923   18.873626      0.3151387 0.08609118  3.660523 2.517010e-04 2.083926e-02
# FBgn0029003 5364.243241     -0.2483548 0.06902117 -3.598241 3.203765e-04 2.425962e-02
# FBgn0030669  399.244296     -0.2763943 0.08155220 -3.389170 7.010455e-04 4.777625e-02
# FBgn0030912   53.831358     -0.3155388 0.09309440 -3.389450 7.003299e-04 4.777625e-02
# FBgn0027546 2228.528333     -0.2210875 0.06572049 -3.364058 7.680533e-04 4.985031e-02

#Note: FBgn0003345 is sd, which is the gene that is mutated in our non-control treatment group
# in this sample dataset, so it makes sense that it is at the top of our list of differentially
# expressed genes

#We can make a volcano plot
plot(x=genotype.results$log2FoldChange, y=-log10(genotype.results$padj), pch=16, col=ifelse(genotype.results$padj <= 0.05, "red", "black"))
#This one's not very interesting because we don't have too many points on the graph, since we used
# a reduced dataset


#Suppose we wanted to do this while accounting for genetic background; we could use
# more complex models. Let's re-load the data specifying a more complex model

load.model.2 <- formula(~ background + genotype)

all.data.2 <- DESeqDataSetFromHTSeqCount(sampleTable=rna.design, directory="./counts", design=load.model.2)

all.data.2 <- DESeq(all.data.2)

#Let's compare these new dispersion estimates to the old ones

#Divide the plot window into two side-by-side panels
par(mfrow=c(1, 2))

#Plot the original dispersion estimates on the left first, then the new dispersion estimates on the right
plotDispEsts(all.data)
plotDispEsts(all.data.2)

#Although the overall fit looks similar, it's pretty clear that for some genes, the dispersion
# estimates have changed from using this more complex model

#Now, let's compare this full model, which includes the effects of both genetic background and genotype,
# to a reduced model, which has only genetic background
#Then we can ask: for which genes is expression better explained by the more complex model?
# (using a likelihood ratio test)

reduced.model <- formula(~ background)

#Note that this line estimates dispersions again using our full and reduced models
all.data.2 <- DESeq(all.data.2, full=load.model.2, reduced=reduced.model, test="LRT")

#And now for the results
genotype.results.2 <- as.data.frame(results(all.data.2, alpha=0.05))

sig.genotype.results.2 <- genotype.results.2[(genotype.results.2$padj <= 0.05) & !is.na(genotype.results.2$padj),]

#Sort them
sig.genotype.results.2 <- sig.genotype.results.2[order(sig.genotype.results.2$padj, decreasing=F),]

#View them
sig.genotype.results.2

#Note that when we are comparing two models using a likelihood ratio test like this,
# the meaning of the log2FoldChange column is not immediately obvious
#We will need to use contrasts to compare two specific treatment groups
#Here, let's look at the overall effect of genotype
contrast.results.2 <- as.data.frame(results(all.data.2, contrast=c("genotype", "wt", "sdE3")))
contrast.results.2 <- contrast.results.2[(contrast.results.2$padj <= 0.05) & !is.na(contrast.results.2$padj),]
contrast.results.2 <- contrast.results.2[order(contrast.results.2$padj, decreasing=F),]




#We could also look at models with interaction effects
# (e.g., in this case, which genes are affected by the sdE3 mutation in a background-
#  dependent way?)

load.model.3 <- formula(~ background + genotype + background:genotype)
all.data.3 <- DESeqDataSetFromHTSeqCount(sampleTable=rna.design, directory="./counts", design=load.model.3)
all.data.3 <- DESeq(all.data.3)

#When we have complex models like this and we're interested in comparing two specific
# treatment groups, we need to use "contrasts"
#This can get rather complicated. We'll do one or two examples here, but see the DESeq2
# tutorial for more details

#Compare ORE-wt to ORE-sdE3

contrast.results.3 <- as.data.frame(results(all.data.3, contrast=list(	c("genotypewt", "backgroundORE.genotypewt") ,
																		c("genotypesdE3", "backgroundORE.genotypesdE3")
						)))
						
contrast.results.3 <- contrast.results.3[(contrast.results.3$padj <= 0.05) & !is.na(contrast.results.3$padj),]
contrast.results.3 <- contrast.results.3[order(contrast.results.3$padj, decreasing=F),]