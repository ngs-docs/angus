## Multivariate tests with NGS data

The tests I will present here can be used universally across multivariate data.

First, either download the data or get it directly from the datasets directory using a function from the RCurl package.

Let's download and install the package first
```R
install.packages("RCurl")
library(RCurl)
```
Using a combination of functions from R, we will create an object called dataset that holds all of the fly data
```R
URL<-("https://raw.githubusercontent.com/ryanjw/ngs-3rdweek/master/datasets/fly_data.txt")
dataset<-read.table(textConnection(getURL(URL)),header=T,check.names=F,sep="\t")
```
If you are unable to do download RCurl, [run this code instead](https://github.com/ryanjw/ngs-3rdweek/blob/master/multivariate-tests/alternative-download.md)

Let's look at the first few rows and columns of the dataset.
```R
head(dataset[,1:10])
```
The output should look like this
```R
  fly type                           file FBgn0000003 FBgn0000008 FBgn0000014 FBgn0000015 FBgn0000017 FBgn0000018 FBgn0000022
1 HYB sdE3 HYB_sdE3_rep1_htseq_counts.txt           0           0           0           0           0           0         200
2 HYB sdE3 HYB_sdE3_rep2_htseq_counts.txt           0           0           0           0           0           0         319
3 HYB   wt   HYB_wt_rep1_htseq_counts.txt           0           0           0           0           0           0         509
4 HYB   wt   HYB_wt_rep2_htseq_counts.txt           0           0           0           0           0           0         331
5 ORE sdE3 ORE_sdE3_rep1_htseq_counts.txt           0           0           0           0           0           0         385
6 ORE sdE3 ORE_sdE3_rep2_htseq_counts.txt           0           0           0           0           0           0         312
```
We can see several columns describing what fly strain, mutant type, and other metadata regarding these samples.  Starting in column three, we can see the counts of specific contigs from this RNAseq experiment.

### Multivariate analysis of variance

There are several options to determine whether or not experimental designs produce significantly different results.  A classic way of doing this is using multivariate analysis of variance or MANOVA.  A major assumption of MANOVA is that variables within each treatment level come from a multivariate normal distribution.

To run a MANOVA, we set up a linear model first and then nest a few functions together.  Here we demonstrate a Wilk's lambda, which is common in the literature.  Note that there are several options; run `?summary.manova` to explore them.  **This code may be too intense for your laptop, consider using EC2**. 
```R
test<-lm(as.matrix(dataset[,-c(1:3)])~dataset$fly)
summary(manova(test),test="Wilks")
``` 
If you run this code, it may not work...see output below
```R
Error in qr.default(D %*% ss[[nt]] %*% D, tol = tol) :
  NA/NaN/Inf in foreign function call (arg 1)
 ```
This code ends up failing, but why may this be?  Is our data violating some assumptions?  Let's grab a few random variables from the dataset and see how they are distributed.
```R
library(ggplot2)
ggplot(dataset)+geom_density(aes(x=FBgn0000022,fill=fly,alpha=0.5))
```
![alt text](https://raw.githubusercontent.com/ryanjw/ngs-3rdweek/master/multivariate-tests/fly-density-plot.jpg)
##Challenge
*Try looking at another variable by replacing `FBgn0000022` with another variable name and coloring by `type` instead of `fly`*

##Trying nonparametric tests

A problem with these types of data is that they are often not normally distributed and it is common to have very few samples.  Multivariate tests that run in a nonparametric framework (i.e. using permutations) are a powerful alternative.  *Note that these can be run on your laptop, while MANOVA may not!*

We will start with a permutational multivariate analysis of variance (PERMANOVA)
```R
# install.packages("vegan")
library(vegan)
?adonis
adonis(dataset[,-c(1:3)]~dataset$fly*dataset$type)
```
The output should look like
```R
Call:
adonis(formula = dataset[, -c(1:3)] ~ dataset$fly * dataset$type) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                         Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)
dataset$fly               2  0.031265 0.0156327 1.10720 0.21511  0.390
dataset$type              1  0.001778 0.0017776 0.12590 0.01223  0.959
dataset$fly:dataset$type  2  0.027591 0.0137957 0.97709 0.18983  0.459
Residuals                 6  0.084715 0.0141191         0.58284       
Total                    11  0.145349                   1.00000    
```
Check your neighbor to see if you have the same numbers.  While there is not a significant effect, it is important to note that `fly` explains the most variation, while the interaction between `fly` and `type` is also strong.  With more samples, these relationships may resolve.

This function, `adonis`, is very strong.  It allows you to incorporate block effects by passing arguments to `strata`.  

We can also use distance based approaches to ask more specific multivariate questions like, *"Does the composition of transcripts differ between flys and types?"* or *"Are different transcripts present or absent between flys and types?"*

We will do this by nesting another function, `decostand`, within the `adonis` function

```R
adonis(decostand(dataset[,-c(1:3)],method="total")~dataset$fly*dataset$type)

Call:
adonis(formula = decostand(dataset[, -c(1:3)], method = "total") ~      dataset$fly * dataset$type) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                         Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)  
dataset$fly               2 0.0086981 0.0043491  1.9775 0.30315  0.020 *
dataset$type              1 0.0022970 0.0022970  1.0444 0.08006  0.409  
dataset$fly:dataset$type  2 0.0045017 0.0022508  1.0235 0.15690  0.447  
Residuals                 6 0.0131955 0.0021993         0.45990         
Total                    11 0.0286923                   1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Note that that by taking a distance based approach we can actually gain resolution.  Further, the amount of variation explained is qualitatively similar to the non-distance based approach.

##Challenge

Try to do the same thing based on presence/absence by changing `total` to `pa` and adding `method="jaccard"` to the `adonis` function.  Be careful, and watch the placement of `()` when nesting these functions!

[*You can find the solution here*](https://github.com/ryanjw/ngs-3rdweek/blob/master/multivariate-tests/pa-soln.md)

#Comparing matrices with mantel tests

Another useful method is the Mantel Test.  Here, two distance matrices are compared by looking at correlations between entries and then permuting entries to determine significance.  
```R
?mantel
unique(dataset$fly)

HYB_subset<-subset(dataset, fly=="HYB")[,-c(1:3)]
ORE_subset<-subset(dataset, fly=="ORE")[,-c(1:3)]

HYB_dist<-vegdist(decostand(HYB_subset,"pa"),method="jaccard")
ORE_dist<-vegdist(decostand(ORE_subset,"pa"),method="jaccard")

mantel(ORE_dist,HYB_dist,method="spearman",permutations=9999)
```

##Challenge

Try to do all pairwise comparisons between fly types.  Are the results as sensitive as the PERMANOVA?
[*You can find the solution here*](https://github.com/ryanjw/ngs-3rdweek/blob/master/multivariate-tests/pairwise-mantel.md)



