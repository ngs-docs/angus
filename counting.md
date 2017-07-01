# RNAseq expression analysis

Some older tutorials and presentations:

* [Degust: RNA-seq exploration, analysis, and visualisation](http://degust.erc.monash.edu/)
* [Routes through differential expression](https://angus.readthedocs.io/en/2016/_static/RoutesThroughDESeqAnalysis.pdf)
* [2016 tutorial!](https://angus.readthedocs.io/en/2016/rob_quant/tut.html)

During this lesson, you’ll learn how to use salmon to rapidly quantify
transcript-level expression from RNA-seq data.

## Make sure R & RStudo are installed:

```
sudo apt-get update && sudo apt-get -y install gdebi-core r-base
```

Try to connect to a running RStudio Web server instance – you can get the Web address by running this command:
```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

If you cannot connect, download and install RStudio.
```
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi -n rstudio-server-1.0.143-amd64.deb 
```
And, finally, change the password to something you can remember. If your username is different than the one below (i.e. `diblions`), you'll need to change that.
```
sudo passwd tx160085
```      

## Install edgeR

Use [this script](https://github.com/ngs-docs/angus/blob/change_link_for_edgeR_script/_static/install-edgeR.R):

```
cd
curl -O -L https://github.com/ngs-docs/angus/raw/change_link_for_edgeR_script/_static/install-edgeR.R
sudo Rscript --no-save install-edgeR.R
```

## Install [salmon](https://salmon.readthedocs.io):

```
cd ~/
curl -L -O https://github.com/COMBINE-lab/salmon/releases/download/v0.8.2/Salmon-0.8.2_linux_x86_64.tar.gz
tar xzf Salmon-0.8.2_linux_x86_64.tar.gz
export PATH=$PATH:$HOME/Salmon-0.8.2_linux_x86_64/bin
```

## Change to the appropriate directory:

```
mkdir yeast
cd yeast
```

## Download some data

From [Schurch et al, 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/) massive study on yeast!

This may take a little while...

```
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458500/ERR458500.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458501/ERR458501.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458502/ERR458502.fastq.gz

curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458493/ERR458493.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458494/ERR458494.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458495/ERR458495.fastq.gz
```

## Download the yeast reference transcriptome:

```
curl -O http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
```

## Index the yeast transcriptome:

```
salmon index --index yeast_orfs --type quasi --transcripts orf_coding.fasta.gz
```

## Run salmon on all the samples:

```
for i in *.fastq.gz
do
    salmon quant -i yeast_orfs --libType U -r $i -o $i.quant --seqBias --gcBias
done
```

Here is a description of the parameters used, for more information refere to the [salmon webpage](http://salmon.readthedocs.io/en/latest/salmon.html)
  - -i tells salmon where to look for the index
  - -p tells salmon how many threads to use
  - -l tells salmon the type of the read library (here, inward facing, unstranded reads). For a more in-depth description of the library types and how to specify them in salmon, have a look here in the docs.
  - -1 this tells salmon where to find the first reads of the pair
  - -2 tells salmon where to find the second reads of the pair
  - -o tells salmon where (the directory) to write the output for this sample. The directory (and the path to it) will be created if it doesn’t exist.
   
What do you think all this stuff with the bias is about?

Read up on [libtype, here](https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype).

## Collect all of the sample counts

Salmon outputs things into subdirectories and in a format that is inconvenient
for R; use [this Python script](https://github.com/ngs-docs/angus/blob/change_link_for_edgeR_script/_static/gather-counts.py) to collect them all!

```
curl -L -O https://raw.githubusercontent.com/ngs-docs/angus/change_link_for_edgeR_script/_static/gather-counts.py
python2 gather-counts.py
```

## Run edgeR (in R)

USe [this script](https://github.com/ngs-docs/angus/blob/change_link_for_edgeR_script/_static/yeast.salmon.R) and take a look at the output:

```
curl -L -O https://github.com/ngs-docs/angus/raw/change_link_for_edgeR_script/_static/yeast.salmon.R
Rscript --no-save yeast.salmon.R
```

This will produce two plots, `yeast-edgeR-MA-plot.pdf` and
`yeast-edgeR-MDS.pdf`. You can view them by going to your RStudio Server
console and looking in the directory `yeast`.

The `yeast-edgeR.csv` file contains the fold expression & significance information in a spreadsheet.

## Extra plotting in R

The plots we made above are nice, but what if we want something a bit more informative?

We can use the skills we learned from the packages dplyr and ggplot2 in order to make a colored plot with gene names. Let's start by installing and loading the necessary packages. 
```r
# set our working directory
setwd("~/yeast/")

# install packages
install.packages("ggplot2")
install.packages("dplyr")

# load the packages
library(ggplot2)
library(dplyr)
```
Next, we can read in our data and get a feel for what the data look like.
```r
results = read.csv("yeast-edgeR.csv")
head(results)
dim(results)
```

We want to produce a plot that will have points that are colored by significance, and that has the most significant genes labeled. First, let's add a column with the mutate() function in dplyr that specifies the level of significance using the false discovery rate (FDR) column.
```r
# mutate with dplyr, adding a column named "significance")
results = mutate(results, significance=ifelse(results$FDR<0.05, "FDR < 0.05", "FDR > 0.05"))
```
Then we will construct a plot with ggplot. We will plot log fold change against -log10 of P value.
```r
# construct a plot
plot = ggplot(results, aes(logFC, -log10(PValue))) +
  geom_point(aes(col=significance)) +
  scale_color_manual(values=c("red", "black"))
 
# View the plot
plot
```
This is a great start! This is a classic ggplot graph. We see that we have a legend indicating what our colors mean, and that we can start to get a feel for the distribution of significance of differentially expressed genes. 

But this plot can still be improved. What if we wanted to know the names of the *most* significantly expressed genes? We could do that by adding a layer to our ggplot using the function geom_text(). Our gene column was named "X", and here we are choosing to label the genes that have a significance value less than 1e-200. You can see that in order to do this, we use the dplyr command filter(). 
```r
# Add gene labels to points that are highly significant
plot + geom_text(data=filter(results, FDR<1e-200), aes(label=X))
```

Alas, this is R, so we still have a problem. Our text is overlapping and fairly difficult to read. Luckily, other people have ran into this problem and they have solved it! There is a package called ggrepel (gg-repel) that "Provides text and label geoms for 'ggplot2' that help to avoid overlapping text labels. Labels repel away from each other and away from the data points. This is exactly what we need! Let's install the package and use it to stagger our labels. We use geom_text_repel instead of geom_text, but we still use the dplyr filter() function.
```r
# Install ggrepel package
      # Provides text and label geoms for 'ggplot2' that help to avoid overlapping text labels.   
      # Labels repel away from each other and away from the data points.
install.packages("ggrepel")
library(ggrepel)

# add labels for genes that are significant
plot + geom_text_repel(data=filter(results, FDR<1e-200), aes(label=X))
```

## Questions to ask/address

1. What is the point or value of the [multidimensional scaling (MDS)](https://en.wikipedia.org/wiki/Multidimensional_scaling) plot?

2. Why does the MA-plot have that shape?

   (This is really important to address - basic statistical noise gives it
   that shape!)

   Related: Why can't we just use fold expression to select the things we're interested in?
   
   (Two reasons: multiple comparison problems; and standard error in measurement. If avg expression is low, then "large" relative deviations may occur frequently; but this won't happen as much for high expression genes.)

   Related: How do we pick the FDR (false discovery rate) threshold?
   
   (Somewhat arbitrarily.  Your FDR threshold essentially gives you a
   tradeoff between sensitivity (you miss more when you use a lower FDR
   threshold) vs specificity (if you raise the FDR threshold, you will
   accept genes that look differentially expressed due solely to measurement
   error/noise.)
   
3. How do we know how many replicates (bio and/or technical) to do?

   GOOD QUESTION. Apart from statistical issues, what biology features in to
   this?

   Related: what confounding factors are there for RNAseq analysis?

   Related: what is our false positive/false negative rate? How would you
   imagine calculating those?

## More reading

"How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?" [Schurch et al., 2016](http://rnajournal.cshlp.org/content/22/6/839).

"Salmon provides accurate, fast, and bias-aware transcript expression estimates using dual-phase inference" [Patro et al., 2016](http://biorxiv.org/content/early/2016/08/30/021592).

Also see [seqanswers](http://seqanswers.com/) and [biostars](https://www.biostars.org/) for general discussions :).

