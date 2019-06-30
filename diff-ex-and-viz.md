# Differential Expression and Visualization in R

Learning objectives:

* Create a gene-level count matrix of Salmon quantification using tximport
* Perform differential expression of a single factor experiment in DESeq2
* Perform quality control and exploratory visualization of RNA-seq data in R

## Getting started on Jetstream

[Start up an m1.medium instance running Ubuntu 18.04 on Jetstream.](jetstream/boot.html)

You should have the salmon counts on your instance in the `~/quant/` folder. If
not, you can get this data by running:

```
curl -L -o
tar xvf ...
```

## Importing gene-level counts into R using tximport

When we finished running salmon, we had a `quant` folder that contained 
transcript-level read counts. So far in our RNA-seq pipeline, we have been
working from the command line in bash. For the rest of the RNA-seq tutorial, 
we will be working in R. As we saw in the R introduction lesson, R is really
good and working with data in table format. When we produced counts for our
reads, we essentially transformed our data to this format. Most of the software
tools written to analyze RNA-seq data in this format are written in R. 

We first need to read our data into R. To do that, we will use a package called
[`tximport`](https://bioconductor.org/packages/release/bioc/html/tximport.html).
This package takes transcript-level counts and summarizes them to the gene 
level. It is compatible with many count input formats, including salmon. 

First, launch RStudio 

ADD LINK TO INSTRUCTIONS ON HOW TO LAUNCH RSTUDIO

Next, load the library `tximport` so we have access to the functions. 

```
library(tximport)
```

We need a couple files to import our counts. 

First, we need a file that tells `tximport` where our files are. We've added
other information to this file about sample names and condition, which we will
also use later for differential expression. 

```
download.file(url = , dest.file = "")
```

Second, we need a file that specifies which transcripts are associated with
which genes. In this case, we have associated transcript names with yeast
ORF names. 

```
download.file(url = , dest.file = "")
```

`
## Further Notes

### Making a tx2gene file

Making a tx2gene file is often a little different for each organism. If your 
organism has a transcriptome (or \*rna\_from\_genomic.fnai.gz file) on RefSeq,
you can often get the information from the `*feature_table.txt.gz`. You 
might also be able to parse a gtf or gff file to produce the information you
need. This information is sometimes also found in the fasta headers of the 
transcriptome file itself.
