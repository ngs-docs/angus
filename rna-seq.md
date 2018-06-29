# # RNAseq

TODO:

* Put learning objectives

Learning objectives:

* ;
  
* ;

* ;

* ;

## Boot up a Jetstream

[Boot an m1.medium Jetstream instance](jetstream/boot.md) and log in.

## Install software

We will be using salmon

```
conda install -y salmon
```

## Change to a new working directory and link the original data

We will be using the same data as before ([Schurch et al, 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/)), so the following commands will create a new folder `rnaseq` and link the data in:

```
cd ~/
mkdir -p rnaseq
cd rnaseq

ln -fs ~/data/*.fastq.gz .
ls
```

## Download the yeast reference transcriptome:


```
curl -O https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
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

Read up on [libtype, here](https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype).

##  Collect all of the sample counts using [this Python script](https://github.com/ngs-docs/2018-ggg201b/blob/master/lab6-rnaseq/gather-counts.py):

```
curl -L -O https://raw.githubusercontent.com/ngs-docs/2018-ggg201b/master/lab6-rnaseq/gather-counts.py
python2 gather-counts.py
```

##  Install edgeR

```
curl -L -O https://raw.githubusercontent.com/ngs-docs/2018-ggg201b/master/lab6-rnaseq/install-edgeR.R
sudo Rscript --no-save install-edgeR.R
```

##  Run edgeR (in R) using [this script](https://github.com/ngs-docs/2018-ggg201b/blob/master/lab6-rnaseq/yeast.salmon.R) and take a look at the output:

```
curl -L -O https://raw.githubusercontent.com/ngs-docs/2018-ggg201b/master/lab6-rnaseq/yeast.salmon.R
Rscript --no-save yeast.salmon.R
```

This will produce two plots, `yeast-edgeR-MA-plot.pdf` and
`yeast-edgeR-MDS.pdf`. You can view them by going to your RStudio server file viewer, changing to  the directory `yeast`, and then clicking on them.

 The `yeast-edgeR.csv` file contains the fold expression & significance information in a spreadsheet.

## Questions to ask/address

1. What is the point or value of the [multidimensional scaling (MDS)](https://en.wikipedia.org/wiki/Multidimensional_scaling) plot?

2. Why does the MA-plot have that shape?

   Related: Why can't we just use fold expression to select the things we're interested in?

   Related: How do we pick the FDR (false discovery rate) threshold?

3. How do we know how many replicates (bio and/or technical) to do?

   Related: what confounding factors are there for RNAseq analysis?

   Related: what is our false positive/false negative rate?
   
4. What happens when you add new replicates?

## More reading

"How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?" [Schurch et al., 2016](http://rnajournal.cshlp.org/content/22/6/839).

"Salmon provides accurate, fast, and bias-aware transcript expression estimates using dual-phase inference" [Patro et al., 2016](http://biorxiv.org/content/early/2016/08/30/021592).

Also see [seqanswers](http://seqanswers.com/) and [biostars](https://www.biostars.org/).