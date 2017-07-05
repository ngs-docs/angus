# Analyzing ChIP-seq data
_(the material here is based on the *ChIP­‐seq Hands­‐on Exercise*, by Remco Loos and Myrto Kostadima at EMBL‐EBI, available [here](https://www.ebi.ac.uk/training/online/course/ebi-next-generation-sequencing-practical-course/chip-seq-analysis/chip-seq-practical)_

## What is ChIP-seq?



## Our goal

The goal of this lesson is to perform some basic tasks in the analysis of ChIP-seq data. The first step includes an unspliced alignment for a small subset of raw reads. We will align raw sequencing data to the mouse genome using Bowtie and then we will manipulate the SAM output in order to visualize the alignment on the IGV browser. Then based on these aligned reads we will find immuno-enriched areas using the peak caller MACS. We will then perform functional annotation and motif analysis on the predicted binding regions.

## Get some sample data

We will download the bundled data directly from the EMBL-EBI exercise [here](https://www.ebi.ac.uk/~emily/Online%20courses/NGS/ChIP-seq.zip)

The actual data that we will use for the ChIP-seq workflow are reported in Chen, X et al. (2008), [Integration of external signaling pathways with the core transcriptional network in embryonic stem cells](http://www.sciencedirect.com/science/article/pii/S009286740800617X). Cell. Jun 13;133(6):1106-17. It contains a few reads from the Mouse genome, and we will try to identify potential transcription factor binding sites of Oct4 in mouse embryonic stem cells.

Lets download these data using `wget` and uncompress using `unzip` as follows:

```
cd ~
wget https://www.ebi.ac.uk/~emily/Online%20courses/NGS/ChIP-seq.zip
unzip ChIP-seq.zip
```

You should be able to see a new folder now in your home directory, with the name `ChIP-seq`, which contains both our reads (files `Oct4.fastq` and `gfp.fastq`), as well as some other files that will be useful down the road.

> Question: How many reads are in these files? How can you find out?


## Let's set our tools





## Reference

* Bailey T, Krajewski P, Ladunga I, et al. [Practical Guidelines for the Comprehensive Analysis of ChIP-seq Data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3828144/). Lewitter F, ed. _PLoS Computational Biology_. 2013;9(11):e1003326. doi:10.1371/journal.pcbi.1003326.
