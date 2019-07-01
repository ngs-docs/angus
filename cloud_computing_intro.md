# Intro to cloud computing

## Learning Objectives:
+ Understand what cloud computing is?
+ Get to know some popular cloud compute resources
+ Perform command line blast
+ 

## Rationale

Bioinformatics worldwide involves researchers connect daily to cloud computing services to perform their analyses.

There are several reasons for this. Firstly biology - like most other areas of science - is dealing with a deluge of data due to the rapid advancement of data collection methods. It is now common that data collected for an experiment doesn't fit on a researcher's laptop and that the resources needed for running an analysis far exceed a desktop computer's computing power.

Secondly the vast majority of research software are developed and released for linux. Most people run MacOS or Windows on their desktop computers and laptop, which makes the installation of some software difficult or at the very least inconvenient.

## What is the Cloud?

The cloud is basically lots of servers (thing big big computers) stacked together in a giant, powerful infrastructure. You can lend part of this infrastructure for your computing needs. While it is not cheap, it is generally scalable and guarantees a stable environment.

In research there are two approaches to lend computing time and power: either (a) you lend computing time and resources from a commercial provider or you obtain access to a research computing infrastructure. Some countries have built national infrastructures where you can apply for computing time for your research projects. Most academic institutions or departments also have their own computing resources.

![](static/cloud_compute.png)

## Popular Cloud/HPC resources
+ Academic
    + [Jetstream](https://www.jetstream-cloud.org/)
    + [Atmosphere](https://atmo.cyverse.org/)

+ Commercial
    + [AWS](https://aws.amazon.com/getting-started/tutorials/launch-a-virtual-machine/)
    + [Google Cloud](https://cloud.google.com/)

## Let's connect to the cloud


open jupyter notebook
http://149.165.171.170:8888/lab

open terminal 


To install blast:
`conda install -c  bioconda blast`

Note: takes 5-10 minutes


# Command Line BLAST
## Basic Local Alignment Search Tool

Given one or more query sequences (usually in FASTA format), BLAST looks for matching sequence regions between them and a subject set.
![](https://i.imgur.com/YPVewm4.png)
A sufficiently close match between subsequences (denoted by arrows in the figure above, though matches are usually longer than illustrated here) is called a high-scoring pair (HSP), while a query sequence is said to hit a target sequence if they share one or more HSPs

# Running command-line BLAST

The goal of this tutorial is to run you through a demonstration of the
command line, which you may not have seen or used much before.

## Running BLAST

First! We need some data (protein sequences) and an Rscript to later visualize blast results.  Let's grab the mouse and zebrafish RefSeq
protein data sets from NCBI, and put them in our home directory. If you've just logged
in, you should be there already, but if you're unsure, run `cd` and hit enter. Now,
we'll use `curl` to download the files; these originally came from the NCBI Web site: [ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot](ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot).

```
curl -o mouse.1.protein.faa.gz -L https://osf.io/v6j9x/download
curl -o mouse.2.protein.faa.gz -L https://osf.io/j2qxk/download
curl -o zebrafish.1.protein.faa.gz -L https://osf.io/68mgf/download
curl -o blastviz.R -L https://osf.io/e548g/download
```

To look at the files in your current directory:

```
ls -l
```

All three of the files are FASTA protein files (that's what the .faa
suggests) that are compressed with `gzip` (that's what the .gz means).

Uncompress them:

```
gunzip *.faa.gz
```

and let's look at the first few sequences in the file:

```
head mouse.1.protein.faa 
```

These are protein sequences in FASTA format.  FASTA format is something
many of you have probably seen in one form or another -- it's pretty
ubiquitous.  It's a text file, containing records; each record
starts with a line beginning with a '>', and then contains one or more
lines of sequence text.

Let's take those first two sequences and save them to a file.  We'll
do this using output redirection with '>', which says "take
all the output and put it into this file here."

```
head -11 mouse.1.protein.faa > mm-first.fa
```

So now, for example, you can do `cat mm-first.fa` to see the contents of
that file (or `less mm-first.fa`).

Now let's BLAST these two sequences against the entire zebrafish
protein data set. First, we need to tell BLAST that the zebrafish
sequences are (a) a database, and (b) a protein database.  That's done
by calling 'makeblastdb':

```
makeblastdb -in zebrafish.1.protein.faa -dbtype prot
```

Next, we call BLAST to do the search:

```
blastp -query mm-first.fa -db zebrafish.1.protein.faa
```

This should run pretty quickly, but you're going to get a lot of output!!
To save it to a file instead of watching it go past on the screen,
ask BLAST to save the output to a file that we'll name `mm-first.x.zebrafish.txt`:

```
blastp -query mm-first.fa -db zebrafish.1.protein.faa -out mm-first.x.zebrafish.txt
```

and then you can 'page' through this file at your leisure by typing:

```
less mm-first.x.zebrafish.txt
```

(Type spacebar to move down, and 'q' to get out of paging mode.)

-----

Let's do some more sequences (this one will take a little longer to run):

```
head -500 mouse.1.protein.faa > mm-second.fa

blastp -query mm-second.fa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.txt
```

will compare the first 83 sequences.  You can look at the output file with:

```
less mm-second.x.zebrafish.txt
```

(and again, type 'q' to get out of paging mode.)

Notes:

* you can copy/paste multiple commands at a time, and they will execute in order;

* why did it take longer to BLAST ``mm-second.fa`` than ``mm-first.fa``?


----

Last, but not least, let's generate a more machine-readable version of that
last file --

```
blastp -query mm-second.fa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.tsv -outfmt 6
```

See [this link](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) for a description of the possible BLAST output formats.


Run the following R script to visualize the blast results as below

```
Rscript blastviz.R
```

A pdf will be generated with the results. Double click to open the pdf 

![](static/blastviz.png)

Things to mention and discuss:

* `blastp` options and -help.
* command line options, more generally - why so many?
