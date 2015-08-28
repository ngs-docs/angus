================================================
Bacterial genome annotation using Prokka
================================================

After you have de novo assembled your genome sequencing reads into contigs,
it is useful to know what genomic features are on those contigs. The process
of identifying and labelling those features is called genome annotation.

In this tutorial you will:

1. Download and install Prokka
2. Annotate a FASTA file of contigs
3. Visualize the annotation using Artemis

The instructions below will work on a Linux server (eg. EC2 instance),
a Linux desktop, and even directly on your Mac OS X Laptop.

Download and install Prokka
===========================

Prokka is simple to install because it comes bundled with all its dependencies:

::

  git clone https://github.com/tseemann/prokka.git
  export PATH=$PWD/prokka/bin:$PATH
  prokka --setupdb
  prokka --version

Get some contigs
================

We will download a genome from NCBI, decompress it, and rename it to something shorter:

::

  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000144955.1_ASM14495v1/GCA_000144955.1_ASM14495v1_genomic.fna.gz
  gunzip GCA_000144955.1_ASM14495v1_genomic.fna.gz
  mv GCA_000144955.1_ASM14495v1_genomic.fna contigs.fasta

We can use the "grep" command to look at the FASTA sequence names in the file:

::

  grep '>' contigs.fasta

How many sequences/contigs are in this file?

We can use the "wc" (word count) command to get a rough idea of the number of basepairs in the FASTA file too
by counting how many characters (bytes) are in the file, as it uses 1 charactert (A,G,T,C) per nucleotide.

::

  wc -c contigs.fasta

How big is the genome in Mbp

Run Prokka on the contigs
=========================

Prokka is a pipeline script which coordinates a series of genome feature predictor tools and sequence similarity
tools to annotate the genome sequence (contigs).

::

  prokka --outdir anno --prefix prokka contigs.fasta

::

  cat ./anno/prokka.txt

How many genes did Prokka find in the contigs?


Install Artemis
===============

Artemis is a graphical Java program to browse annotated genomes.
It is a a bit like IGV but sepcifically designed for bacteria.
You will need to install this on your desktop computer.
You could run it remotely over SSH using X11 forwarding from Amazon
but it is probably too slow to be useful.

Download: https://www.sanger.ac.uk/resources/software/artemis/#download


Viewing the annotated genome
============================

* Start Artemis
* Go to "File" -> "SSH File Manager"
* Type in the IP number of your Amazon EC2 instance
* Browse to the "anno/prokka.gff" file
