=======================
Analyzing nanopore data
=======================

Last week in Woods Hole, MA we used our `lab's <http://ivory.idyll.org/lab/>`__ `MinION <https://www.nanoporetech.com/>`__ to sequence a new bacterial species isolated by `Rebecca Mickol <https://news.uark.edu/articles/27669/earth-organisms-survive-under-low-pressure-martian-condition>`__ in the `Microbial Diversity Course at the Marine Biological Lab <http://www.mbl.edu/microbialdiversity/>`__.

The goals of this tutorial are to:

*  convert oxford nanopore data in .fast5 format to .fastq
*  assess a nanopore run
*  assemble data
*  evaluate the assembly

Starting an AWS instance and installing software:
==================================================

Start a blank Amazon instance (m3.xlarge) with 50 GB and `log in <http://angus.readthedocs.io/en/2016/amazon/index.html>`__.

Copy/paste to update and install software on your new instance:
::
    sudo apt-get update && \
    sudo apt-get -y install build-essential ruby screen git curl gcc make g++ python-dev unzip \
        default-jre pkg-config libncurses5-dev r-base-core \
        r-cran-gplots python-matplotlib sysstat python-virtualenv \
        python-setuptools cmake cython libhdf5-serial-dev \
        python-numpy python-scipy python-pandas python-pandas-lib python-biopython parallel

We will now install several software packages that are specific for analyzing long reads data, as comes from the Oxford Nanopore MinION.

`poretools <http://poretools.readthedocs.io/en/latest/content/installation.html#basic-installation>`__
::
    sudo pip install poretools
    poretools

You should see output like this:
::
    usage: poretools [-h] [-v]
                     {combine,fastq,fasta,stats,hist,events,readstats,tabular,nucdist,metadata,index,qualdist,qualpos,winner,squiggle,times,yield_plot,occupancy,organise}
                     ...
    poretools: error: too few arguments

(Ignore the error, we're expecting it because we have not given it any arguments!)

To install the rest of the software, we will use Linux brew: https://github.com/Linuxbrew/brew
::
    sudo mkdir /home/linuxbrew
    sudo chown $USER:$USER /home/linuxbrew
    git clone https://github.com/Linuxbrew/brew.git /home/linuxbrew/.linuxbrew
    export PATH=/home/linuxbrew/.linuxbrew/bin:$PATH
    brew tap homebrew/science
    
Now install `canu <http://canu.readthedocs.io/en/stable/tutorial.html>`__, `samtools <https://github.com/samtools/samtools/>`__, `bwa mem <http://bio-bwa.sourceforge.net/>`__, `nanopolish <https://github.com/jts/nanopolish>`__:
::
    brew install jdk canu bwa samtools
    
We will install the *new* (thanks `@jts <https://github.com/jts>`__!) R9 version of `nanopolish <https://github.com/jts/nanopolish>`__:
::
    git clone --recursive https://github.com/jts/nanopolish.git
    cd nanopolish/
    make -j 4
    ./nanopolish_test

You should see:
::
    ===============================================================================
    All tests passed (421 assertions in 4 test cases)

Get Oxford Nanopore MinION data
===============================

Last week we collected about 46k reads from three flowcells. Download a subset of these reads:
::
    wget https://s3.amazonaws.com/ngs2016/ectocooler_onp_subset.zip
    mkdir ectocooler_subset/
    unzip ectocooler_onp_subset.zip ectocooler_subset/
    ls ectocooler_subset/
    
You should see a bunch of .fast5 files.

Download the fastq:
::
    wget https://s3.amazonaws.com/ngs2016/ectocooler_onp_all.fastq.gz
    gunzip ectocooler_onp_all.fastq.gz

Convert ONP data in .fast5 to .fastq and .fasta
===============================================

As the MinION instrument is collecting raw data, it is uploaded to the Metrichor server which runs the basecalling software. Reads are then downloaded as .fast5 files. Let's assess the run.
::
    cd
    directory="ectocooler_subset/"
    poretools stats $directory

Here are the 2D reads:
::
    poretools stats --type 2D $directory

How many reads are there? How many 2D? What is the longest read? 

This is only a subset of the reads from the whole run. All of the .fast5 files from the three flowcells we used was 30GB! (`This is a report I generated last week. <https://github.com/ljcohen/dib_ONP_MinION/blob/master/Ectocooler/Ectocooler_read_stats_all3runs.ipynb>`__)

Convert your .fast5 to .fastq and/or .fasta files:
::
    cd ~/
    poretools fastq $directory > ectocooler_subset.fastq
    poretools fasta $directory > ectocooler_subset.fasta

Take a look at a few reads with web blastn. Try to identify what species or closest taxa these data came from. What do you come up with?

Assemble the data
==================

We will use the canu assembler on the full dataset:
::
    canu \
        -p ecto -d ectocooler_assembly \
        genomeSize=3.0m \
        -nanopore-raw ectocooler_onp_all.fastq

Or the subset of data:
::
    canu \
        -p ecto_subset -d ectocooler_assembly \
        genomeSize=3.0m \
        -nanopore-raw ectocooler_subset.fastq

Try both! Compare with your neighbor. 

From the output files, you are interested in the ``ecto.contigs.fasta`` (or ``ecto_subset.contigs.fast``) file. How many contigs do you have? How many contigs are you expecting? How many do you have?

Annotate with prokka:
=====================
Run this command to run prokka:
::
    prokka --outdir anno --prefix prokka contigs.fasta

Check the output:
::
    cat ./anno/prokka.txt

How many genes did Prokka find in the contigs?

Does this meet your expectations?

Evaluate the assembly:
======================

Here is the command:
::
    bwa mem -t 4 -x ont2d ecto.contigs.fasta ectocooler_onp_all.fastq | samtools sort > ectocooler_align.sorted.bam

This will give you a ectocooler_align.sorted.bam.bai
::
    samtools index mapped_reads.sorted

Download the resulting ectocooler_align.sorted.bam, ectocooler_align.sorted.bam.bai, ecto.contigs.fasta to your local computer.
::
    scp -i amazon.pem ubuntu@xxx.amazon.com:/home/ubuntu/ectocooler_align.sorted.bam .
    scp -i amazon.pem ubuntu@xxx.amazon.com:/home/ubuntu/ectocooler_align.sorted.bam.bai
    scp -i amazon.pem ubuntu@xxx.amazon.com:/home/ubuntu/ecto.contigs.fasta

Download this closely-related species:
::
    wget https://github.com/ljcohen/dib_ONP_MinION/blob/master/Ectocooler/Tenacibaculum_dicentrarchi_CP013671.fasta

Open all of these in IGV.

What does it look like? What is the coverage like? Can you spot any problems? What is the Oxford Nanopore error profile? Does it do badly in any regions, which ones? Why?

Fix the assembly with nanopolish
================================

Run this command using your reads and your assembly:
::
    

4. Evaluation of the assembly with alignment of reads to the assembled contigs

   * indexing the reference genome - in this case the reference genome is our de novo assembly
   * aligning, converting SAM to BAM, then sorting the BAM file
   * indexing the BAM file

We will first use the screen command so that we can start the program and then walk away. You can close your computer and the program will keep running. Type Ctrl-A-D to detach and then again Ctrl-A-D to return to the screen later. This is a good time to get a cup of coffee or have lunch!
::
    screen


References:
===========

https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies

Acknowledgements
================

This is a modified lesson by `Nick Loman <http://angus.readthedocs.io/en/2015/analyzing_nanopore_data.html>`__ from 2015, contributions by Torsten Seeman, Harriet Alexander, and Lisa Cohen.
