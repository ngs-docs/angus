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

Start a blank Amazon instance (m3.xlarge) and `log in <http://angus.readthedocs.io/en/2016/amazon/index.html>`__. VERY IMPORTANT: Make sure to select "us-east-1d" when you create your instance.

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
    
We will install the *new* (thanks @jts!) R9 version of `nanopolish <https://github.com/jts/nanopolish>`__:
::
    wget https://github.com/jts/nanopolish/archive/v0.5.0.tar.gz
    tar -xzf v0.5.0.tar.gz
    cd nanopolish-0.5.0/
    make -j 4

Get Oxford Nanopore MinION data
===============================

Last week we got about 46k reads. They are saved in an AWS snapshot: `snap-c4d9a35c`. First, we will create a new volume of this image then attach and mount it to the instance you just created.

1. In your web browser, go to the 'Volumes' tab in the AWS EC2 web page. Select 'Create Volume'.
2. Enter the information. VERY IMPORTANT: Make sure to select the Availability Zone 'us-east-1d'. (Notice a pattern?)
3. Attach volume with your instance ID. (Your instance should say "Running".) Note the 'Device': /dev/sdf
4. Go to your terminal
5. type these commands: (Here is a tricky bit, that the device above is `/dev/sdf` shoudl actually be `/dev/xvdf` below. If your drive letter is slightly different, `a` or `d` etc, change only the last letter below.)
::
        df -h
        sudo mount /dev/xvdf /mnt
        sudo chown -R ubuntu:ubuntu /mnt
        df -h
        cd /mnt/
        ls

You should see a directory called `ectocooler/`. This directory contains >46,000 reads. DO NOT use `ls` in this directory, because there are SO many files!

Now we will work with these files.

Convert ONP data in .fast5 to .fastq and .fasta
===============================================

As the MinION instrument is collecting data, it is uploaded to the Metrichor server which runs the basecalling software. Reads are then downloaded as .fast5 files. Let's assess the run.
::
    directory="/mnt/ectocooler"
    poretools stats -q $directory
    poretools stats -q --type 2D $directory

How many reads are there? How many 2D? What is the longest read? Write these down or save this information.

A directory of ~30 GB of .fast5 files is useless! Convert these to .fastq and/or .fasta files:
::
    poretools fastq $directory > ectocooler.fastq
    poretools fasta $directory > ectocooler.fasta

Take a look at a few reads with web blastn. Try to identify what species or closest taxa these data came from. What do you come up with?

Find the closest complete genome and download. (Need more instructions here.)

Assemble the data
==================

We will use canu.
::
    canu \
        -p ecto -d ectocooler_assembly \
        genomeSize=3.0m \
        -nanopore-raw ectocooler.fastq

This will give you a series of files output. You are interested in the ``ecto.contigs.fasta`` file. How many contigs do you have? How many contigs are you expecting? How many do you have? Is this a good assembly?

Where are the discontinuities? (Hint: find and look at the diagonal plot.)

https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies

Fix the assembly with nanopolish
================================

Run this command using your reads and your assembly:
::
    make -f /home/ubuntu/.linuxbrew/Cellar/nanopolish/0.4.0/scripts/consensus.make READS=/mnt/Ectocooler/Ectocooler_all.fasta ASSEMBLY=/mnt/Ectocooler/Ectocooler_assembly/canu_3m_er08/ecto.contigs.fasta

4. Evaluation of the assembly with alignment of reads to the assembled contigs

   * indexing the reference genome - in this case the reference genome is our de novo assembly
   * aligning, converting SAM to BAM, then sorting the BAM file
   * indexing the BAM file

We will first use the screen command so that we can start the program and then walk away. You can close your computer and the program will keep running. Type Ctrl-A-D to detach and then again Ctrl-A-D to return to the screen later. This is a good time to get a cup of coffee or have lunch!
::
    screen

Here is the command:
::
    /home/ubuntu/bwa-0.7.15/bwa mem -t 4 -x ont2d ecto.contigs.fasta ../Ectocooler/Ectocooler_all.fastq | /home/ubuntu/samtools-1.3.1/samtools sort > ectocooler_align.sorted.bam

This will give you a mapped_reads.sorted.bam.bai
::
    samtools index mapped_reads.sorted

Download the resulting mapped_reads.sorted.bam, mapped_reads.sorted.bam.bai and nanopore-ecoli-sc/scaffolds.fasta files and open in IGV.

What does it look like? What's the coverage like? Can you spot any problems? What is the Oxford Nanopore error profile? Does it do badly in any regions, which ones? Why?

Acknowledgements
================

This is a modified lesson by (http://angus.readthedocs.io/en/2015/analyzing_nanopore_data.html)[Nick Loman] from 2015, contributions by Torsten Seeman, Harriet Alexander, and Lisa Cohen.
