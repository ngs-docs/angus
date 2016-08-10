=======================
Analyzing nanopore data
=======================

Last week in Woods Hole, MA we used our `lab's <http://ivory.idyll.org/lab/>`__ `MinION <https://www.nanoporetech.com/>`__ to sequence a new bacterial species isolated by `Rebecca Mickol <https://news.uark.edu/articles/27669/earth-organisms-survive-under-low-pressure-martian-condition>`__ in the `Microbial Diversity Course at the Marine Biological Lab <http://www.mbl.edu/microbialdiversity/>`__.

The goals of this tutorial are to:
   *  convert oxford nanopore data in .fast5 format to .fastq
   *  assess a nanopore run
   *  assemble data
   *  evaluate the assembly

Starting and AWS instance and installing software:
==================================================

Start a blank Amazon instance (m3.xlarge) and `log in <http://angus.readthedocs.io/en/2016/amazon/index.html>`__.

Copy/paste to update and install software on your new instance:
::
    sudo apt-get update && \
    sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
      default-jre pkg-config libncurses5-dev r-base-core \
      r-cran-gplots python-matplotlib sysstat python-virtualenv \
      python-setuptools cmake cython libhdf5-serial-dev

We will now install several software packages that are specific for analyzing long reads data, as comes from the Oxford Nanopore MinION.

`poretools <http://poretools.readthedocs.io/en/latest/content/installation.html#basic-installation>`__

This requires installing R 3.0:
::
    deb http://www.stats.bris.ac.uk/R/bin/linux/ubuntu trusty/

Now we will install poretools:
::
    git clone https://github.com/arq5x/poretools
    cd poretools
    sudo python setup.py install
    poretools

You should see output like this:
::
    usage: poretools [-h] [-v]
                     {combine,fastq,fasta,stats,hist,events,readstats,tabular,nucdist,metadata,index,qualdist,qualpos,winner,squiggle,times,yield_plot,occupancy,organise}
                     ...
    poretools: error: too few arguments

(Ignore the error, we're expecting it because we have not given it any arguments!)

`canu <http://canu.readthedocs.io/en/stable/tutorial.html>`__

Install:
::
    git clone https://github.com/marbl/canu.git
    cd canu/src
    make -j 4
    canu

You should see output like this:
::
    canu \
       -d <working-directory> \
       -p <file-prefix> \
      [-s specifications] \
      [-correct | -trim | -assemble] \
       errorRate=<fraction-error> \
       genomeSize=<genome-size>\
      [parameters] \
      [-pacbio-raw         <read-file>]
      [-pacbio-corrected   <read-file>]
      [-nanopore-raw       <read-file>]
      [-nanopore-corrected <read-file>]

`samtools <http://www.htslib.org/download/>`__

Install:
::
    wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
    tar -xf samtools-1.3.1.tar.bz2
    cd samtools-1.3.1/
    make
    /home/ubuntu/samtools-1.3.1/samtools/samtools

`bwa mem <http://bio-bwa.sourceforge.net/>`__

Install:
::
    wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2
    tar -xf bwa-0.7.15.tar.bz2 
    cd bwa-0.7.15/
    make
    /home/ubuntu/bwa-0.7.15/bwa mem

`Nanopolish <https://github.com/jts/nanopolish>`__

Has dependencies, `libhdf5 <https://www.hdfgroup.org/HDF5/release/obtain5.html>`__
and gcc-4.8

Install:
::
    git clone --recursive https://github.com/jts/nanopolish.git
    cd nanopolish
    make

Acquiring nanopore data
===============================

Last week we got about 46k reads. You can download them and take a look:
::
    (insert link to data)

Exercise
=========

1.  Evaluation of the run with poretools. How many reads are there? How many 2D? What is the longest read?

Can we identify what species these data came from? Why or why not?

2.  Assembly with canu. What is the N50? Where are the discontiguities (hint: find and look at the diagonal plot).

https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies

3.  Fix the assembly with nanopolish

Edit and run this command using your reads and your assembly:
::
    make -f scripts/consensus.make READS=reads.fa ASSEMBLY=draft.fa

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
