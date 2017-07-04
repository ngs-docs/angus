# Assessing and Assembling Nanopore data

Last year (2016) in Woods Hole, MA we used our [lab's](http://ivory.idyll.org/lab/) [MinION](https://www.nanoporetech.com/) to sequence a new bacterial species isolated by [Rebecca Mickol](https://news.uark.edu/articles/27669/earth-organisms-survive-under-low-pressure-martian-condition) in the [Microbial Diversity Course at the Marine Biological Lab](http://www.mbl.edu/microbialdiversity/).

If you're interested, you can read a [blog post](https://monsterbashseq.wordpress.com/2016/08/) about it.

The goals of this tutorial are to:

*  Assess an Oxford Nanopore Technologies (ONT) sequencing run on the MinION
*  Create an assembly from raw fastq files
*  Evaluate the assembly

## Start a Jetstream instance and install software:

[Start a blank Jetstream instance](https://angus.readthedocs.io/en/2017/jetstream/boot.html) (m1.xlarge) with CPU: 24, Mem: 60, Disk 60 GB RAM:

Copy/paste to update and install software on your new instance:

```
    sudo apt-get update && \
    sudo apt-get -y install build-essential ruby screen git curl gcc make g++ python-dev unzip \
        default-jre pkg-config libncurses5-dev r-base-core \
        r-cran-gplots python-matplotlib sysstat python-virtualenv \
        python-setuptools cmake cython libhdf5-serial-dev \
        python-numpy python-scipy python-pandas python-pandas-lib \
        python-biopython parallel python-h5py python-tornado \
        bioperl libxml-simple-perl default-jre gdebi-core r-base gnuplot
```

To install some of the software, we will use [Linux brew](https://github.com/Linuxbrew/brew):

```
    sudo mkdir /home/linuxbrew
    sudo chown $USER:$USER /home/linuxbrew
    git clone https://github.com/Linuxbrew/brew.git /home/linuxbrew/.linuxbrew
    export PATH=/home/linuxbrew/.linuxbrew/bin:$PATH
    brew tap homebrew/science
```

Now install [canu](http://canu.readthedocs.io/en/stable/tutorial.html), [samtools](https://github.com/samtools/samtools/), [bwa mem](http://bio-bwa.sourceforge.net/):

```
    brew install jdk canu bwa samtools
```

Install poretools:

```
sudo pip install poretools
```

Install prokka:

```
    git clone https://github.com/tseemann/prokka.git
    export PATH=$PWD/prokka/bin:$PATH
    prokka --setupdb
    prokka --version
```

Install assembly-stats:

```
   git clone https://github.com/sanger-pathogens/assembly-stats.git
   cd assembly-stats/
   mkdir build
   cd build
   cmake ..
   make
   make test
   sudo make install
```

Install RStudio:

```
    wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
    sudo gdebi -n rstudio-server-1.0.143-amd64.deb
```
Change your password for RStudio:

```
    sudo passwd <account name>
```

Install [miniasm and minimap](https://github.com/lh3/miniasm):

```
    git clone https://github.com/lh3/minimap && (cd minimap && make)
    git clone https://github.com/lh3/miniasm && (cd miniasm && make)
```

Install mummer:

```
    wget https://github.com/mummer4/mummer/releases/download/v3.9.4alpha/mummer-3.9.4alpha.tar.gz
    tar xvzf mummer-3.9.4alpha.tar.gz
    cd mummer-3.9.4alpha
    ./configure
    make
    sudo make install
```

## Get Oxford Nanopore MinION data and convert it

Our data were collected from three flowcells in 2016. Download a subset of these reads:

```
    cd
    wget https://s3.amazonaws.com/ngs2016/ectocooler_subset.zip
    unzip ectocooler_subset.zip 
    ls ectocooler_subset/
```

You should see a bunch of .fast5 files.
  
This is only a subset of the reads from the whole run. ([Click here for stats from the full data set.](https://github.com/ljcohen/dib_ONP_MinION/blob/master/Ectocooler/Ectocooler_read_stats_all3runs.ipynb))

The MinION instrument collects raw data in .fast5 format. The local basecalling software, [Albacore (sorry, link requires ONT MAP login access)](https://community.nanoporetech.com/downloads), converts .fast5 files into .fastq or .fasta files. Poretools is another method for converting .fast5 files into .fastq files. 

Convert your .fast5 to .fastq and .fasta files:

```
    cd
    directory="ectocooler_subset/"
    poretools fastq $directory > ectocooler_subset.fastq
    poretools fasta $directory > ectocooler_subset.fasta
```

Take a look at the reads:

```
    head ectocooler_subset.fastq
```

Copy a few reads and use the [web blastn](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) to try to identify what species or closest taxa these data came from. What do you come up with?

Download the full dataset, fastq and fasta files:

```
    cd
    wget https://s3.amazonaws.com/ngs2016/ectocooler_all_2D.fastq
    wget https://s3.amazonaws.com/ngs2016/ectocooler_all_2D.fasta
```

## Assess the Data

```
    assembly-stats ectocooler_subset.fastq
```

1. How many reads are there total? 
2. What is the longest read?

Assess the full data set:

```
    assembly-stats ectocooler_all_2D.fastq
```

Run this to make a file with all read lengths :
```
    cat ectocooler_all_2D.fastq | paste - - - - | awk -F"\t" '{print length($2)}' > lengths.txt
```

Start RStudio server:
```
    echo My RStudio Web server is running at: http://$(hostname):8787/
```

## Assemble the data

We will use the program canu to assemble the reads. The full data set will take several hours. So, we will only assemble the subset. Which data are better to use, 2D or a mixture of template and compliment? Pick one, assemble, and compare with your neighbor.

```
    canu \
    -p ecto_subset -d ectocooler_assembly \
    genomeSize=3.0m \
    -nanopore-raw ectocooler_subset.fastq
```

From the output files, you are interested in the ``ecto_subset.contigs.fasta`` file. Let's copy that file to the home directory:

```
    cd
    cp ectocooler_assembly/ecto_subset.contigs.fasta .
```

Assess the assembly:
```
    assembly-stats ecto_subset.contigs.fasta
```

How many contigs do you have? 

Download the pre-assembled contigs from the full data set:

```
    wget https://raw.githubusercontent.com/ljcohen/dib_ONP_MinION/master/Ectocooler/ecto.contigs.fasta
    assembly-stats ecto.contigs.fasta
```

Compare this with your assembly. How are they different?

## All-by-all comparison

```
    /home/ljcohen/mummer-3.9.4alpha/nucmer -maxmatch -c 100 -p ecotcooler ecto.contigs.fasta ecto_subset.contigs.fasta
    /home/ljcohen/mummer-3.9.4alpha/mummerplot -fat -filter -png -large -p ectocooler ectocooler.delta
```

Edit nucmer.gp before running gnuplot
```
    gnuplot ectocooler.gp #edit nucmer.gp before running gnuplot
```

Annotate with prokka:
=====================

Yesterday, you used Torsten's program, `prokka <http://angus.readthedocs.io/en/2016/prokka_genome_annotation.html>`__ to annotate a bacterial genome. We will use this to annotate these new contigs we have assembled.
::
    prokka --outdir anno_subset --prefix ecto_subset_prokka ecto_subset.contigs.fasta

Check the output:
::
    cat ./anno_subset/ecto_subset_prokka.txt

1. How many genes did Prokka find in the contigs?
2. Does this meet your expectations?

Use this command to run prokka on the contigs assembled with the full data set:
::
    prokka --outdir anno_full --prefix ecto_full_prokka ecto.contigs.fasta

Check the output:
::
    cat ./anno_full/ecto_full_prokka.txt

Evaluate the assembly:
======================

Align the reads to the assembled subset of contigs. (Or use the contigs assembled from full data set. Pick one and compare with your neighbor!)

* index the reference genome - in this case the reference genome is our de novo assembly
* align, converting SAM to BAM, then sorting the BAM file
* index the BAM file
   
Index:
::
    bwa index ecto_subset.contigs.fasta

Align
::
    bwa mem -x ont2d -t 8 ecto_subset.contigs.fasta ectocooler_subset_2D.fasta | samtools sort > ecto_subset.sorted.bam
    
This will give you an indexed bam file:
::
    samtools index ecto_subset.sorted.bam

Download the resulting ectocooler_align.sorted.bam, ectocooler_align.sorted.bam.bai, ecto.contigs.fasta to your local computer.
::
    scp -i amazon.pem ubuntu@xxx.amazon.com:/home/ubuntu/ecto_subset.sorted.bam .
    scp -i amazon.pem ubuntu@xxx.amazon.com:/home/ubuntu/ecto_subset.sorted.bam.bai .
    scp -i amazon.pem ubuntu@xxx.amazon.com:/home/ubuntu/ecto_subset.contigs.fasta .

In IGV, open ecto_subset.contigs.fasta as "Genome" and ecto_subset.sorted.bam.

1. What does the alignment look like? 
2. What is the coverage? 
3. Can you spot any problems? 
4. What is the Oxford Nanopore error profile? 
5. Does it do badly in any regions, which ones? Why?


Can fix the assembly using nanopolish


The program `nanopolish <https://github.com/jts/nanopolish>`__ will align your reads to the assembly and compute a consensus. This will take some time.

References:
===========

* https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies
* http://www.nature.com/nmeth/journal/v12/n8/full/nmeth.3444.html
* https://github.com/ljcohen/dib_ONP_MinION
* http://nbviewer.jupyter.org/github/arq5x/poretools/blob/master/poretools/ipynb/test_run_report.ipynb
* http://porecamp.github.io/2015/timetable.html
* http://porecamp.github.io/2016/
* http://porecamp.github.io/texas/

Acknowledgements
================

This is a modified lesson by `Nick Loman <http://angus.readthedocs.io/en/2015/analyzing_nanopore_data.html>`__ from 2015, contributions by Torsten Seemann, Harriet Alexander, and Lisa Cohen.
