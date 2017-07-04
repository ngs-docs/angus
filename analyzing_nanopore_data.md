# Assessing and Assembling Nanopore data

Last year (2016) in Woods Hole, MA we used our [lab's](http://ivory.idyll.org/lab/) [MinION](https://www.nanoporetech.com/) to sequence a new bacterial species isolated by [Rebecca Mickol](https://news.uark.edu/articles/27669/earth-organisms-survive-under-low-pressure-martian-condition) in the [Microbial Diversity Course at the Marine Biological Lab](http://www.mbl.edu/microbialdiversity/).

If you're interested, you can read a [blog post](https://monsterbashseq.wordpress.com/2016/08/13/adventures-with-ont-minion-at-mbls-microbial-diversity-course/) about it.

The goals of this tutorial are to:

*  Assess an Oxford Nanopore Technologies (ONT) sequencing run on the MinION
*  Create an assembly from raw fastq files
*  Evaluate the assembly

## Start a Jetstream instance and install software:

[Start a blank Jetstream instance](https://angus.readthedocs.io/en/2017/jetstream/boot.html) (m1.medium) with CPU: 6, Mem: 16, Disk 60 GB RAM:

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
cd
sudo mkdir /home/linuxbrew
sudo chown $USER:$USER /home/linuxbrew
git clone https://github.com/Linuxbrew/brew.git /home/linuxbrew/.linuxbrew
export PATH=/home/linuxbrew/.linuxbrew/bin:$PATH
brew tap homebrew/science
```

Now install [canu](http://canu.readthedocs.io/en/stable/tutorial.html), [samtools](https://github.com/samtools/samtools/), [bwa](http://bio-bwa.sourceforge.net/):

```
brew install jdk canu bwa samtools
```

Install [poretools](https://poretools.readthedocs.io/en/latest/):

```
sudo pip install poretools
```

Install [prokka](http://www.vicbioinformatics.com/software.prokka.shtml):

```
cd
git clone https://github.com/tseemann/prokka.git
export PATH=$PWD/prokka/bin:$PATH
prokka --setupdb
prokka --version
```

Install [assembly-stats](https://github.com/sanger-pathogens/assembly-stats):

```
cd
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
cd
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi -n rstudio-server-1.0.143-amd64.deb
```
Change your password for RStudio:

```
cd
sudo passwd <account name>
```

Install [mummer](http://mummer.sourceforge.net/):

```
cd
wget https://github.com/mummer4/mummer/releases/download/v3.9.4alpha/mummer-3.9.4alpha.tar.gz
tar xvzf mummer-3.9.4alpha.tar.gz
cd mummer-3.9.4alpha
./configure
make
sudo make install
```

## Get Oxford Nanopore MinION data and convert it

Our data were collected from three R9.4 flowcells in 2016. Download a subset of these reads:

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
How does the full data set compare to the subset?

How does it compare to these results from three R9.5 flowcells of killifish (Fundulus olivaceus) data collected at [Porecamp](http://www.txgen.tamu.edu/porecamp_usa/) in 2017?
```
[ljcohen@globus-00 fastq2]$ /mnt/home/ljcohen/bin/assembly-stats/assembly-stats porecamp_killifish2.fastq
stats for porecamp_killifish2.fastq
sum = 4962626713, n = 740248, ave = 6704.01, largest = 973552
N50 = 12726, n = 117202
N60 = 10357, n = 160433
N70 = 8098, n = 214460
N80 = 5724, n = 286845
N90 = 3229, n = 400661
N100 = 5, n = 740248
N_count = 0
Gaps = 0
```

Run this to make a file with all read lengths:

```
cat ectocooler_all_2D.fastq | paste - - - - | awk -F"\t" '{print length($2)}' > lengths.txt
```

Start RStudio server:

```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

Run these commands in RStudio:

```
lengths <- read.table("lengths.txt")[,1]
hist(lengths, xlim=c(0,30000), breaks=100, col="red")
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
~/mummer-3.9.4alpha/nucmer --maxmatch -c 100 -p ectocooler ecto.contigs.fasta ecto_subset.contigs.fasta
~/mummer-3.9.4alpha/mummerplot --fat --filter --png --large -p ectocooler ectocooler.delta
```

Edit nucmer.gp before running gnuplot to remove the three lines that have the word "mouse".
```
gnuplot ectocooler.gp #edit nucmer.gp before running gnuplot
```

Annotate with prokka:
=====================

This week, you have used Torsten's program, [prokka](http://angus.readthedocs.io/en/2016/prokka_genome_annotation.html)to annotate a bacterial genome. We will use this to annotate these new contigs we have assembled.

```
prokka --outdir anno_subset --prefix ecto_subset_prokka ecto_subset.contigs.fasta
```

Check the output:

```
cat ./anno_subset/ecto_subset_prokka.txt
```

1. How many genes did Prokka find in the contigs?
2. Does this meet your expectations?

Use this command to run prokka on the contigs assembled with the full data set:

```
prokka --outdir anno_full --prefix ecto_full_prokka ecto.contigs.fasta
```

Check the output:

```
cat ./anno_full/ecto_full_prokka.txt
```

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

This is a modified lesson by [Nick Loman](http://angus.readthedocs.io/en/2015/analyzing_nanopore_data.html) from 2015, contributions by Torsten Seemann, Harriet Alexander, Mick Watson, Danny Miller, Jon Badalamenti, and Lisa Cohen.
