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
echo 'export PATH=/home/linuxbrew/.linuxbrew/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
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
echo 'export PATH=$PWD/prokka/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
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

Install [mummer](http://mummer.sourceforge.net/) and gnuplot v4:

mummer
```
cd
wget https://github.com/mummer4/mummer/releases/download/v3.9.4alpha/mummer-3.9.4alpha.tar.gz
tar xvzf mummer-3.9.4alpha.tar.gz
cd mummer-3.9.4alpha
./configure
make
sudo make install
echo 'export LD_LIBRARY_PATH="/usr/local/lib"' >> ~/.bashrc
source ~/.bashrc
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

Grab lengths from the killifish runs to compare:

```
wget https://raw.githubusercontent.com/ngs-docs/angus/2017/killifish_lengths.txt
```

Start RStudio server:

```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

Run these commands in RStudio:

```
setwd("~/")
lengths <- read.table("lengths.txt")[,1]
hist(lengths, xlim=c(0,30000), breaks=100, col="red")
killifish_lengths <- read.table("killifish_lengths.txt")[,1]
hist(killifish_lengths, xlim=c(0,90000), breaks=1000, col="blue")

```
## Assemble the data

We will use the program canu to assemble the reads. The full data set will take several hours. So, we will only assemble the subset. Which data are better to use, 2D or a mixture of template and compliment? Pick one, assemble, and compare with your neighbor.

```
canu \
-p ecto_subset -d ectocooler_assembly \
genomeSize=3.0m \
-nanopore-raw ectocooler_subset.fastq
```

This will take ~15 min. So, go ahead and grab some coffee!

After the assembly has finished, the output file you are interested is `ecto_subset.contigs.fasta`. Let's copy that file to the home directory:

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

## All-by-all comparisons

Before running `mummerplot`, run this:

```
sudo awk -i inplace '/\$P_FORMAT .=/ { print "#" $0; next } { print }' $(which mummerplot)
```

because [reasons](https://sourceforge.net/p/mummer/mailman/message/34939032/)


Quick assembly vs. assembly of full reads:

```
nucmer --maxmatch -c 100 -p ectocooler ecto.contigs.fasta ecto_subset.contigs.fasta
mummerplot --fat --filter --png --large -p ectocooler ectocooler.delta
```

Grab closely-related reference genome of [Tenacibaculum dicentrarchi](https://www.ncbi.nlm.nih.gov/genome/?term=txid669041[orgn]):

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/483/385/GCF_001483385.1_ASM148338v1/GCF_001483385.1_ASM148338v1_genomic.fna.gz
gunzip GCF_001483385.1_ASM148338v1_genomic.fna.gz
```
Do all-by-all comparison of reference genome vs. ecocooler full reads assembly

```
nucmer --maxmatch -c 100 -p ectocooler_Tdicent ecto.contigs.fasta GCF_001483385.1_ASM148338v1_genomic.fna
mummerplot --fat --filter --png --large -p ectocooler_Tdicent ectocooler_Tdicent.delta
```

Annotate with prokka:
=====================

This week, you have used Torsten's program, [prokka](http://angus.readthedocs.io/en/2016/prokka_genome_annotation.html) to annotate a bacterial genome. We will use this to annotate these new contigs we have assembled.

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

## References:

* https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies
* http://www.nature.com/nmeth/journal/v12/n8/full/nmeth.3444.html
* https://github.com/ljcohen/dib_ONP_MinION
* http://nbviewer.jupyter.org/github/arq5x/poretools/blob/master/poretools/ipynb/test_run_report.ipynb
* http://porecamp.github.io/2015/timetable.html
* http://porecamp.github.io/2016/
* http://porecamp.github.io/texas/

## Acknowledgements

This is a modified lesson by [Nick Loman](http://angus.readthedocs.io/en/2015/analyzing_nanopore_data.html) from 2015, contributions by Torsten Seemann, Harriet Alexander, Mick Watson, Danny Miller, Jon Badalamenti, and Lisa Cohen.

## canu.report stats from Fundulus olivaceus reads

If you're interested in read length distributions from R9.5 data:
```
[ljcohen@dev-intel14 killifish_assembly]$ cat killifish.report
[CORRECTION/READS]
--
-- In gatekeeper store 'correction/killifish.gkpStore':
--   Found 614587 reads.
--   Found 4883910982 bases (4.43 times coverage).
--
--   Read length histogram (one '*' equals 4297.12 reads):
--        0   4999 300799 **********************************************************************
--     5000   9999 145553 *********************************
--    10000  14999  81269 ******************
--    15000  19999  40302 *********
--    20000  24999  20366 ****
--    25000  29999  11161 **
--    30000  34999   6305 *
--    35000  39999   3560
--    40000  44999   2043
--    45000  49999   1217
--    50000  54999    734
--    55000  59999    463
--    60000  64999    297
--    65000  69999    165
--    70000  74999    103
--    75000  79999     61
--    80000  84999     33
--    85000  89999     30
--    90000  94999     15
--    95000  99999     17
--   100000 104999     13
--   105000 109999      5
--   110000 114999      5
--   115000 119999      5
--   120000 124999      8
--   125000 129999      2
--   130000 134999      8
--   135000 139999      6
--   140000 144999      4
--   145000 149999      0
--   150000 154999      2
--   155000 159999      3
--   160000 164999      2
--   165000 169999      2
--   170000 174999      0
--   175000 179999      1
--   180000 184999      0
--   185000 189999      4
--   190000 194999      0
--   195000 199999      2
--   200000 204999      0
--   205000 209999      1
--   210000 214999      1
--   215000 219999      1
--   220000 224999      0
--   225000 229999      0
--   230000 234999      0
--   235000 239999      1
--   240000 244999      1
--   245000 249999      1
--   250000 254999      0
--   255000 259999      1
--   260000 264999      0
--   265000 269999      1
--   270000 274999      0
--   275000 279999      0
--   280000 284999      1
--   285000 289999      1
--   290000 294999      1
--   295000 299999      0
--   300000 304999      0
--   305000 309999      0
--   310000 314999      0
--   315000 319999      0
--   320000 324999      0
--   325000 329999      0
--   330000 334999      1
--   335000 339999      1
--   340000 344999      0
--   345000 349999      0
--   350000 354999      0
--   355000 359999      0
--   360000 364999      0
--   365000 369999      0
--   370000 374999      0
--   375000 379999      2
--   380000 384999      0
--   385000 389999      0
--   390000 394999      0
--   395000 399999      0
--   400000 404999      0
--   405000 409999      0
--   410000 414999      0
--   415000 419999      1
--   420000 424999      1
--   425000 429999      0
--   430000 434999      0
--   435000 439999      0
--   440000 444999      0
--   445000 449999      0
--   450000 454999      0
--   455000 459999      0
--   460000 464999      0
--   465000 469999      0
--   470000 474999      1
--   475000 479999      0
--   480000 484999      0
--   485000 489999      0
--   490000 494999      0
--   495000 499999      0
--   500000 504999      1
--   505000 509999      0
--   510000 514999      0
--   515000 519999      0
--   520000 524999      0
--   525000 529999      0
--   530000 534999      0
--   535000 539999      0
--   540000 544999      0
--   545000 549999      0
--   550000 554999      0
--   555000 559999      0
--   560000 564999      0
--   565000 569999      0
--   570000 574999      0
--   575000 579999      0
--   580000 584999      0
--   585000 589999      0
--   590000 594999      0
--   595000 599999      0
--   600000 604999      0
--   605000 609999      0
--   610000 614999      0
--   615000 619999      0
--   620000 624999      0
--   625000 629999      0
--   630000 634999      0
--   635000 639999      0
--   640000 644999      0
--   645000 649999      0
--   650000 654999      0
--   655000 659999      0
--   660000 664999      0
--   665000 669999      0
--   670000 674999      0
--   675000 679999      0
--   680000 684999      0
--   685000 689999      0
--   690000 694999      0
--   695000 699999      0
--   700000 704999      0
--   705000 709999      0
--   710000 714999      0
--   715000 719999      0
--   720000 724999      0
--   725000 729999      0
--   730000 734999      0
--   735000 739999      0
--   740000 744999      0
--   745000 749999      0
--   750000 754999      0
--   755000 759999      0
--   760000 764999      0
--   765000 769999      0
--   770000 774999      0
--   775000 779999      0
--   780000 784999      0
--   785000 789999      0
--   790000 794999      0
--   795000 799999      0
--   800000 804999      0
--   805000 809999      0
--   810000 814999      0
--   815000 819999      0
--   820000 824999      0
--   825000 829999      0
--   830000 834999      0
--   835000 839999      0
--   840000 844999      0
--   845000 849999      0
--   850000 854999      0
--   855000 859999      0
--   860000 864999      0
--   865000 869999      0
--   870000 874999      0
--   875000 879999      0
--   880000 884999      0
--   885000 889999      0
--   890000 894999      0
--   895000 899999      1
--   900000 904999      1
--   905000 909999      0
--   910000 914999      0
--   915000 919999      0
--   920000 924999      0
--   925000 929999      0
--   930000 934999      0
--   935000 939999      0
--   940000 944999      0
--   945000 949999      0
--   950000 954999      0
--   955000 959999      0
--   960000 964999      0
--   965000 969999      0
--   970000 974999      1

[CORRECTION/MERS]
--
--  16-mers                                                                                           Fraction
--    Occurrences   NumMers                                                                         Unique Total
--       1-     1 513587462 *******************************************************************--> 0.3659 0.1054
--       2-     2 302853120 ********************************************************************** 0.5817 0.2296
--       3-     4 297949496 ********************************************************************   0.7116 0.3419
--       5-     7 166204607 **************************************                                 0.8485 0.5152
--       8-    11  72366333 ****************                                                       0.9314 0.6771
--      12-    16  28375480 ******                                                                 0.9700 0.7901
--      17-    22  11291731 **                                                                     0.9861 0.8576
--      23-    29   4915132 *                                                                      0.9929 0.8966
--      30-    37   2381975                                                                        0.9960 0.9200
--      38-    46   1279754                                                                        0.9975 0.9351
--      47-    56    750151                                                                        0.9983 0.9454
--      57-    67    471014                                                                        0.9988 0.9529
--      68-    79    316003                                                                        0.9992 0.9587
--      80-    92    217116                                                                        0.9994 0.9633
--      93-   106    152316                                                                        0.9995 0.9671
--     107-   121    108860                                                                        0.9996 0.9701
--     122-   137     79676                                                                        0.9997 0.9726
--     138-   154     60884                                                                        0.9998 0.9747
--     155-   172     47146                                                                        0.9998 0.9765
--     173-   191     36929                                                                        0.9998 0.9780
--     192-   211     29814                                                                        0.9999 0.9794
--     212-   232     23866                                                                        0.9999 0.9806
--     233-   254     19555                                                                        0.9999 0.9817
--     255-   277     15830                                                                        0.9999 0.9826
--     278-   301     13153                                                                        0.9999 0.9835
--     302-   326     11139                                                                        0.9999 0.9843
--     327-   352      9446                                                                        0.9999 0.9850
--     353-   379      7738                                                                        0.9999 0.9856
--     380-   407      6591                                                                        1.0000 0.9862
--     408-   436      5681                                                                        1.0000 0.9867
--     437-   466      5128                                                                        1.0000 0.9872
--     467-   497      4612                                                                        1.0000 0.9877
--     498-   529      3954                                                                        1.0000 0.9882
--     530-   562      3413                                                                        1.0000 0.9886
--     563-   596      3177                                                                        1.0000 0.9890
--     597-   631      2742                                                                        1.0000 0.9893
--     632-   667      2515                                                                        1.0000 0.9897
--     668-   704      2268                                                                        1.0000 0.9900
--     705-   742      1991                                                                        1.0000 0.9903
--     743-   781      1873                                                                        1.0000 0.9906
--     782-   821      1654                                                                        1.0000 0.9909
--
--      620808 (max occurrences)
--  4361104715 (total mers, non-unique)
--   890053109 (distinct mers, non-unique)
--   513587462 (unique mers)
```
