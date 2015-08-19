========================================
Assembling E. coli sequences with SPAdes
========================================

The goal of this tutorial is to show you the basics of assembly using
`the SPAdes assembler <http://bioinf.spbau.ru/spades>`__.

We'll be using data from `Efficient de novo assembly of single-cell
bacterial genomes from short-read data sets, Chitsaz et al., 2011
<http://www.ncbi.nlm.nih.gov/pubmed/21926975>`__.

Booting an Amazon AMI
~~~~~~~~~~~~~~~~~~~~~

Start up an Amazon computer (m3.large or m3.xlarge) running
Ubuntu 14.04, as in :doc:`amazon/index`, and log in.

Logging in
==========

Log in and type::

   sudo apt-get update && \
   sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
              default-jre pkg-config libncurses5-dev r-base-core \
              r-cran-gplots python-matplotlib sysstat python-virtualenv \
              python-setuptools cmake

to update the computer with all the bundled software you'll need.

At this time, you might also make /mnt writeable::

   sudo chmod a+rwxt /mnt

Packages to install
===================

Install `khmer <http://khmer.readthedocs.org/>`__::

   cd
   python -m virtualenv env
   source env/bin/activate
   pip install -U setuptools
   pip install khmer==1.4.1

and download and compile the SPAdes assembler::

   cd
   curl -O http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0.tar.gz
   tar xvf SPAdes-3.5.0.tar.gz
   cd SPAdes-3.5.0
   ./spades_compile.sh
   export PATH="$PATH:$(pwd)/bin"

as well as `Quast <http://quast.bioinf.spbau.ru/manual.html>`__,
software for evaluating the assembly against the known reference: ::

   cd
   curl -L http://sourceforge.net/projects/quast/files/quast-3.0.tar.gz/download > quast-3.0.tar.gz
   tar xvf quast-3.0.tar.gz

Getting the data
================

Now, let's create a working directory::

   cd /mnt
   mkdir assembly
   cd assembly

Download some E. coli data.  This data set
(ecoli_ref-5m-trim.fastq.gz) is the trimmed data from the Chitsaz
paper, E. coli reference sequencing. ::

   curl -O https://s3.amazonaws.com/public.ged.msu.edu/ecoli_ref-5m-trim.fastq.gz

Now, pull out the paired reads::

   extract-paired-reads.py ecoli_ref-5m-trim.fastq.gz
   mv ecoli_ref-5m-trim.fastq.gz.se ecoli_ref-5m-trim.se.fq
   mv ecoli_ref-5m-trim.fastq.gz.pe ecoli_ref-5m-trim.pe.fq

Running an assembly
===================

Now, let's run an assembly::

   spades.py --12 ecoli_ref-5m-trim.pe.fq -s ecoli_ref-5m-trim.se.fq -o spades.d

This will take about 15 minutes; it should end with::


   * Corrected reads are in /mnt/assembly/spades.d/corrected/
   * Assembled contigs are in /mnt/assembly/spades.d/contigs.fasta (contigs.fastg)
   * Assembled scaffolds are in /mnt/assembly/spades.d/scaffolds.fasta (scaffolds.fastg)

Looking at the assembly
=======================

Run QUAST::

   ~/quast-3.0/quast.py spades.d/scaffolds.fasta -o report

and then look at the report::

   less report/report.txt

You should see::

   All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

   Assembly                   scaffolds
   # contigs (>= 0 bp)        160      
   # contigs (>= 1000 bp)     84       
   Total length (>= 0 bp)     4571783  
   Total length (>= 1000 bp)  4551354  
   # contigs                  93       
   Largest contig             264754   
   Total length               4557807  
   GC (%)                     50.75    
   N50                        132618   
   N75                        64692    
   L50                        12       
   L75                        24       
   # N's per 100 kbp          0.00     

Comparing and evaluating assemblies - QUAST
===========================================

Download the true reference genome::

   cd /mnt/assembly
   curl -O https://s3.amazonaws.com/public.ged.msu.edu/ecoliMG1655.fa.gz
   gunzip ecoliMG1655.fa.gz

and run QUAST again::

   ~/quast-3.0/quast.py -R ecoliMG1655.fa spades.d/scaffolds.fasta -o report

Note that here we're looking at *all* the assemblies we've generated.

Now look at the results::

   less report/report.txt

and now we have a lot more information!

A second assembler - MEGAHIT
============================

Let's try out the `MEGAHIT assembler
<http://www.ncbi.nlm.nih.gov/pubmed/25609793>`__.  MEGAHIT is
primarily intended for metagenomes but works well on microbial genomes
in general.

The MEGAHIT source code is on GitHub, here:
https://github.com/voutcn/megahit.  Let's go grab it and build it!
::

   cd
   git clone https://github.com/voutcn/megahit.git
   cd megahit
   make

Now, let's go run an assembly -- ::

   cd /mnt/assembly
   ~/megahit/megahit --12 *.pe.fq -r *.se.fq

This will take about a minute, and the output will be placed in
``megahit_out/final.contigs.fa``.  Let's evaluate it against the SPAdes
assembly with QUAST::

   cp spades.d/scaffolds.fasta spades-assembly.fa
   cp megahit_out/final.contigs.fa megahit-assembly.fa
   ~/quast-3.0/quast.py -R ecoliMG1655.fa spades-assembly.fa \
            megahit-assembly.fa -o report

Let's look at the report! ::

   less report/report.txt

Reference-free comparison
=========================

Above, we've been using the genome reference to do assembly
comparisons -- but often you don't have one. What do you do to
evaluate and compare assemblies without a reference?

One interesting trick is to just run QUAST with one assembly as a reference,
and the other N assemblies against it.  My only suggestion is to first
eliminate short, fragmented contigs from the assembly you're going to use
as a reference.

Let's try that, using ``extract-long-sequences.py`` from `khmer
<http://khmer.readthedocs.org>`__::

   extract-long-sequences.py -l 1000 spades-assembly.fa > spades-long.fa

and then re-run QUAST and put the output in ``report-noref/report.txt``::

   ~/quast-3.0/quast.py -R spades-long.fa spades-assembly.fa \
            megahit-assembly.fa -o report-noref

When you look at the report, ::

   less report-noref/report.txt

take particular note of the following -- ::

   Assembly                     spades-assembly  megahit-assembly
   ...
   Misassembled contigs length  0                814643          
   # local misassemblies        0                9               
   # unaligned contigs          9 + 0 part       7 + 14 part     
   Unaligned length             6453             7404            
   Genome fraction (%)          100.000          99.833          

Challenge exercise
==================

Take your assembled genome, and:

* Install BLAST;
* Grab a FASTA sequence from NCBI for an E. coli protein (e.g. `CRP <http://athyra.idyll.org/~t/crp.fa>`__);
* Save it to a file;
* TBLASTN that protein against your newly assembled genome.

See :doc:`running-command-line-blast` for the basics.
Hint -- you'll need to format your assembly as a BLAST database.
