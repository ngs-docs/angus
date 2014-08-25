========================================
Assembling E. coli sequences with Velvet
========================================

The goal of this tutorial is to show you the basics of assembly using
`the Velvet assembler
<http://en.wikipedia.org/wiki/Velvet_assembler>`__.

We'll be using data from `Efficient de novo assembly of single-cell
bacterial genomes from short-read data sets, Chitsaz et al., 2011
<http://www.ncbi.nlm.nih.gov/pubmed/21926975>`__.

Booting an Amazon AMI
~~~~~~~~~~~~~~~~~~~~~

Start up an Amazon computer (m1.large or m1.xlarge) using AMI
ami-7607d01e (see :doc:`amazon/start-up-an-ec2-instance` and
:doc:`amazon/starting-up-a-custom-ami`).

Log in `with Windows <amazon/log-in-with-ssh-win.html>`__ or
`from Mac OS X <amazon/log-in-with-ssh-mac.html>`__.

Logging in
==========

Log in and type::

   sudo bash

to change into superuser mode.

Updating the operating system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy and paste the following two commands
::

   apt-get update
   apt-get -y install screen git curl gcc make g++ python-dev unzip \
           default-jre pkg-config libncurses5-dev r-base-core \
           r-cran-gplots python-matplotlib sysstat

to update the computer with all the bundled software you'll need.

Packages to install
===================

Install `khmer <http://khmer.readthedocs.org/en/v1.1/>`__::

   cd /usr/local/share
   git clone https://github.com/ged-lab/khmer.git
   cd khmer
   git checkout v1.1
   make install

and install the Velvet assembler::

   cd /root
   curl -O http://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
   tar xzf velvet_1.2.10.tgz
   cd velvet_1.2.10
   make MAXKMERLENGTH=51
   cp velvet? /usr/local/bin

as well as `Quast <http://quast.bioinf.spbau.ru/manual.html>`__,
software for evaluating the assembly against the known reference: ::

   cd /root
   curl -O -L https://downloads.sourceforge.net/project/quast/quast-2.3.tar.gz
   tar xzf quast-2.3.tar.gz

Getting the data
================

Now, let's create a working directory::

   cd /mnt
   mkdir assembly
   cd assembly

Download some E. coli data.  The first data set
(ecoli_ref-5m-trim.fastq.gz) is the trimmed PE data sets from the
other day (see :doc:`short-read-quality-evaluation`), and the second
data set is a specially processed data set using `digital
normalization <http://ged.msu.edu/papers/2012-diginorm/>`__ that will
assemble quickly. ::

   curl -O https://s3.amazonaws.com/public.ged.msu.edu/ecoli_ref-5m-trim.fastq.gz
   curl -O https://s3.amazonaws.com/public.ged.msu.edu/ecoli-reads-5m-dn-paired.fa.gz

Running an assembly
===================

Now... assemble the small, fast data sets, using the Velvet assembler.  Here
we will set the required parameter k=21::

   velveth ecoli.21 21 -shortPaired -fasta.gz ecoli-reads-5m-dn-paired.fa.gz
   velvetg ecoli.21 -exp_cov auto

Check out the stats for the assembled contigs for a cutoff of 1000::

   python /usr/local/share/khmer/sandbox/assemstats3.py 1000 ecoli.*/contigs.fa

Also try assembling with k=23 and k=25::

   velveth ecoli.23 23 -shortPaired -fasta.gz ecoli-reads-5m-dn-paired.fa.gz
   velvetg ecoli.23 -exp_cov auto

   velveth ecoli.25 25 -shortPaired -fasta.gz ecoli-reads-5m-dn-paired.fa.gz
   velvetg ecoli.25 -exp_cov auto

Now check out the stats for the assembled contigs for a cutoff of 1000::

   python /usr/local/share/khmer/sandbox/assemstats3.py 1000 ecoli.*/contigs.fa

(Also read: `What does k control in de Bruijn graph assemblers? <http://ivory.idyll.org/blog/the-k-parameter.html>`__.)

Comparing and evaluating assemblies - QUAST
===========================================

Download the true reference genome::

   cd /mnt/assembly
   curl -O https://s3.amazonaws.com/public.ged.msu.edu/ecoliMG1655.fa.gz
   gunzip ecoliMG1655.fa.gz

and run QUAST::

   /root/quast-2.3/quast.py -R ecoliMG1655.fa ecoli.*/contigs.fa

Note that here we're looking at *all* the assemblies we've generated.

Now look at the results::

   more quast_results/latest/report.txt

The first bits to look at are Genome fraction (%) and # misassembled contigs,
I think.

Searching assemblies -- BLAST
=============================

Install BLAST::

   cd /root

   curl -O ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.24/blast-2.2.24-x64-linux.tar.gz
   tar xzf blast-2.2.24-x64-linux.tar.gz
   cp blast-2.2.24/bin/* /usr/local/bin
   cp -r blast-2.2.24/data /usr/local/blast-data

Build BLAST databases for the assemblies you've done::

   cd /mnt/assembly

   for i in 21 23 25
   do
      extract-long-sequences.py -o ecoli-$i.fa -l 500 ecoli.$i/contigs.fa
      formatdb -i ecoli-$i.fa -o T -p F
   done

and then let's search for a specific gene -- first, download a file containing
your protein sequence of interest::

   curl -O http://athyra.idyll.org/~t/crp.fa

and now search::

   blastall -i crp.fa -d ecoli-21.fa -p tblastn -b 1 -v 1

Questions and Discussion Points
===============================

Why do we use a lower cutoff of 1kb for the assemstats3 links, above?  Why
not 0?

Followup work
=============

Try running an assembly of the larger read data set::

   velveth ecoli-full.31 31 -short -fastq.gz ecoli_ref-5m-trim.fastq.gz
   velvetg ecoli-full.31 -exp_cov auto

.. @@(You might want to do this in screen.)

