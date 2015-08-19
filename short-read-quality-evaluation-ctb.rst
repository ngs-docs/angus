=============================================================
Evaluating the quality of your short reads, and trimming them
=============================================================

As useful as BLAST is, we really want to get into sequencing data,
right?  One of the first steps you must do with your data is
evaluate its quality and throw away bad sequences.

Before you can do that, though, you need to install a bunch o' software.

Logging in
==========

Log in and type::

   sudo bash

to change into superuser mode.

Packages to install
===================

Install `khmer <http://khmer.readthedocs.org/en/v1.1/>`__::

   cd /usr/local/share
   git clone https://github.com/ged-lab/khmer.git
   cd khmer
   git checkout v1.1
   make install

Install `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`__::

   cd /root
   curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
   unzip Trimmomatic-0.32.zip 
   cp Trimmomatic-0.32/trimmomatic-0.32.jar /usr/local/bin

Install `FastQC <http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc>`__::

   cd /usr/local/share
   curl -O http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip
   unzip fastqc_v0.11.2.zip
   chmod +x FastQC/fastqc

Install libgtextutils and fastx::

   cd /root
   curl -O http://hannonlab.cshl.edu/fastx_toolkit/libgtextutils-0.6.1.tar.bz2
   tar xjf libgtextutils-0.6.1.tar.bz2 
   cd libgtextutils-0.6.1/
   ./configure && make && make install

   cd /root
   curl -O http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit-0.0.13.2.tar.bz2
   tar xjf fastx_toolkit-0.0.13.2.tar.bz2
   cd fastx_toolkit-0.0.13.2/
   ./configure && make && make install

In each of these cases, we're downloading the software -- you can use
google to figure out what each package is and does if we don't discuss
it below.  We're then unpacking it, sometimes compiling it (which we
can discuss later), and then installing it for general use.

Getting some data
=================

Start at your EC2 prompt, then type ::

   cd /mnt

Now, grab the 5m E. coli reads from our data storage (originally from
`Chitsaz et al. <http://www.ncbi.nlm.nih.gov/pubmed/21926975>`__)::

   curl -O https://s3.amazonaws.com/public.ged.msu.edu/ecoli_ref-5m.fastq.gz

You can take a look at the file contents by doing::
 
   gunzip -c ecoli_ref-5m.fastq.gz | less

(use 'q' to quit the viewer).  This is what raw FASTQ looks like!

Note that in this case we've given you the data *interleaved*, which
means that paired ends appear next to each other in the file.  Most of
the time sequencing facilities will give you data that is split out
into s1 and s2 files.  We'll need to split it out into these files for
some of the trimming steps, so let's do that -- ::

   split-paired-reads.py ecoli_ref-5m.fastq.gz
   mv ecoli_ref-5m.fastq.gz.1 ecoli_ref-5m_s1.fq
   mv ecoli_ref-5m.fastq.gz.2 ecoli_ref-5m_s2.fq

This uses the khmer script 'split-paired-reads' `(see documentation)
<http://khmer.readthedocs.org/en/v1.1/scripts.html#scripts-read-handling>`__
to break the reads into left (/1) and right (/2). (This takes a long time!
5m reads is a lot of data...)

We'll also need to get some Illumina adapter information -- here::

   curl -O https://s3.amazonaws.com/public.ged.msu.edu/illuminaClipping.fa

These sequences are (or were) "trade secrets" so it's hard to find 'em.
Don't ask me how I got 'em.

Trimming and quality evaluation of your sequences
=================================================

Start at the EC2 login prompt.  Then, ::

  cd /mnt

Make a directory to store all your trimmed data in, and go there::

  mkdir trim
  cd trim

Now, run `Trimmomatic <http://www.usadellab.org/cms/index.php?page=trimmomatic>`__ to eliminate Illumina adapters from your sequences -- ::

  java -jar /usr/local/bin/trimmomatic-0.32.jar PE ../ecoli_ref-5m_s1.fq ../ecoli_ref-5m_s2.fq s1_pe s1_se s2_pe s2_se ILLUMINACLIP:../illuminaClipping.fa:2:30:10

Next, let's take a look at data quality using `FastQC <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__ ::

   mkdir /root/Dropbox/fastqc
   /usr/local/share/FastQC/fastqc s1_* s2_* --outdir=/root/Dropbox/fastqc

This will dump the FastQC output into your Dropbox folder, under the
folder 'fastqc'.  Go check it out on your local computer in Dropbox --
you're looking for folders named <filename>_fastqc, for example
's1_pe_fastqc'; then double click on 'fastqc_report.html'.
(You can't look at these on the dropbox.com Web site -- it won't interpret
the HTML for you.)

It looks like a lot of bad data is present after base 70, so let's just trim
all the sequences that way.  Before we do that, we want to interleave the
reads again (don't ask) -- ::

   interleave-reads.py s1_pe s2_pe > combined.fq

`(interleave-reads is another khmer scripts) <http://khmer.readthedocs.org/en/v1.1/scripts.html#scripts-read-handling>`__

Now, let's use the FASTX toolkit to trim off bases over 70, and
eliminate low-quality sequences.  We need to do this both for our
combined/paired files and our remaining single-ended files::

   fastx_trimmer -Q33 -l 70 -i combined.fq | fastq_quality_filter -Q33 -q 30 -p 50 > combined-trim.fq

   fastx_trimmer -Q33 -l 70 -i s1_se | fastq_quality_filter -Q33 -q 30 -p 50 > s1_se.filt

Let's take a look at what we have -- ::

   ls -la

You should see::

  drwxr-xr-x 2 root root       4096 2013-04-08 03:33 .
  drwxr-xr-x 4 root root       4096 2013-04-08 03:21 ..
  -rw-r--r-- 1 root root  802243778 2013-04-08 03:33 combined-trim.fq
  -rw-r--r-- 1 root root 1140219324 2013-04-08 03:26 combined.fq
  -rw-r--r-- 1 root root  570109662 2013-04-08 03:23 s1_pe
  -rw-r--r-- 1 root root     407275 2013-04-08 03:23 s1_se
  -rw-r--r-- 1 root root     319878 2013-04-08 03:33 s1_se.filt
  -rw-r--r-- 1 root root  570109662 2013-04-08 03:23 s2_pe
  -rw-r--r-- 1 root root          0 2013-04-08 03:22 s2_se

Let's run FastQC on things again, too::

   mkdir /root/Dropbox/fastqc.filt
   /usr/local/share/FastQC/fastqc combined-trim.fq s1_se.filt --outdir=/root/Dropbox/fastqc.filt

Now go look in your Dropbox folder under 'fastqc.filt', folder
'combined-trim.fq_fastqc' -- looks a lot better, eh?

