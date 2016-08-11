.. highlight:: rst

Aligning (and Mapping) RNA-seq reads
====================================

During this lab, you'll learn how to *align* RNA-seq reads to the genome using `STAR <https://github.com/alexdobin/STAR>`_, and
how to *map* RNA-seq reads directly to the transcriptome using `RapMap <https://github.com/COMBINE-lab/RapMap>`_.

``> ssh -i ~/Downloads/?????.pem ubuntu@XX.XX.XX.XX``

Update the package list
-----------------------

``> sudo apt-get update``

Install some base packages
--------------------------

First, install the "build tools" (compilers etc. that may be needed)

``> sudo apt-get install build-essential``

To make use of linuxbrew, we'll need Ruby and Git:

``> sudo apt-get install git ruby``

Do the necessary Linuxbrew incantations
---------------------------------------

Get linuxbrew with the following command:

``> ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"``

when prompted, hit `RETURN`.  To make the software we'll install via linuxbrew accessable we have 
to place the directory where linuxbrew installs programs into our ``PATH``.  The following sequence of 
commands will do this:

``> echo 'export PATH="/home/ubuntu/.linuxbrew/bin:$PATH"' >>~/.bashrc``

``> echo 'export MANPATH="/home/ubuntu/.linuxbrew/share/man:$MANPATH"' >>~/.bashrc``

``> echo 'export INFOPATH="/home/ubuntu/.linuxbrew/share/info:$INFOPATH"' >>~/.bashrc``

To make our new paths active, we have to ``source`` them:

``> source ~/.bashrc``

The programs we're interested in installing are part of hombrew-science.  We can "tap" the science keg (;P) as follows:

  ``> brew tap homebrew/science``
  
Now, install samtools:

  ``> brew install homebrew/science/samtools``

Mounting the reads
------------------

We have prepared (thanks; `@monsterbashseq! <https://ljcohen.github.io/>`_) an Amazon volume from which you can load the reads directly.  When we created our AWS instance, we attached the volume with the reads to ``/dev/xvdf``.  We have to *mount* this device.  First, we'll create a place to mount it:

``> sudo mkdir -p /mnt/reads``

Now, we actualy mount the device at the mount point:

``> sudo mount /dev/xvdf /mnt/reads``

When this command finishes (should only take a few seconds) we're good to go, but just need to change the permissions on this folder.

``> sudo chown -R ubuntu:ubuntu /mnt/reads``

Now all of the read files should be available in ``/mnt/reads``.  Check this out with:

``> ls -lha /mnt/reads``

You should see something similar to::


  > ls -lha /mnt/reads/
  total 72G
  drwxr-xr-x 3 ubuntu ubuntu 4.0K Aug 10 22:54 .
  drwxr-xr-x 3 root   root   4.0K Aug 11 03:08 ..
  drwx------ 2 ubuntu ubuntu  16K Aug 10 20:56 lost+found
  -rw-rw-r-- 1 ubuntu ubuntu 3.2G Aug  9 15:18 OREf_SAMm_sdE3_ATTCCT_L002_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 3.3G Aug  9 15:19 OREf_SAMm_sdE3_ATTCCT_L002_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.5G Aug  9 15:20 OREf_SAMm_w_GTCCGC_L006_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.5G Aug  9 15:21 OREf_SAMm_w_GTCCGC_L006_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 3.3G Aug  9 15:08 ORE_sdE3_r1_GTGGCC_L004_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 3.3G Aug  9 15:10 ORE_sdE3_r1_GTGGCC_L004_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 3.4G Aug  9 15:11 ORE_sdE3_r2_TGACCA_L005_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 3.5G Aug  9 15:13 ORE_sdE3_r2_TGACCA_L005_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.8G Aug  9 15:14 ORE_w_r1_ATCACG_L001_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.8G Aug  9 15:15 ORE_w_r1_ATCACG_L001_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 3.5G Aug  9 15:16 ORE_w_r2_GTTTCG_L002_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 3.5G Aug  9 15:17 ORE_w_r2_GTTTCG_L002_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.8G Aug  9 15:29 SAMf_OREm_sdE3_TAGCTT_L001_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.8G Aug  9 15:30 SAMf_OREm_sdE3_TAGCTT_L001_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.9G Aug  9 15:31 SAMf_OREm_w_CAGATC_L005_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 3.0G Aug  9 15:32 SAMf_OREm_w_CAGATC_L005_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 3.7G Aug  9 15:22 SAM_sdE3_r1_ATGTCA_L006_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 3.7G Aug  9 15:23 SAM_sdE3_r1_ATGTCA_L006_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.9G Aug  9 15:24 SAM_sdE3_r2_GCCAAT_L007_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.9G Aug  9 15:25 SAM_sdE3_r2_GCCAAT_L007_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.5G Aug  9 15:26 SAM_w_r1_ACTTGA_L003_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.6G Aug  9 15:27 SAM_w_r1_ACTTGA_L003_R2_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.6G Aug  9 15:28 SAM_w_r2_GAGTGG_L004_R1_001.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 2.6G Aug  9 15:29 SAM_w_r2_GAGTGG_L004_R2_001.fastq.gz


Obtaining the refernece data
----------------------------

We'll be aligning our reads to the Drosophila genome with STAR, and against the Drosophila transcriptome with RapMap.  So, we'll need these files (we'll also need the GTF file corresponding to the genome, as STAR uses this to index known splice sites).

Grab the genome:

``> wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/fasta/dmel-all-chromosome-r6.12.fasta.gz``

and the annotation for the genome

``> wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/gtf/dmel-all-r6.12.gtf.gz``

and the transcriptome

``> wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/fasta/dmel-all-transcript-r6.12.fasta.gz``

We'll put all of these in a folder called ``ref``, and, since they're fairly small, we'll unzip them all::

  > mkdir ref
  > mv *.gz ref
  > cd ref
  > gunzip *.gz
  > cd ..
  
Great; now, we're ready to grab our aligner and align some reads!

Using STAR
--------------

""""""""""""""
Obtaining STAR
""""""""""""""

We'll grab what was, at the time this tutorial was created, the latest version of `STAR <https://github.com/alexdobin/STAR>`_ (v.2.5.2a).  Alex Dobin, the author and maintainer of STAR updates the tool fairly regularly, so you'll ususally want to check for a recent version and the documented changes before you start a new batch of analyses::

  > wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.5.2a.tar.gz
  > tar xzvf 2.5.2a.tar.gz

"""""""""""""""""""""
Building STAR's index
"""""""""""""""""""""

In order to align reads efficiently, STAR has to build an index (in this case, a suffix array), over the genome.  First, we'll create the directory where the index will live:

``> mkdir star_index``

Now, we have to tell STAR to build the index using our reference genome and the GTF annotation.  That command looks like::

  > ~/STAR-2.5.2a/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate \
        --genomeDir star_index --genomeFastaFiles ref/dmel-all-chromosome-r6.12.fasta \
        --sjdbGTFfile ref/dmel-all-r6.12.gtf --sjdbOverhang 99

Here, we're telling STAR that it can use up to 8 threads, and it should build the index on the genome and using the annotation we provide.  The ``sjdbOverhang`` parameter is helpful in setting some internal options, and is recommended to be set as read_lenght - 1.  Once you execute this command, ``STAR`` should run for ~3 minutes before finishing and placing the index in the requested directory.  We can check the contents of the index file:

``> ls -lha star_index``

And we'll see a bunch of files related to the index built by STAR.  The total size of this folder should be ~3.3G --- quite a bit larger than the input reference genome (140M).  This is one of the trade-offs that STAR makes; to provide very fast alignment speeds (and STAR is *very* fast), it uses a large amount of memory.  When aligning to e.g. the human genome, you should be prepared to have 20-30G of RAM available for STAR.

""""""""""""""""""""""""""""
Aligning the reads with STAR
""""""""""""""""""""""""""""

STAR has a *dizzying* array of options. You can find most of them explained in detail in the `STAR manual <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>`_.  For the purposes of this lesson, we'll keep things simple and I'll explain the main options we're setting here.  First, let's create the output directory where our alignments will live:

``> mkdir alignments``

Now, we'll run STAR to align the reads using the following command::

  > /usr/bin/time ~/STAR-2.5.2a/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir star_index \
      --readFilesIn /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R1_001.fastq.gz /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R2_001.fastq.gz \
      --readFilesCommand gunzip -c --outFileNamePrefix alignments/ORE_sdE3_r1_GTGGCC_L004 --outSAMtype BAM Unsorted

* **--runThreadN** tells STAR to use the specified number of threads
* **--genomeDir** tells STAR where to look for the index
* **--readFilesIn** tells STAR where to find the files it will be aligning.  When aligning paired-end reads, the left and right reads are separated by a space.  If there are multiple files for the left and right reads, these are separated by commas.  It is important, if you have multiple left and right files, that they are given in the *same order* in their respective lists.
* **--readFilesCommand** tells star what command it should use to coerce the input into standard FASTA/FASTQ format.  Here, since our reads are gzipped, we tell STAR to use ``gzip -c`` to produce a ``FASTQ`` format file from the input ``fastq.gz`` format files.
* **--outFileNamePrefix** tells STAR how it should name its output files.  There are defaults, but here we override them with the name of the library we're aligning
* **--outSAMtype** tells STAR the format in which the output should be written.  Here, we're telling STAR that the output should be in BAM format (binary and compressed), and that it's OK for the alignments to be unsorted by position / name / etc.

This command will take a little while to run.  While STAR is doing it's thing, we can monitor it's progress with this nifty little command:

``> tail -f alignments/ORE_sdE3_r1_GTGGCC_L004Log.progress.out``

The ``tail -f`` command will *follow* the file, and will write the end (tail) of the  file to the console when it's updated.  At this point, while we wait, it would be an ideal time to discuss what STAR is doing, or to answer any questions you might have.



commands::
  
  > wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.5.2a.tar.gz
  > tar xzvf 2.5.2a.tar.gz
  > mkdir star_index
  > ~/STAR-2.5.2a/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate \
        --genomeDir star_index --genomeFastaFiles ref/dmel-all-chromosome-r6.12.fasta \
        --sjdbGTFfile ref/dmel-all-r6.12.gtf --sjdbOverhang 99
  > mkdir alignments
  > /usr/bin/time ~/STAR-2.5.2a/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir star_index \
        --readFilesIn /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R1_001.fastq.gz /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R2_001.fastq.gz \
        --readFilesCommand gunzip -c --outFileNamePrefix alignments/ORE_sdE3_r1_GTGGCC_L004 --outSAMtype BAM Unsorted
  > wget --no-check-certificate https://github.com/COMBINE-lab/RapMap/releases/download/v0.3.0/RapMap-v0.3.0_linux_x86-64.tar.gz
  > tar xzvf RapMap-v0.3.0_linux_x86-64.tar.gz
  > ~/RapMap-v0.3.0_CentOS5/bin/rapmap quasiindex -t ref/dmel-all-transcript-r6.12.fasta -i rapmap_index
  > mkdir mappings
  > ~/RapMap-v0.3.0_CentOS5/bin/rapmap quasimap -i rapmap_index -t 8 -1 <(gunzip -c /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R1_001.fastq.gz) -2 <(gunzip -c /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R2_001.fastq.gz) | samtools -Sb -@4 - > mappings/mapped_reads.bam
  
Now, we can get linuxbrew with the following command:

``> ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"``

when prompted, hit `RETURN`.  To make the software we'll install via linuxbrew accessable we have 
to place the directory where linuxbrew installs programs into our "path".  The following sequence of 
commands will do this:

``> echo 'export PATH="/home/ubuntu/.linuxbrew/bin:$PATH"' >>~/.bash_profile``

``> echo 'export MANPATH="/home/ubuntu/.linuxbrew/share/man:$MANPATH"' >>~/.bash_profile``

``> echo 'export INFOPATH="/home/ubuntu/.linuxbrew/share/info:$INFOPATH"' >>~/.bash_profile``

The programs we're interested in installing are part of hombrew-science.  We can "tap" the science keg (;P) as follows:

  ``> brew tap homebrew/science``
  ``> brew install homebrew/science/samtools``

Now, we can install the STAR aligner like so:

``> brew install homebrew/science/rna-star``

Since the latest (pre-release) salmon is not yet a binary available in linuxbrew, we'll grab a pre-compiled binary directly.
We can download it using `wget` like so:

``> wget --no-check-certificate 'https://drive.google.com/uc?export=download&id=0B3iS9-xjPftjaFQwYUlvQnN0UFU' -O Salmon-v0.7.0.tgz``

and we can untar and unzip the resulting file with the following command:

``> tar xzf Salmon-v0.7.0.tgz``

Finally, so that we can simply type `salmon` to execute salmon, we'll add the appropriate directory to our path variable again.

``> echo 'export PATH="SalmonBeta-0.7.0-pre-july27_CentOS5/bin:$PATH"' >>~/.bash_profile``

The path will be automatically set when you login, but we want the changes to take effect now, so we must "source" the 
file containing all of the commands that we created:

``> source ~/.bash_profile``
