.. highlight:: rst


Quantifying transcript expression with Salmon
=============================================
	       
During this lab, you'll learn how to use `salmon <https://github.com/COMBINE-lab/salmon>`_ to rapidly quantify transcript-level expression
from RNA-seq data.

Log into your instance
-----------------------

For this tutorial, we'll use a **c4.2xlarge** instance.  Make sure you create the instance with the
volume containing the reads attached!::

   > ssh -i ~/Downloads/?????.pem ubuntu@XX.XX.XX.XX``

Update the package list
-----------------------

::

   > sudo apt-get update

Install some base packages
--------------------------

First, install the "build tools" (compilers etc. that may be needed)::

  > sudo apt-get install build-essential


Mounting the reads
------------------

We have prepared (thanks; `@monsterbashseq! <https://ljcohen.github.io/>`_) an Amazon volume from which you can load the reads directly.  When we created our AWS instance, we attached the volume with the reads to ``/dev/xvdf``.  We have to *mount* this device.  Since we're using the volume from yesterday, a place for the volume `/mnt/reads` already exists. Here, we just mount the device at the mount point::

  >sudo mount /dev/xvdf /mnt/reads

When this command finishes (should only take a few seconds) we're good to go, but just need to change the permissions on this folder.::

  > sudo chown -R ubuntu:ubuntu /mnt/reads

Now all of the read files should be available in ``/mnt/reads``.  Check this out with::

  > ls -lha /mnt/reads

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

.. note::

   If you're on an instance that already has the reference transcriptome from the mapping lab yesterday, then
   you can skip this step
   
We'll be quantifying against the Drosophila transcriptome, so let's grab that file again:::
  
  > wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/fasta/dmel-all-transcript-r6.12.fasta.gz

We'll put this in a folder called ``ref``, and unzip it there:::

  > mkdir ref
  > mv dmel-all-transcript-r6.12.fasta.gz ref
  > cd ref
  > gunzip dmel-all-transcript-r6.12.fasta.gz 
  > cd ..
  
Great; now, let's run ``salmon``. 

  
Installing Salmon
-----------------

The latest release of Salmon is available either as a pre-compiled binary from GitHub, or via linuxbrew (thanks @sjackman!), we'll grab a pre-compiled binary directly, or you can install via linuxbrew if you want. We can download it using ``wget`` like so::

  > wget --no-check-certificate 'https://github.com/COMBINE-lab/salmon/releases/download/v0.7.0/Salmon-0.7.0_linux_x86_64.tar.gz' 

and we can untar and unzip the resulting file with the following command::

  > tar xzf Salmon-0.7.0_linux_x86_64.tar.gz

Finally, so that we can simply type ``salmon`` to execute salmon, we'll add the appropriate directory to our path variable again.::

  > echo 'export PATH="/home/ubuntu/SalmonBeta-0.7.0_linux_x86_64/bin:$PATH"' >>~/.bashrc

Running Salmon
--------------

"""""""""""""""""""""""""
Creating the Salmon index
"""""""""""""""""""""""""

Since Salmon uses `quasi-mapping <http://bioinformatics.oxfordjournals.org/content/32/12/i192.abstract>`_ behind the scenes, we'll need to build an index on the transcriptome.  Building the ``salmon`` index is relatively quick, we do it with the following command::

  > salmon index -t ref/dmel-all-transcript-r6.12.fasta -i salmon_index

The ``-t`` option tells ``salmon`` where to look for the transcript sequences and ``-i`` tells it where to write the index.


"""""""""""""""""""""""
Quantifying with Salmon
"""""""""""""""""""""""

Now, we'll run Salmon on all of our samples.  We're let salmon use defaults for almost all parameters, but I'll explain the
options and their arguments below.  It will be rather burdensome to run salmon by hand for each sample, so we'll write a small
shell script to run each of the samples one-by-one.  Here's the shell script we'll use::

  #!/bin/bash
  
  for fn in /mnt/reads/*R1_001.fastq.gz
  do

  # get the path to the file
  dir=`dirname $fn`;
  
  # get just the file (without the path)
  base=`basename $fn`;

  # the read filename, without the _R1_001.fastq.gz suffix
  rf=${base%_R1_001.fastq.gz};
  
  # Do whatever we want with it
  salmon quant -i salmon_index -p 8 -l IU -1 <(gunzip -c ${dir}/${rf}_R1_001.fastq.gz) -2 <(gunzip -c ${dir}/${rf}_R2_001.fastq.gz) -o quants/${rf}

  done

  
The call to ``salmon`` takes a few arguments; almost all of them required:

* **-i** tells ``salmon`` where to look for the index
* **-p** tells ``salmon`` how many threads to use
* **-l** tells ``salmon`` the type of the read library (here, inward facing, unstranded reads).  For a more in-depth description of the library types
  and how to specify them in ``salmon``, have a look `here <http://salmon.readthedocs.io/en/develop/library_type.html>`_ in the docs.
* **-1** similar to RapMap, this tells ``salmon`` where to find the first reads of the pair
* **-2** tells ``salmon`` where to find the second reads of the pair
* **-o** tells ``salmon`` where (the directory) to write the output for this sample.  The directory (and the path to it) will be created if it doesn't exist.

  
.. attention::

   We are quantifying *all* 12 samples here.  This totals ~400 -- 500 million read pairs (~800M --- 1B individual reads).
   Salmon will take ~4 minutes per sample, so this process should take 40 - 50 minutes.  This is a good time for us
   to chat, or for you to ask questions you may have thought of during the lecture or up until this point in the
   practical.


""""""""""""""""""""""""""""""""""""
Taking a look at the quantifications
""""""""""""""""""""""""""""""""""""

For Ian's lecture on differential expression, you'll need the quantification results on your local machine, so let's pull them down::

  > scp -i ?????????.pem -r ubuntu@XX.XX.XX.XX:~/quants .

This will copy the ``quants`` directory, recursively, from the server to your local machine.  Let's take a quick peek at some of the quantification results (we'll use R).  Open up RStudio, and set the current directory as the working directory.  We'll do some "sanity checks" using the commands `here <https://github.com/ngs-docs/angus/blob/2016/rob_quant/sanity_check.R>`_ (*please don't make fun of my lack of R-fu --- I'm a Pythonista*).

Note: A pre-computed version of the quantification results is available `here <http://angus.readthedocs.io/en/2016/_static/quantification.tgz>`_.

**TERMINATE YOUR INSTANCE!!!**
