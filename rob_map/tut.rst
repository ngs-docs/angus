.. highlight:: rst

Aligning (and Mapping) RNA-seq reads
====================================

During this lab, you'll learn how to *align* RNA-seq reads to the genome using `STAR <https://github.com/alexdobin/STAR>`_, and
how to *map* RNA-seq reads directly to the transcriptome using `RapMap <https://github.com/COMBINE-lab/RapMap>`_.


Preparing the instance (creating and attaching a volume)
--------------------------------------------------------

Log into your AWS console, and go to the `Volumes` tab.  Create a new volume (set the size equal to 100GiB).  Set the Snapshot ID to `snap-996e1b03`, and make sure to set the zone appropriately.  

Now, create a new instance (**c4.2xlarge**), and navigate to the `storage` tab.  Click `Edit storage`, and then `Add New Volume`.  From the `device` column, select `/dev/sdf`, and in the `Snapshot` column, enter our snapshot id.  When you're ready, click `Review and Launch`.

Log into your instance
-----------------------

For this tutorial, we'll use a **c4.2xlarge** instance.  Make sure you create the instance with the
volume containing the reads attached!::

  > ssh -i ~/Downloads/?????.pem ubuntu@XX.XX.XX.XX

Update the package list
-----------------------

::
   
   > sudo apt-get update

Install some base packages
--------------------------

First, install the "build tools" (compilers etc. that may be needed)::

  > sudo apt-get install build-essential

To make use of linuxbrew, we'll need Ruby and Git:::

  > sudo apt-get install git ruby

Install ``tmux`` so that we can have multiple terminals in the same session::

  > sudo apt-get install tmux

Start your tmux session::

  > tmux new -s quant
  
  
Do the necessary Linuxbrew incantations
---------------------------------------

Get linuxbrew with the following command:::

  > ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"

  
when prompted, hit `RETURN`.  To make the software we'll install via linuxbrew accessable we have 
to place the directory where linuxbrew installs programs into our ``PATH``.  The following sequence of 
commands will do this:::

  > echo 'export PATH="/home/ubuntu/.linuxbrew/bin:$PATH"' >>~/.bashrc

  > echo 'export MANPATH="/home/ubuntu/.linuxbrew/share/man:$MANPATH"' >>~/.bashrc

  > echo 'export INFOPATH="/home/ubuntu/.linuxbrew/share/info:$INFOPATH"' >>~/.bashrc

To make our new paths active, we have to ``source`` them:::

  > source ~/.bashrc

The programs we're interested in installing are part of hombrew-science.  We can "tap" the science keg (;P) as follows:::

  > brew tap homebrew/science
  
Now, install samtools:::

  > brew install homebrew/science/samtools

Mounting the reads
------------------

We have prepared (thanks; `@monsterbashseq! <https://ljcohen.github.io/>`_) an Amazon volume from which you can load the reads directly.  When we created our AWS instance, we attached the volume with the reads to ``/dev/xvdf``.  We have to *mount* this device.  First, we'll create a place to mount it:::

  > sudo mkdir -p /mnt/reads

Now, we actualy mount the device at the mount point:::

  > sudo mount /dev/xvdf /mnt/reads

When this command finishes (should only take a few seconds) we're good to go, but just need to change the permissions on this folder.::

  > sudo chown -R ubuntu:ubuntu /mnt/reads

Now all of the read files should be available in ``/mnt/reads``.  Check this out with:::

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

We'll be aligning our reads to the Drosophila genome with STAR, and against the Drosophila transcriptome with RapMap.  So, we'll need these files (we'll also need the GTF file corresponding to the genome, as STAR uses this to index known splice sites).

Grab the genome:::

  > wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/fasta/dmel-all-chromosome-r6.12.fasta.gz

and the annotation for the genome::

  > wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/gtf/dmel-all-r6.12.gtf.gz

and the transcriptome::

  > wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/fasta/dmel-all-transcript-r6.12.fasta.gz

We'll put all of these in a folder called ``ref``, and, since they're fairly small, we'll unzip them all::

  > mkdir ref
  > mv *.gz ref
  > cd ref
  > gunzip *.gz
  > cd ..
  
Great; now, we're ready to grab our aligner and align some reads ... almost!  For a reason that I'd be happy to explain if anyone is actually curious, we have to "warm up" the file's we'll be aligning.  We can do this as follows::

  > wc -l <(gunzip -c /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R1_001.fastq.gz)
  > wc -l <(gunzip -c /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R2_001.fastq.gz)

Now, we're *actually* ready to do our alignment / mapping.

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

In order to align reads efficiently, STAR has to build an index (in this case, a suffix array), over the genome.  First, we'll create the directory where the index will live:::

  > mkdir star_index

Now, we have to tell STAR to build the index using our reference genome and the GTF annotation.  That command looks like::

  > ~/STAR-2.5.2a/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate \
        --genomeDir star_index --genomeFastaFiles ref/dmel-all-chromosome-r6.12.fasta \
        --sjdbGTFfile ref/dmel-all-r6.12.gtf --sjdbOverhang 99

Here, we're telling STAR that it can use up to 8 threads, and it should build the index on the genome and using the annotation we provide.  The ``sjdbOverhang`` parameter is helpful in setting some internal options, and is recommended to be set as read_lenght - 1.  Once you execute this command, ``STAR`` should run for ~3 minutes before finishing and placing the index in the requested directory.  We can check the contents of the index file:::

  > ls -lha star_index

And we'll see a bunch of files related to the index built by STAR.  The total size of this folder should be ~3.3G --- quite a bit larger than the input reference genome (140M).  This is one of the trade-offs that STAR makes; to provide very fast alignment speeds (and STAR is *very* fast), it uses a large amount of memory.  When aligning to e.g. the human genome, you should be prepared to have 20-30G of RAM available for STAR.

""""""""""""""""""""""""""""
Aligning the reads with STAR
""""""""""""""""""""""""""""

STAR has a *dizzying* array of options. You can find most of them explained in detail in the `STAR manual <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>`_.  For the purposes of this lesson, we'll keep things simple and I'll explain the main options we're setting here.  First, let's create the output directory where our alignments will live:::

  > mkdir alignments

Now, we'll run STAR to align the reads.  Before we do this, we'll create another tmux window that we'll use to monitor STAR's progress.  Create a new
``tmux`` window with `ctrl + b, c`.  Now, we'll start running STAR in this window using the following command::

  > /usr/bin/time ~/STAR-2.5.2a/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir star_index \
      --readFilesIn /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R1_001.fastq.gz /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R2_001.fastq.gz \
      --readFilesCommand gunzip -c --outFileNamePrefix alignments/ORE_sdE3_r1_GTGGCC_L004 --outSAMtype BAM Unsorted

The relevant options are:

* **--runThreadN** tells STAR to use the specified number of threads
* **--genomeDir** tells STAR where to look for the index
* **--readFilesIn** tells STAR where to find the files it will be aligning.  When aligning paired-end reads, the left and right reads are separated by a space.  If there are multiple files for the left and right reads, these are separated by commas.  It is important, if you have multiple left and right files, that they are given in the *same order* in their respective lists.
* **--readFilesCommand** tells star what command it should use to coerce the input into standard FASTA/FASTQ format.  Here, since our reads are gzipped, we tell STAR to use ``gzip -c`` to produce a ``FASTQ`` format file from the input ``fastq.gz`` format files.
* **--outFileNamePrefix** tells STAR how it should name its output files.  There are defaults, but here we override them with the name of the library we're aligning
* **--outSAMtype** tells STAR the format in which the output should be written.  Here, we're telling STAR that the output should be in BAM format (binary and compressed), and that it's OK for the alignments to be unsorted by position / name / etc.

This command will take a little while to run.  While STAR is doing it's thing, we can monitor it's progress.  First, move over to the other window with
`ctrl + b, p`.  Then we'll monitor STAR's progress with this nifty little command:::

  > tail -f alignments/ORE_sdE3_r1_GTGGCC_L004Log.progress.out

The ``tail -f`` command will *follow* the file, and will write the end (tail) of the  file to the console when it's updated.  At this point, while we wait, it would be an ideal time to discuss what STAR is doing, or to answer any questions you might have.

Using RapMap
------------

""""""""""""""""
Obtaining RapMap
""""""""""""""""

Just as with STAR, we'll grab the latest RapMap binary from GitHub.  I'm the maintainer of RapMap, and I update the software on a somewhat regular basis (though not as often as Alex updates STAR).  If you're going to start a new analysis using RapMap, it's probably worth checking for the latest version.  We can grab the latest pre-compiled binary like so:::

  > wget --no-check-certificate 'https://drive.google.com/uc?export=download&id=0B3iS9-xjPftjbkJoeVFscnUwc00' -O RapMap-v0.3.1_CentOS5.tar.gz

and then we expand the tarball::

  > tar xzvf RapMap-v0.3.1_CentOS5.tar.gz 

add RapMap to the path so we can just call the executable.::
  
  > export PATH=/home/ubuntu/RapMap-0.3.1_CentOS5/bin:$PATH

  
"""""""""""""""""""""""""
Building the RapMap Index
"""""""""""""""""""""""""

Like STAR, RapMap will need an index in order to map reads efficiently.  Unlike STAR, however, RapMap will build an index over the transcriptome rather than the entire genome.  We can build RapMap's index as follows:::

  > rapmap quasiindex -t ref/dmel-all-transcript-r6.12.fasta -i rapmap_index

Unlike STAR, RapMap will create the index folder if it doesn't exist yet, so there's no need to create it first.


"""""""""""""""""""""""""""""
Mapping the reads with RapMap
"""""""""""""""""""""""""""""

Now, we'll map the same set of reads we aligned above with STAR, but we'll map them to the transcriptome using RapMap.  The commands we'll use for this is::

  > mkdir mappings
  > rapmap quasimap -i rapmap_index -t 8 -1 <(gunzip -c /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R1_001.fastq.gz) \
  -2 <(gunzip -c /mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R2_001.fastq.gz) | samtools view -Sb -@4 - > mappings/mapped_reads.bam

Though RapMap itself has far fewer options than STAR, there's still quite a bit going on above.  Let's unpack this command; first the RapMap options:

* **-i** tells RapMap where to find the index
* **-t** tells RapMap how many threads it can use
* **-1** tells RapMap where to find the first set of reads for a paired-end library
* **-2** tells RapMap where to find the second set of reads for a paired-end library

By default, RapMap will write it's output, in SAM format, to ``stdout``.  Here, this is what we want, but if you want the output to be redirected to a file, that can be done with the ``-o`` option.  One thing to make note of here is the ``<()`` syntax we're using to provide the read files.  Actually, this version of RapMap can accept compressed files directly, but I just wanted to show off the process substitution syntax.  However, we achieve this differently than STAR.  Here, we use `process substitution <http://tldp.org/LDP/abs/html/process-sub.html>`_ to directly create a pipe from which RapMap will read the decompressed sequences.  The process substitution syntax runs the command within the ``<()`` in a separate process, and writes the output to a file descriptor (e.g. ``/dev/fd/<n>``).  RapMap reads the input from this file descriptor as if it were a normal file.  Generally, this process substitution syntax is *insanely* useful.  I highly recommend you learn to become comfortable with it, as it can make many processing tasks much easier.

After the RapMap command we are piping the output to ``samtools``.  Since RapMap does not (yet) have the built-in ability to write to BAM format, we're using ``samtools`` to convert the SAM output to BAM format on-the-fly.  The command we're using tells ``samtools`` to read input as SAM and write output as BAM, to use up to 4 threads for compression (don't worry that 8 + 4 > 8).  The ``-`` tells ``samtools`` to read from stdin and it, by default, writes its output to stdout.  We then pipe this output directly to the file we wish to create.  Now, we wait for RapMap to finish.  It should be a bit faster than STAR, though in my testing on the AWS instance, it's bottlenecked by the SAM -> BAM conversion, and so won't be able to make use of all the cores we're allowing it to.


Looking at the results
----------------------

Now, we've aligned the reads to the genome with STAR, and mapped them to the transcriptome with RapMap.  We can do a quick comparison of these BAM files::

  > ls -lha alignments/ORE_sdE3_r1_GTGGCC_L004Aligned.out.bam
  -rw-rw-r-- 1 ubuntu ubuntu 7.3G Aug 11 03:57 alignments/ORE_sdE3_r1_GTGGCC_L004Aligned.out.bam
  > ls -lha mappings/mapped_reads.bam
  -rw-rw-r-- 1 ubuntu ubuntu 7.2G Aug 11 04:14 mappings/mapped_reads.bam

The files are about the same size (this won't always be the case), but actually contain *very* different information.

.. raw:: html
	 
	 <details> 
	 <summary>Q: Why is the information so different?</summary>
	 A: STAR is aligning to the <b>genome</b> while RapMap is aligning to the <b>transcriptome</b>.  Thus,
	    STAR's BAM file will contained spliced alignments, while RapMap's won't.  Conversely, We expect
	    RapMap's BAM file to encode many more <i>multimapping</i> reads than STAR's BAM file.
	 </details>


Let's get some details about the mappings from the BAM files.  We can use samtools' ``flagstat`` command for this::

  > samtools flagstat alignments/ORE_sdE3_r1_GTGGCC_L004Aligned.out.bam
  76410744 + 0 in total (QC-passed reads + QC-failed reads)
  7985862 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  76410744 + 0 mapped (100.00% : N/A)
  68424882 + 0 paired in sequencing
  34212441 + 0 read1
  34212441 + 0 read2
  68424882 + 0 properly paired (100.00% : N/A)
  68424882 + 0 with itself and mate mapped
  0 + 0 singletons (0.00% : N/A)
  0 + 0 with mate mapped to a different chr
  0 + 0 with mate mapped to a different chr (mapQ>=5)

  > samtools flagstat mappings/mapped_reads.bam
  182443558 + 0 in total (QC-passed reads + QC-failed reads)
  116396448 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  178493051 + 0 mapped (97.83% : N/A)
  66047110 + 0 paired in sequencing
  33023555 + 0 read1
  33023555 + 0 read2
  63599806 + 0 properly paired (96.29% : N/A)
  63599806 + 0 with itself and mate mapped
  1223652 + 0 singletons (1.85% : N/A)
  0 + 0 with mate mapped to a different chr
  0 + 0 with mate mapped to a different chr (mapQ>=5)


From the number of alignments, you can see that the multimapping rate of RapMap is higher than that of STAR. If we assume that they aligned the same number of reads (they *didn't*, but the numbers are close), then there are, on average, 2.34 RapMap mappings for every STAR alignment --- multimapping in the transcriptome is *prevalent*.

.. raw:: html

	 <details>
	 <summary>Q: How may reads were there in the input?</summary>
	 A: We can compute this with
	 <p><code class="docutils literal"><span class="pre">&gt;bc</span> <span class="pre">-l</span> <span class="pre">&lt;&lt;&lt;</span> <span class="pre">&quot;$(gunzip</span> <span class="pre">-c</span> <span class="pre">/mnt/reads/ORE_sdE3_r1_GTGGCC_L004_R1_001.fastq.gz</span> <span class="pre">|</span> <span class="pre">wc</span> <span class="pre">-l)</span> <span class="pre">/</span> <span class="pre">4&quot;</span></code></p>
	 There are 36,968,390 reads in the original file.


"""""""""""""""""""""""""""""""
Preparing for tomorrow's lesson
"""""""""""""""""""""""""""""""

Tomorrow, rather than just dealing with a single pair of read files, we'll quantify *all* of our samples using the ultra-fast quantification tool, Salmon.  In order to prepare your Amazon instance, you'll have to "warm up" the volume with the data.  To do this, you'll have to run the following command, and let it execute until it completes::
  
  >sudo fio --filename=/dev/xvdf --rw=randread --bs=128k --iodepth=32 --ioengine=libaio --direct=1 --name=volume-initialize

This will tell the instance to "touch" all of the blocks on this volume, which will make sure they are physically present, and not just "accessible".  This command make take 1.5 - 2 hours to complete.  Once this is done, we will **stop** the instance.  *Note:*  **Do not terminate your instance, stop it instead**, otherwise you'll loose all of the work that was just done to warm the volume.
