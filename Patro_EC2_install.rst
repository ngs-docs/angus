.. highlight:: rst

Aligning (and Mapping) RNA-seq reads
====================================

During this lab, you'll learn how to *align* RNA-seq reads to the genome using `STAR <https://github.com/alexdobin/STAR>`_, and
how to *map* RNA-seq reads directly to the transcriptome using `RapMap <https://github.com/COMBINE-lab/RapMap>`_.

``> ssh -i ~/Downloads/?????.pem ubuntu@XX.XX.XX.XX``

First, install the "build tools" (compilers etc. that may be needed)

``> sudo apt-get update``

``> sudo apt-get install build-essential``

To make use of linuxbrew, we'll need Ruby:

``> sudo apt-get install git ruby``

commands::
  
  > sudo mkdir -p /mnt/reads
  > sudo mount /dev/xvdf /mnt/reads
  > chown -R ubuntu:ubuntu /mnt/reads
  > mkdir mapping && cd mapping
  > wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/fasta/dmel-all-chromosome-r6.12.fasta.gz
  > wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/gtf/dmel-all-r6.12.gtf.gz
  > wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/fasta/dmel-all-transcript-r6.12.fasta.gz
  > mkdir ref
  > mv *.gz ref
  > cd ref
  > gunzip *.gz
  > cd ..
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
