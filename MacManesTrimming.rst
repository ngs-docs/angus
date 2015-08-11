================================================ 
Quality Trimming
================================================

During this lab, we will acquaint ourselves with the software packages
Trimmomatic, khmer and Jellyfish. Your objectives are:

1. Familiarize yourself with software, how to install and execute it and optionally how to
   visualize results.
2. Characterize sequence quality.

The Trimmomatoc manual: http://www.usadellab.org/cms/?page=trimmomatic

The JellyFish manual: http://www.genome.umd.edu/jellyfish.html

The Khmer webpage: http://khmer.readthedocs.org/en/v1.4.1/

--------------

Step 1: Launch and AMI. For this exercise, we will use a **c4.2xlarge** with a 500Gb EBS volume. Remember to change the permission of your key code ``chmod 400 ~/Downloads/????.pem`` (change ????.pem to whatever you named it)

::

    ssh -i ~/Downloads/?????.pem ubuntu@ec2-???-???-???-???.compute-1.amazonaws.com

--------------

**Update Software**

::

    sudo bash
    apt-get update

--------------

**Install updates**

::

    apt-get -y upgrade

--------------

**Install other software** Note that you can install a large amount of software from the Ubuntu "App Store" using a single command. Some of this software we will not use for this tutorial, but...

::

    apt-get -y install tmux git gcc make g++ python-dev unzip default-jre build-essential libcurl4-openssl-dev zlib1g-dev python-pip fastqc

--------------

**Mount hard drive** The EBS volume we asked for is not automatically mounted - we need to do that. 

::

    df -h
    mkfs -t ext4 /dev/xvdb  
    mount /dev/xvdb /mnt  
    chown -R ubuntu:ubuntu /mnt  
    df -h

--------------

**Install Trimmomatic**: Trimmomatic is my favorite tool for adapter and quality trimming, there are several others available. 

::

    cd $HOME
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
    unzip Trimmomatic-0.33.zip
    cd Trimmomatic-0.33
    chmod +x trimmomatic-0.33.jar

--------------

**Install Jellyfish**: Jellyfish is a kmer counter that is very popular. Note that most (all?) of what Jellyfish can do, khmer can do (better/faster?). I'll introduce it simply to let you know it exists.

::

    cd $HOME
    wget ftp://ftp.genome.umd.edu/pub/jellyfish/jellyfish-2.1.3.tar.gz
    tar -zxf jellyfish-2.1.3.tar.gz
    cd jellyfish-2.1.3/
    ./configure
    make
    PATH=$PATH:$(pwd)/bin

--------------

**Install khmer**: This is The CTB lab software. 

::

    pip install --upgrade setuptools
    pip install khmer

--------------

**Download data**: For this lab, we'll be using files from Jack Gilbert's Merlot wine study (http://mbio.asm.org/content/6/2/e02527-14.full). The details are not important right now, but this is a metagenomic sample from root of the grape vine.

You are downloading from MG-RAST, which is a popular metagenomics analysis package. There are a lot of places to get raw data.

::

   mkdir /mnt/reads 
   cd /mnt/reads/

   curl http://api.metagenomics.anl.gov//download/mgm4520306.3?file=050.1 > root_S13.R1.fq

   curl http://api.metagenomics.anl.gov//download/mgm4520307.3?file=050.1 > root_S13.R2.fq

--------------

**Do 2 different trimming levels -- Phred=2 and Phred=30**: One of these is very harsh, the other is probably more appropriate.  Which one is which?

Look at the output from this command, which should start with ``Input Read Pairs:``

::

    mkdir /mnt/trimming
    cd /mnt/trimming

    #paste the below lines together as 1 command

    java -Xmx10g -jar $HOME/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
    -threads 8 -baseout root_S13.Phred2.fq \
    /mnt/reads/root_S13.R1.fq \
    /mnt/reads/root_S13.R2.fq \
    ILLUMINACLIP:$HOME/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
    SLIDINGWINDOW:4:2 \
    LEADING:2 \
    TRAILING:2 \
    MINLEN:25

    #and

    java -Xmx10g -jar $HOME/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
    -threads 8 -baseout root_S13.Phred30.fq \
    /mnt/reads/root_S13.R1.fq \
    /mnt/reads/root_S13.R2.fq \
    ILLUMINACLIP:$HOME/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
    SLIDINGWINDOW:4:30 \
    LEADING:30 \
    TRAILING:30 \
    MINLEN:25


--------------

**Run khmer and Jellyfish**

::

  interleave-reads.py root_S13.Phred30_1P.fq root_S13.Phred30_2P.fq > root_S13.Phred30.interleaved.fq

  interleave-reads.py root_S13.Phred2_1P.fq root_S13.Phred2_2P.fq > root_S13.Phred2.interleaved.fq

  mkdir /mnt/jelly
  cd /mnt/jelly


  jellyfish count -m 25 -s 200M -t 8 -C -o trim2.jf /mnt/trimming/root_S13.Phred2.interleaved.fq
  jellyfish histo trim2.jf -o trim2.histo

  #and

  jellyfish count -m 25 -s 200M -t 8 -C -o trim30.jf /mnt/trimming/root_S13.Phred30.interleaved.fq
  jellyfish histo trim30.jf -o trim30.histo

--------------


**Look at the 2 histograms**

::

  head *histo

--------------

**Run FastQC on your data**

::

  mkdir /mnt/fastqc
  cd /mnt/fastqc

  fastqc -t 8 /mnt/reads/root_S13.R1.fq /mnt/reads/root_S13.R2.fq
  fastqc -t 8 /mnt/trimming/root_S13.Phred30_1P.fq /mnt/trimming/root_S13.Phred30_2P.fq
  fastqc -t 8 /mnt/trimming/root_S13.Phred2_1P.fq /mnt/trimming/root_S13.Phred2_2P.fq
  ls -lth
  
**Download FastQC .zip file to your computer**

Open up a new terminal window using the buttons command-t, then unzip as per normal. 

::

  scp -i ~/Downloads/????.pem ubuntu@ec2-??-???-???-??.compute-1.amazonaws.com:/mnt/reads/*zip ~/Downloads/

  scp -i ~/Downloads/????.pem ubuntu@ec2-??-???-???-??.compute-1.amazonaws.com:/mnt/trimming/*zip ~/Downloads/


--------------


**WON'T COVER THE STUFF BELOW, THOUGH YOU SHOULD TRY TO DO IT**

Now look at the ``.histo`` file, which is a kmer distribution. I want you to plot the distribution using R and RStudio.

**OPEN RSTUDIO**: Google and install locally. There are OSX and Windows versions. 

Open up a new terminal window using the buttons command-t

::

  scp -i ~/Downloads/????.pem ubuntu@ec2-??-???-???-??.compute-1.amazonaws.com:/mnt/jelly/*histo ~/Downloads/


Import and visualize the 2 histogram datasets:

::

    trim2 <- read.table("~/Downloads/trim2.histo", quote="\"")
    trim30 <- read.table("~/Downloads/trim30.histo", quote="\"")

    #Plot: Make sure and change the names to match what you import.
    #What does this plot show you?? 

    barplot(c(trim2$V2[1],trim30$V2[1]),
        names=c('Phred2', 'Phred30'),
        main='Number of unique kmers')

    # plot differences between non-unique kmers

    plot(trim2$V2[2:30] - trim30$V2[2:30], type='l',
        xlim=c(1,5), xaxs="i", yaxs="i", frame.plot=F,
        ylim=c(0,20000000), col='red', xlab='kmer frequency',
        lwd=4, ylab='count',
        main='Diff in 25mer counts of freq 1 to 5 \n Phred2 vs. Phred30')
