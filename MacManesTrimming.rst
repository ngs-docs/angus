================================================ 
Quality Trimming
================================================

--

--------------

During this lab, we will acquaint ourselves with the software packages
Trimmomatic and JellyFish. Your objectives are:

-  

1. Familiarize yourself with the software, how to execute it, how to
   visualize results.

2. Regarding your dataset. Characterize sequence quality.

The Trimmomatoc manual: http://www.usadellab.org/cms/?page=trimmomatic

The JellyFish manual: http://www.genome.umd.edu/jellyfish.html

The Khmer webpage: http://khmer.readthedocs.org/en/v1.4.1/

--------------

    Step 1: Launch and AMI. For this exercise, we will use a c4.2xlarge
    with a 500Gb EBS volume. Remember to change the permission of your
    key code ``chmod 400 ~/Downloads/????.pem`` (change ????.pem to
    whatever you named it)

::

    ssh -i ~/Downloads/?????.pem ubuntu@ec2-???-???-???-???.compute-1.amazonaws.com

--------------

    Update Software

::

    sudo bash
    apt-get update

--------------

    Install updates

::

    apt-get -y upgrade

--------------

    Install other software

::

    apt-get -y install tmux git gcc make g++ python-dev unzip default-jre build-essential libcurl4-openssl-dev zlib1g-dev python-pip

--------------

    Mount hard drive

::

    df -h
    mkfs -t ext4 /dev/xvdb  
    mount /dev/xvdb /mnt  
    chown -R ubuntu:ubuntu /mnt  
    df -h

--------------

    Install Trimmomatic

::

    cd $HOME
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
    unzip Trimmomatic-0.33.zip
    cd Trimmomatic-0.33
    chmod +x trimmomatic-0.33.jar

--------------

    Install Jellyfish

::

    cd $HOME
    wget ftp://ftp.genome.umd.edu/pub/jellyfish/jellyfish-2.1.3.tar.gz
    tar -zxf jellyfish-2.1.3.tar.gz
    cd jellyfish-2.1.3/
    ./configure
    make
    PATH=$PATH:$(pwd)/bin

--------------

    Install khmer

::

    pip install --upgrade setuptools
    pip install khmer

--------------

    Download data. For this lab, we'll be using files from Jack
    Gilbert's Merlot wine study
    (http://mbio.asm.org/content/6/2/e02527-14.full). The details are
    not important right now, but this is a metagenomic sample from root
    of the grape vine.

-  mkdir /mnt/reads cd /mnt/reads/

   curl
   http://api.metagenomics.anl.gov//download/mgm4520306.3?file=050.1 >
   root\_S13.R1.fq

   curl
   http://api.metagenomics.anl.gov//download/mgm4520307.3?file=050.1 >
   root\_S13.R2.fq

--------------

    Do 2 different trimming levels -- 2 and 30.

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


    java -Xmx10g -jar $HOME/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
    -threads 8 -baseout root_S13.Phred30.fq \
    /mnt/reads/root_S13.R1.fq \
    /mnt/reads/root_S13.R2.fq \
    ILLUMINACLIP:$HOME/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
    SLIDINGWINDOW:4:30 \
    LEADING:30 \
    TRAILING:30 \
    MINLEN:25

**WON'T COVER THE STUFF BELOW, THOUGH YOU SHOULD TRY TO DO IT**

    Now look at the ``.histo`` file, which is a kmer distribution. I
    want you to plot the distribution using R and RStudio.

    OPEN RSTUDIO

::

    #Import all 3 histogram datasets: this is the code for importing 1 of them..

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
