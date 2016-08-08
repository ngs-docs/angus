================================================
Quality Trimming
================================================

During this lab, we will acquaint ourselves with the software packages
Trimmomatic, khmer and Jellyfish. Your objectives are:

1. Familiarize yourself with software, how to install and execute it and optionally how to
   visualize results.
2. Characterize sequence quality.

The Skewer manual: https://github.com/relipmoc/skewer

The JellyFish manual: http://www.genome.umd.edu/jellyfish.html




Step 1: Launch and AMI. For this exercise, we will use a **c4.2xlarge** with a 500Gb EBS volume. Remember to change the permission of your key code ``chmod 400 ~/Downloads/????.pem`` (change ????.pem to whatever you named it)

::

    ssh -i ~/Downloads/?????.pem ubuntu@XX.XX.XX.XX



**Update Software**

::

    sudo apt-get update



**Install updates**

::

    sudo apt-get -y upgrade


**Install other software** Note that you can install a large amount of software from the Ubuntu "App Store" using a single command. Some of this software we will not use for this tutorial, but...

::

    sudo apt-get -y install tmux git gcc make g++ python-dev unzip default-jre build-essential libcurl4-openssl-dev \
    zlib1g-dev python-pip

**Install Ruby**  Ruby is a computer language like Python or Perl.

::

    cd
    wget https://keybase.io/mpapis/key.asc
    gpg --import key.asc
    \curl -sSL https://get.rvm.io | bash -s stable --ruby
    source /home/ubuntu/.rvm/scripts/rvm

**Install Brew** Brew is a piece of software the serves as a 'package manager'. It makes installing software easy! You can use it for lots of things, but not everything. Knowing it's limitations will come with time.

::

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"

    # press return


    echo 'export PATH="/home/ubuntu/.linuxbrew/bin:$PATH"' >>~/.profile
    echo 'export MANPATH="/home/ubuntu/.linuxbrew/share/man:$MANPATH"' >>~/.profile
    echo 'export INFOPATH="/home/ubuntu/.linuxbrew/share/info:$INFOPATH"' >>~/.profile
    source ~/.profile


**Install Bioinformatics Packages via Brew** These are the packages that we will use to do real work!!! YAY!!!

::

    brew tap homebrew/science
    brew install jellyfish
    brew install skewer
    brew install fastqc


**Install khmer**

::

    easy_install --user -U setuptools
    pip install --user khmer
    echo 'export PATH="/home/ubuntu/.local/bin/:$PATH"' >>~/.profile
    source ~/.profile

**Download data**: For this lab, we'll be using files from Jack Gilbert's Merlot wine study (http://mbio.asm.org/content/6/2/e02527-14.full). The details are not important right now, but this is a metagenomic sample from root of the grape vine.

You are downloading from MG-RAST, which is a popular metagenomics analysis package. There are a lot of places to get raw data.

::

   mkdir $HOME/reads
   cd $HOME/reads/

   curl http://api.metagenomics.anl.gov//download/mgm4520306.3?file=050.1 > root_S13.R1.fq

   curl http://api.metagenomics.anl.gov//download/mgm4520307.3?file=050.1 > root_S13.R2.fq

--------------

**Do 2 different trimming levels -- Phred=2 and Phred=30**: One of these is very harsh, the other is probably more appropriate.  Which one is which?

Look at the output from this command, which should start with ``Input Read Pairs:``

::

    mkdir $HOME/trimming
    cd $HOME/trimming


    curl -LO https://s3.amazonaws.com/gen711/TruSeq3-PE.fa

    skewer -l 25 -m pe -o skewerQ2 --mean-quality 2 --end-quality 2 -t 16 \
    -x TruSeq3-PE.fa \
    $HOME/reads/root_S13.R1.fq $HOME/reads/root_S13.R2.fq

    #and

    skewer -l 25 -m pe -o skewerQ30 --mean-quality 30 --end-quality 30 -t 16 \
    -x TruSeq3-PE.fa \
    $HOME/reads/root_S13.R1.fq $HOME/reads/root_S13.R2.fq


**Interleave reads**

::

    interleave-reads.py skewerQ2-trimmed-pair1.fastq skewerQ2-trimmed-pair2.fastq > Q2.interleave.fq
    interleave-reads.py skewerQ30-trimmed-pair1.fastq skewerQ30-trimmed-pair2.fastq > Q30.interleave.fq


**Run Jellyfish**

::

  mkdir $HOME/jelly
  cd $HOME/jelly


  jellyfish count -m 25 -s 200M -t 16 -C -o trim30.jf $HOME/trimming/Q30.interleave.fq
  jellyfish histo trim30.jf -o trim30.histo

  #and

  jellyfish count -m 25 -s 200M -t 16 -C -o trim2.jf $HOME/trimming/Q2.interleave.fq
  jellyfish histo trim2.jf -o trim2.histo

--------------


**Look at the 2 histograms**

::

  head *histo

--------------

**Run FastQC on your data**

::

  mkdir $HOME/fastqc
  cd $HOME/fastqc

  fastqc -t 16 $HOME/trimming/Q2.interleave.fq
  fastqc -t 16 $HOME/trimming/Q30.interleave.fq
  ls -lth

**Download FastQC .zip file to your computer**

Open up a new terminal window using the buttons command-t, then unzip as per normal.

::

  scp -i ~/Downloads/????.pem ubuntu@??-???-???-?:/home/ubuntu/trimming/*zip ~/Downloads/


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

    plot(log(trim2$V2[2:100] - trim30$V2[2:100]), type='l',
     xlim=c(0,100), xaxs="i", yaxs="i", frame.plot=F,
     ylim=c(0,20), col='red', xlab='kmer frequency',
     lwd=4, ylab='log diff count',
     main='Log Diff in 25mer counts of freq 1 to 100 \n Phred2 vs. Phred30')
