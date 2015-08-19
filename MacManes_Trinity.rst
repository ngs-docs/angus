================================================ 
Trinity and Transcriptome Evaluation
================================================

Trinity: http://trinityrnaseq.github.io/

Transrate: http://hibberdlab.com/transrate/installation.html

--------------

Step 1: Launch and AMI. For this exercise, we will use a **c4.4xlarge** with a 500Gb EBS volume. Remember to change the permission of your key code ``chmod 400 ~/Downloads/????.pem`` (change ????.pem to whatever you named it)

::

    ssh -i ~/Downloads/?????.pem ubuntu@ec2-???-???-???-???.compute-1.amazonaws.com

--------------

**Update Software**

::

    sudo apt-get update

--------------

**Install updates**

::

    sudo apt-get -y upgrade

--------------

**Install other software** Note that you can install a large amount of software from the Ubuntu "App Store" using a single command. Some of this software we will not use for this tutorial, but...

::

  sudo apt-get -y install build-essential tmux git gcc make g++ python-dev unzip default-jre libcurl4-openssl-dev zlib1g-dev python-pip fastqc samtools bowtie ncbi-blast+ hmmer emboss

--------------

**Mount hard drive** The EBS volume we asked for is not automatically mounted - we need to do that. 

::

    sudo mkfs -t ext4 /dev/xvdb  
    sudo mount /dev/xvdb /mnt  
    sudo chown -R ubuntu:ubuntu /mnt  
    df -h

--------------


**INSTALL TRANSRATE**

::
  
  cd $HOME
  curl -LO https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz
  tar -zxf transrate-1.0.1-linux-x86_64.tar.gz
  cd transrate-1.0.1-linux-x86_64
  PATH=$PATH:$(pwd)

--------------


**INSTALL Augustus**

::

  cd $HOME
  curl -O http://augustus.gobics.de/binaries/augustus.2.5.5.tar.gz
  tar -zxf augustus.2.5.5.tar.gz
  cd augustus.2.5.5/
  make
  PATH=$PATH:$(pwd)/bin
  export AUGUSTUS_CONFIG_PATH=$(pwd)/config

--------------

**INSTALL BUSCO**:

::

  cd $HOME
  curl -O http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
  tar -zxf BUSCO_v1.1b1.tar.gz
  cd BUSCO_v1.1b1/
  PATH=$PATH:$(pwd)

--------------

**Install Trinity**

::

  git clone https://github.com/trinityrnaseq/trinityrnaseq.git
  cd trinityrnaseq
  make -j4
  PATH=$PATH:$(pwd)

--------------

**Install Trimmomatic**

::

  cd $HOME
  wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
  unzip Trimmomatic-0.33.zip
  cd Trimmomatic-0.33
  chmod +x trimmomatic-0.33.jar

--------------

**Download data**: For this lab, we'll be using 
::

   mkdir /mnt/reads 
   cd /mnt/reads/

   curl -LO https://www.dropbox.com/s/4o6eduzcw11gz53/subsamp.2.fq.gz

   curl -LO https://www.dropbox.com/s/i4wst01yz10i9x9/subsamp.1.fq.gz

--------------

**Do 2 different trimming levels -- Phred=2 and Phred=30**: One of these is very harsh, the other is probably more appropriate.  Which one is which?

Look at the output from this command, which should start with ``Input Read Pairs:``

::

    mkdir /mnt/trimming
    cd /mnt/trimming

    #paste the below lines together as 1 command

    java -Xmx10g -jar $HOME/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
    -threads 8 -baseout subsamp.Phred2.fq \
    /mnt/reads/subsamp.1.fq.gz \
    /mnt/reads/subsamp.2.fq.gz \
    ILLUMINACLIP:$HOME/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
    SLIDINGWINDOW:4:2 \
    LEADING:2 \
    TRAILING:2 \
    MINLEN:25

    #and

    java -Xmx10g -jar $HOME/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
    -threads 8 -baseout subsamp.Phred30.fq \
    /mnt/reads/subsamp.1.fq.gz \
    /mnt/reads/subsamp.2.fq.gz \
    ILLUMINACLIP:$HOME/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
    SLIDINGWINDOW:4:30 \
    LEADING:30 \
    TRAILING:30 \
    MINLEN:25


--------------

**Run Trinity**

::

  mkdir /mnt/assembly
  cd /mnt/assembly

  #Open tumx window
  
  tmux new -s trinity

  #Phred30 dataset  

  Trinity --seqType fq --max_memory 10G --left /mnt/trimming/subsamp.Phred30_1P.fq \
  --right /mnt/trimming/subsamp.Phred30_2P.fq --CPU 16

  #Phred2 dataset

  Trinity --seqType fq --max_memory 10G --left /mnt/trimming/subsamp.Phred2_1P.fq \
  --right /mnt/trimming/subsamp.Phred2_2P.fq --CPU 16

**Fix Trinity Headers**

::

  sed -i 's_|_-_g' /mnt/assembly/trinity_out_dir/Trinity.fasta

  Control-b d #to exit tmux

--------------



**Run BUSCO for assemblies**: There are Eukaryote, Metazoa, Arthropod, Vertebrate, Plant references for use with other genomes. 

::
  
  
  mkdir /mnt/busco
  cd /mnt/busco

  #Download busco database

  tmux new -s busco
  
  curl -LO http://busco.ezlab.org/files/vertebrata_buscos.tar.gz
  tar -zxf vertebrata_buscos.tar.gz

  python3 /home/ubuntu/BUSCO_v1.1b1/BUSCO_v1.1b1.py \
  -m trans -in /mnt/assembly/trinity_out_dir/Trinity.fasta \
  --cpu 16 -l vertebrata -o trin.assemblty

  less run*/short*

  Control-b d #to exit tmux


--------------

**Run Transrate**

::

  tmux new -s transrate

  mkdir /mnt/transrate
  cd /mnt/transrate
  $HOME/transrate-1.0.1-linux-x86_64/transrate -a /mnt/assembly/trinity_out_dir/Trinity.fasta -t 16 \  
  --left /mnt/trimming/subsamp.Phred30_1P.fq \  
  --right /mnt/trimming/subsamp.Phred30_2P.fq

  Control-b d #to exit tmux

-----------------------------------------

**CHALLENGE**: Talk to me for details...

What Genus/Species did I sequence?
What tissue?

---------------------------------------------

==================================
Terminate your instance
==================================
