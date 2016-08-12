================================================
Trinity and Transcriptome Evaluation
================================================

Trinity: http://trinityrnaseq.github.io/

Transrate: http://hibberdlab.com/transrate/installation.html



**Update Software**

::

    sudo apt-get update && sudo apt-get -y upgrade

**Install other software** Note that you can install a large amount of software from the Ubuntu "App Store" using a single command. Some of this software we will not use for this tutorial, but...

::

    sudo apt-get -y install build-essential tmux git gcc make g++ python-dev unzip \
                            default-jre libcurl4-openssl-dev zlib1g-dev python-pip ruby2.0

**Install LinuxBrew**

::

    sudo mkdir /home/linuxbrew
    sudo chown $USER:$USER /home/linuxbrew
    git clone https://github.com/Linuxbrew/brew.git /home/linuxbrew/.linuxbrew
    echo 'export PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"' >> ~/.profile
    echo 'export MANPATH="/home/linuxbrew/.linuxbrew/share/man:$MANPATH"' >> ~/.profile
    echo 'export INFOPATH="/home/linuxbrew/.linuxbrew/share/info:$INFOPATH"' >> ~/.profile
    source ~/.profile
    brew tap homebrew/science
    brew update
    brew doctor

**INSTALL TRANSRATE**

::

    gem install transrate

--------------


**INSTALL Augustus, BUSCO, Trinity, RCorrector, Skewer**

::

    brew install augustus emboss busco Trinity Rcorrector Skewer


**Download data**: For this lab, we'll be using
::

    mkdir /mnt/reads
    cd /mnt/reads/

    curl -LO https://www.dropbox.com/s/4o6eduzcw11gz53/subsamp.2.fq.gz

    curl -LO https://www.dropbox.com/s/i4wst01yz10i9x9/subsamp.1.fq.gz


**Correct Reads**

::

    Rcorrector ...



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


==================================
Terminate your instance
==================================
