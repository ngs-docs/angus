================================================
Trinity and Transcriptome Evaluation
================================================


**Launch a BIG maching** Maybe a c4.8xl that has 32 cores and 60GB RAM. Get 100Gb Hard Drive

Trinity: http://trinityrnaseq.github.io/

Transrate: http://hibberdlab.com/transrate/installation.html

RCorrector: https://github.com/mourisl/Rcorrector

BUSCO: http://busco.ezlab.org/


**Update Software**

::

    sudo apt-get update && sudo apt-get -y upgrade

**Install other software** Note that you can install a large amount of software from the Ubuntu "App Store" using a single command. Some of this software we will not use for this tutorial, but...

::

    sudo apt-get -y install build-essential tmux git gcc make g++ python-dev unzip \
                            default-jre libcurl4-openssl-dev zlib1g-dev python-pip


**Install Ruby**  Ruby is a computer language like Python or Perl.

::

    cd
    wget https://keybase.io/mpapis/key.asc
    gpg --import key.asc
    \curl -sSL https://get.rvm.io | bash -s stable --ruby
    source /home/ubuntu/.rvm/scripts/rvm



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

    curl -LO https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz
    tar -zxf transrate-1.0.3-linux-x86_64.tar.gz
    echo 'export PATH=$PATH:"$HOME/transrate-1.0.3-linux-x86_64"' >> ~/.profile
    source ~/.profile



**INSTALL BLAST**

::

    curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-x64-linux.tar.gz
    tar -zxf ncbi-blast-2.3.0+-x64-linux.tar.gz
    echo 'export PATH="$HOME/ncbi-blast-2.3.0+/bin:$PATH"' >> ~/.profile
    source ~/.profile

**INSTALL Augustus, BUSCO, Trinity, RCorrector, Skewer**

::

    brew install gcc augustus emboss Trinity --without-express Rcorrector Skewer busco --without-blast


**Download data**: For this lab, we'll be using
::

    #Open tumx window

    tmux new -s trinity

    mkdir $HOME/reads && cd /$HOME/reads/
    curl -LO https://s3.amazonaws.com/NYGC_August2015/raw_data/382-Kidney_ACTTGA_BC6PR5ANXX_L008_001.R1.fastq.gz
    curl -LO https://s3.amazonaws.com/NYGC_August2015/raw_data/382-Kidney_ACTTGA_BC6PR5ANXX_L008_001.R2.fastq.gz
    zcat 382-Kidney_ACTTGA_BC6PR5ANXX_L008_001.R1.fastq.gz | head -8000000 > read1.fq
    zcat 382-Kidney_ACTTGA_BC6PR5ANXX_L008_001.R2.fastq.gz | head -8000000 > read2.fq


    # control-b d to get out of tmux window.


    tmux attach -t trinity #to get back in tmux window.

**Correct Reads**

::

    run_rcorrector.pl -k 31 -t 30 \
    -1 $HOME/reads/read1.fq \
    -2 $HOME/reads/read2.fq



**Run Skewer**

::

    curl -LO https://s3.amazonaws.com/gen711/TruSeq3-PE.fa

    skewer -l 25 -m pe -o skewerQ2 --mean-quality 2 --end-quality 2 -t 30 \
    -x TruSeq3-PE.fa \
    $HOME/reads/read1.cor.fq \
    $HOME/reads/read2.cor.fq


**Run Trinity**

::

    mkdir $HOME/assembly && cd $HOME/assembly


    Trinity --seqType fq --max_memory 40G --left $HOME/reads/skewerQ2-trimmed-pair1.fastq \
    --right $HOME/reads/skewerQ2-trimmed-pair2.fastq --CPU 30


**Run BUSCO for assemblies**: There are Eukaryote, Metazoa, Arthropod, Vertebrate, Plant references for use with other genomes.

::


    mkdir $HOME/busco && cd $HOME/busco

    export AUGUSTUS_CONFIG_PATH=/home/ubuntu/.linuxbrew/Cellar/augustus/3.2.2_1/libexec/config/

    #Download busco database


    curl -LO http://busco.ezlab.org/files/vertebrata_buscos.tar.gz
    tar -zxf vertebrata_buscos.tar.gz

    busco -m trans -in $HOME/assembly/trinity_out_dir/Trinity.fasta \
    --cpu 30 -l vertebrata -o trin.assem

    less run*/short*

**Run Transrate**

::

    mkdir $HOME/transrate && cd $HOME/transrate
    transrate -a $HOME/assembly/trinity_out_dir/Trinity.fasta -t 30 \
    --left $HOME/reads/read1.cor.fq \
    --right $HOME/reads/read2.cor.fq


==================================
Terminate your instance
==================================
