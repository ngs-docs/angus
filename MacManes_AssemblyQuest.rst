================================================
Assembly Quest
================================================

I am going to give you some data and a genome (well, 500Mb of the genome), I want you to quality trim and map using the Trimming lecture and the Mapping Lecture.

** Launch an AMI. For this exercise, try a **c4.8xlarge** with a 200Gb EBS volume. Wow yes this is a BIG machine!!


You will need to install `abyss`
Installing BUSCO is a bit of a PITA... Here is the way, once you have brew installed.

::

    brew install emboss
    brew install busco

    curl -LO curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-x64-linux.tar.gz
    tar -zxf ncbi-blast-2.3.0+-x64-linux.tar.gz
    PATH=$HOME/ncbi-blast-2.3.0+/bin:$PATH

    export AUGUSTUS_CONFIG_PATH=/home/ubuntu/.linuxbrew/Cellar/augustus/3.2.2_1/libexec/config/

**Download the reads** They have already been trimmed for you :)..

::

    Paired-end data

    curl -LO https://s3.amazonaws.com/macmanes_share/ERR1328557_1.fastq.gz
    curl -LO https://s3.amazonaws.com/macmanes_share/ERR1328557_2.fastq.gz

    Mate-pair data

    curl -LO https://s3.amazonaws.com/macmanes_share/ERR1328558_1.fastq.gz
    curl -LO https://s3.amazonaws.com/macmanes_share/ERR1328558_2.fastq.gz


**Use ABySS to assemble reads** You have 32 threads to use! I suggest you use the all.

0. Pick a random kmer between 31 and 101

**Evaluate your assembly data using busco**  You'll need to download and un-compress the Fungal database before you run BUSCO. I'm giving you that code, but not the code for running BUSCO.. You'll have to do that yourself...

::

    cd $HOME
    curl -LO http://busco.ezlab.org/files/fungi_buscos.tar.gz
    tar -zxf fungi_buscos.tar.gz


**Add your data** about the assembly

https://docs.google.com/spreadsheets/d/1HzVZGn1FAo4GuxJXX3gv99Mi12h64mDMpD_DA3mgWp4/edit?usp=sharing


**Super bonus points to anybody that maps with a different assembler** e.g., SPAdes, AllPaths, Velvet, ...

**TERMINATE YOUR INSTANCE!!!**
