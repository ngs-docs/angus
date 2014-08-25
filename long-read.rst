PacBio Tutorial
---------------

Launch a generic AMI (m3.2xlarge), update and install basic software.
You can use the generic ami-864d84ee or any other Ubuntu machine.

::

    #update stuff
    sudo apt-get update

    #install basic software
    sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
    default-jre pkg-config

    #install Perl modules required by PBcR, paste in to terminal one at a time.. 
    #Will be a couple of prompts (answer YES to both)
    sudo cpan App::cpanminus
    sudo cpanm Statistics::Descriptive

    #Install wgs-assembler
    wget http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.2beta/wgs-8.2beta-Linux_amd64.tar.bz2
    tar -jxf wgs-8.2beta-Linux_amd64.tar.bz2

    #add wgs to $PATH
    PATH=$PATH:$HOME/wgs-8.2beta/Linux-amd64/bin/

Download sample Lambda phage dataset. We are using this only because it
is very small and can be assembled quickly and with limited hardware
requirements. For a more challenging test (read: expert with a big
computer) try one of publicly available PacBio datasets here:
https://github.com/PacificBiosciences/DevNet/wiki/Datasets

::

    #make sure you have the appropriate permissions to read and write.
    sudo chown -R ubuntu:ubuntu /mnt
    mkdir /mnt/data
    cd /mnt/data

    #Download the sample data
    wget http://www.cbcb.umd.edu/software/PBcR/data/sampleData.tar.gz
    tar -zxf sampleData.tar.gz
    cd sampleData/

Convert fastA to faux-fastQ

::

    #This is really old PacBio data, provided in fastA format. Look at the reads - note that they are not actually as long as I just told you they should be. The PacBio tech has improved massively over the past few years. 

    java -jar convertFastaAndQualToFastq.jar \
    pacbio.filtered_subreads.fasta > pacbio.filtered_subreads.fastq

Run the assembly, using wgs, after error-correcting the reads. You could
do the error correction separately, but no need to, here, for our
purposes.

::

    PBcR -length 500 -partitions 200 -l lambda -s pacbio.spec \
    -fastq pacbio.filtered_subreads.fastq genomeSize=50000

Look at the output. The phage genome has been assembled into 2 contigs
(meh). Try a larger dataset for a more difficult (and rewarding
challenge)