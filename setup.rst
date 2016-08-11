===============
Setup your AWS instance
===============

This lesson will show you how to set up you AWS instance for general use. Most lessons will use a similar process, though the specific software installed may vary.

In General, there are a few differnet ways to install software. Installing from source, installing from apt-get, installing from Brew, pip, etc. You'll get to know each of these packages over the course of the next two weeks.

**Launch a c4.2xl instance**

**Update System Software** This command will check for updates, and install them. ``apt-get`` is like the OSX App store, for those of you with Macs.

::

    sudo apt-get update && sudo apt-get -y upgrade


**Install Basic System Software** This command will install various software on your AWS instance.

::

    sudo apt-get -y install build-essential cmake sparsehash valgrind libibnetdisc-dev gsl-bin \
        libgsl0-dev libgsl0ldbl  subversion tmux git curl parallel gcc make g++ zlib1g-dev \
        libncurses5-dev  python-dev unzip dh-autoreconf default-jre python-pip \
        mcl libhdf5-dev r-base pkg-config libpng12-dev libfreetype6-dev libsm6 \
        libxrender1 libfontconfig1 liburi-escape-xs-perl liburi-perl


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


**Install Python Packages** This will take some time...

::

    pip install --user biopython

**Install BLAST**

::

    curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.4.0+-x64-linux.tar.gz
    tar -zxf ncbi-blast-2.4.0+-x64-linux.tar.gz
    PATH=$PATH:$HOME/ncbi-blast-2.4.0+/bin



**Install Bioinformatics Packages via Brew** These are the packages that we will use to do *real* work!!! YAY!!!

::

    brew tap homebrew/science
    brew install hmmer
    brew install augustus
    brew install emboss
    brew install busco --without-blast


**How to tell if something is installed**

::

    which blastp
    which busco


** Most of the lessons** will start by installing these types of software. Practice makes perfect. Try terminating your instance, restarting, and reinstalling...

===============
TERMINATE your instance
===============
