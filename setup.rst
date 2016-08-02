===============
Setup your AWS instance
===============

This lesson will show you how to set up you AWs instance for general use. Most lessons will use a similar process, though the specific software installed may vary.

**Update System Software**

::

    sudo apt-get update && sudo apt-get -y upgrade


**Install Basic System Software**

::

    sudo apt-get -y install cmake sparsehash valgrind libboost-atomic1.55-dev libibnetdisc-dev gsl-bin \
        libgsl0-dev libgsl0ldbl libboost1.55-all-dev libboost1.55-dbg subversion tmux git curl parallel \
        libncurses5-dev gcc make g++ python-dev unzip dh-autoreconf default-jre python-pip zlib1g-dev \
        mcl libhdf5-dev r-base pkg-config libpng12-dev libfreetype6-dev build-essential \
        libsm6 libxrender1 libfontconfig1 liburi-escape-xs-perl emboss liburi-perl


**Install Python Packages**

::
    pip install --user biopython numpy scipy sklearn

**Install Bioinformatics Packages via Brew**

::
    brew tap homebrew/science
    brew install samtools
    brew install hmmer
    brew install busco




    ===============
    TERMINATE your instance
    ===============
