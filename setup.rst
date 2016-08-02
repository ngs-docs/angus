===============
Setup your AWS instance
===============

This lesson will show you how to set up you AWs instance for general use. Most lessons will use a similar process, though the specific software installed may vary.

**Update System Software**

::

    sudo apt-get update && sudo apt-get -y upgrade


**Install Basic System Software**

::

    sudo apt-get -y install build-essential cmake sparsehash valgrind libibnetdisc-dev gsl-bin \
        libgsl0-dev libgsl0ldbl  subversion tmux git curl parallel gcc make g++ zlib1g-dev \
        libncurses5-dev  python-dev unzip dh-autoreconf default-jre python-pip \
        mcl libhdf5-dev r-base pkg-config libpng12-dev libfreetype6-dev libsm6 \
        libxrender1 libfontconfig1 liburi-escape-xs-perl liburi-perl


**Install Ruby**

::

    cd
    wget https://keybase.io/mpapis/key.asc
    gpg --import key.asc
    \curl -sSL https://get.rvm.io | bash -s stable --ruby
    source /home/ubuntu/.rvm/scripts/rvm

**Install Brew**

::

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"

    # press return



    echo 'export PATH="/home/ubuntu/.linuxbrew/bin:$PATH"' >>~/.profile
    echo 'export MANPATH="/home/ubuntu/.linuxbrew/share/man:$MANPATH"' >>~/.profile
    echo 'export INFOPATH="/home/ubuntu/.linuxbrew/share/info:$INFOPATH"' >>~/.profile
    source ~/.profile


**Install Python Packages**

::

    pip install --user biopython numpy scipy sklearn

**Install Bioinformatics Packages via Brew**

::

    brew tap homebrew/science
    brew install samtools
    brew install hmmer
    brew install augustus
    brew install emboss
    brew install busco




===============
TERMINATE your instance
===============
