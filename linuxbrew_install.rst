======================================
Install Linuxbrew on your AWS instance
======================================

These commands will set up you vanilla Ubuntu 14.04 AWS instance for Linuxbrew.

::

    sudo apt-get update

    sudo apt-get -y upgrade

    sudo apt-get -y install build-essential cmake sparsehash valgrind libibnetdisc-dev gsl-bin \
        libgsl0-dev libgsl0ldbl subversion tmux git curl parallel gcc make g++ zlib1g-dev \
        libncurses5-dev python-dev unzip dh-autoreconf default-jre python-pip \
        mcl libhdf5-dev r-base pkg-config libpng12-dev libfreetype6-dev libsm6 \
        libxrender1 libfontconfig1 liburi-escape-xs-perl liburi-perl ruby \
	python-pandas python-scipy python-numpy

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"

    echo 'export PATH="$HOME/.linuxbrew/bin:$PATH"' >> ~/.profile
    echo 'export MANPATH="$HOME/.linuxbrew/share/man:$MANPATH"' >> ~/.profile
    echo 'export INFOPATH="$HOME/.linuxbrew/share/info:$INFOPATH"' >> ~/.profile
    source ~/.profile

    brew tap homebrew/science
    brew update
    brew doctor

   