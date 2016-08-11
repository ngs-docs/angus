==================================
Species identification with Kraken
==================================

In this tutorial you will:

1. Download and install Kraken
2. Run Kraken on a set of FASTQ reads
3. Interpret the results

The instructions below will work on a Ubuntu 14.04 Amazon instance.

Start an instance
=================

You will need an ``m4.xlarge`` instance (4 CPU, 16 GB) and make your root partition a 16 GB SSD (default is 8).

Install Kraken dependencies
===========================

::

   sudo apt-get -y update
   sudo apt-get install bioperl ruby build-essential curl git python-setuptools

Install Kraken
==============

::
  
  ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"
  export PATH="$USER/.linuxbrew/bin:$PATH"
  brew tap homebrew/science
  brew install kraken
  
Install the MiniKraken database
===============================

::

  wget https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz
  tar zxvf minikraken.tgz
  FIXME 

Get some reads
==============

::

  curl FIXME > R1.fq.gz
  curl FIXME > R2.fq.gz

Run Kraken on the reads
=======================

::

  kraken --threads 8 --db FIXME/minikraken --paired R1.fq.fz R2.fq.gz > kraken.out
  kraken-report kraken.out > kraken.tab
  
Examine the output
==================

::

  less kraken.tab

What species are in this data set?

Read more in the Kraken manual: http://ccb.jhu.edu/software/kraken/MANUAL.html

