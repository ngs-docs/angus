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

You will need an ``c3.xlarge`` instance (4 CPU, 8 GB RAM) and at least an 8 GB SSD root partition.

Install Kraken dependencies
===========================

::

   sudo apt-get -y update
   sudo apt-get -y install bioperl ruby build-essential curl git python-setuptools wget

Install Kraken
==============

::
  
  ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"

::

  export PATH="$HOME/.linuxbrew/bin:$PATH"
  brew tap homebrew/science
  brew install kraken
  
Install the MiniKraken database
===============================

::

  cd $HOME
  curl https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz | tar zxvf -
  export KRAKEN_DEFAULT_DB=$HOME/minikraken_20141208
  export KRAKEN_NUM_THREADS=$(getconf _NPROCESSORS_ONLN)

Get some reads
==============

These are paired-end Illumina reads taken from a skin microbiome study: <http://www.ebi.ac.uk/ena/data/view/SRR2423672>
::

  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR242/002/SRR2423672/SRR2423672_1.fastq.gz
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR242/002/SRR2423672/SRR2423672_2.fastq.gz

Run Kraken on the reads
=======================

::

  kraken --paired SRR2423672_1.fastq.gz SRR2423672_2.fastq.gz > kraken.out
  kraken-report kraken.out > kraken.tab
  
Examine the output
==================

::

  less kraken.tab

Which species are in this data set?

What is the dominant species on this person's skin?


More details
============

Kraken manual: http://ccb.jhu.edu/software/kraken/MANUAL.html

