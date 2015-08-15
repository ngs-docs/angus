Automating Kmer Abundance Counts
================================

Adrienne Hoarfrost, aka *lady wizard* was given an assignment to automate this task. She expertly executed and provided her files for you in her github repository.

    **On Fri, Aug 14, 2015 at 07:19:19PM -0400, C. Titus Brown wrote:**
    Some of these files are good shotgun sequencing runs. 
    Others are not. 
    Which is which? 

These are the raw sequencing data files:

http://athyra.idyll.org/~t/transfer/aa.txt  

http://athyra.idyll.org/~t/transfer/bb.txt  

http://athyra.idyll.org/~t/transfer/cc.txt  

http://athyra.idyll.org/~t/transfer/dd.txt  

Link from Matt MacManes:
http://www.howtogeek.com/102420/geeks-versus-non-geeks-when-doing-repetitive-tasks-funny-chart/

Here are the steps to re-create Adrienne's task:

Start AWS instance
==================

Open an AWS instance and install the following packages:
t2.micro is fine

.. code:: bash

    sudo apt-get update &&
    sudo apt-get -y install screen 
        git curl python-dev unzip r-base-core r-cran-gplots 

Make a directory called *wand*, then clone Adrienne's github repository:
.. code:: 

    mkdir wand
    git clone https://github.com/ahoarfrost/wand.git
    cd wand

Her resulting script files are located here:

`_static/kmerhistograms.R`_
`_static/kmer.py`_
`_static/do.sh`_


Resulting pdf is here:
`_static/kmer10hists.pdf`_
