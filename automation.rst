Automating Kmer Abundance Counts
================================

Adrienne Hoarfrost, aka *Lady Wizard* was given an assignment to automate this task. She expertly executed and provided her files for you in her `github repository <https://github.com/ahoarfrost/wand>`__.

    **On Fri, Aug 14, 2015 at 07:19:19PM -0400, C. Titus Brown wrote:**
    Some of these files are good shotgun sequencing runs. 
    Others are not. 
    Which is which? 

These are the raw sequencing data files:

* `http://athyra.idyll.org/~t/transfer/aa.txt <http://athyra.idyll.org/~t/transfer/aa.txt>`  
* `http://athyra.idyll.org/~t/transfer/bb.txt <http://athyra.idyll.org/~t/transfer/bb.txt>`  
* `http://athyra.idyll.org/~t/transfer/cc.txt <http://athyra.idyll.org/~t/transfer/cc.txt>` 
* `http://athyra.idyll.org/~t/transfer/dd.txt <http://athyra.idyll.org/~t/transfer/dd.txt>` 

Cartoons, when to automate:

* http://www.howtogeek.com/102420/geeks-versus-non-geeks-when-doing-repetitive-tasks-funny-chart/  
* https://xkcd.com/1205/  
* https://xkcd.com/1319/  

Here are the steps to re-create Adrienne's task:

Start AWS instance
==================

Open an AWS instance and install the following packages (t2.micro is fine):


.. code:: bash

    sudo apt-get update &&
    sudo apt-get -y install git curl python-dev r-base-core r-cran-gplots 

Make a directory called *wand*, then clone *Lady Wizard's* github repository:

.. code:: 

    mkdir wand
    git clone https://github.com/ahoarfrost/wand.git
    cd wand

Her script files, with comments are located here:

* `kmerhistograms.R <_static/kmerhistograms.R>`__ 
* `kmer.py <_static/kmer.py>`__ 
* `do.sh <_static/do.sh>`__ 

After git clone, these files will be located in your directory /wand. Run do.sh

.. code:: 

      bash do.sh


After running, you should have this same resulting pdf:

* `kmer10hists.pdf <_static/kmer10hists.pdf>`__ 
