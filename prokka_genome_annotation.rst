================================================
Bacterial genome annotation using Prokka
================================================

After you have de novo assembled your genome sequencing reads into contigs,
it is useful to know what genomic features are on those contigs. The process
of identifying and labelling those features is called genome annotation.

In this tutorial you will:

1. Download and install Prokka
2. Annotate a FASTA file of contigs
3. Visualize the annotation using Artemis

Prokka is a tool that facilitates the fast annotation of prokaryotic genomes.

The goals of this tutorial are to:

*  Install Prokka
*  Use Prokka to annotate our genomes

Installing Prokka
=================

Download and extract the latest version of prokka:
::
    cd ~/
    wget http://www.vicbioinformatics.com/prokka-1.11.tar.gz
    tar -xvzf prokka-1.11.tar.gz

We also will need some dependencies such as bioperl:
::
    sudo apt-get install bioperl libdatetime-perl libxml-simple-perl libdigest-md5-perl
    sudo perl -MCPAN -e shell
    sudo perl -MCPAN -e 'install "XML::Simple"'

Now, you should be able to add Prokka to your ``$PATH`` and set up the index for the sequence database:
::
    export PATH=$PATH:$HOME/prokka-1.11/bin
    prokka --setupdb

Prokka should be good to go now-- you can check to make sure that all is well by typing ``prokka``. This should print the help screen with all available options.

Running Prokka
==============

Make a new directory for the annotation:
::
    cd 
    mkdir annotation
    cd annotation

Link the metagenome assembly file into this directory:
::
    ln -fs ~/work/ecoli-assembly.fa

Now it is time to run Prokka! There are tons of different ways to specialize the running of Prokka. We are going to keep it simple for now, though. It will take a little bit to run.
::
    prokka ecoli-assembly.fa --outdir prokka_annotation --prefix metagG

This will generate a new folder called ``prokka_annotation`` in which will be a series of files, which are detailed `here <https://github.com/tseemann/prokka/blob/master/README.md#output-files>`__.


References
===========

* http://www.vicbioinformatics.com/software.prokka.shtml
* https://www.ncbi.nlm.nih.gov/pubmed/24642063
* https://github.com/tseemann/prokka/blob/master/README.md
