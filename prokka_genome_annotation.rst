================================================
Bacterial genome annotation using Prokka
================================================

(You should run this after running `Genome Assembly
<genome-assembly.html>`__.)

After you have de novo assembled your genome sequencing reads into contigs,
it is useful to know what genomic features are on those contigs. The process
of identifying and labelling those features is called genome annotation.

In this tutorial you will:

1. Download and install Prokka
2. Annotate a FASTA file of contigs
3. Search the resulting annotated genes with BLAST.

`Prokka <http://www.vicbioinformatics.com/software.prokka.shtml>`__ is a tool that facilitates the fast annotation of prokaryotic genomes.

The goals of this tutorial are to:

* Install Prokka
* Use Prokka to annotate our genomes
* Take a brief look at some downstream uses of the annotations.

Installing Prokka
=================

Download and extract the latest version of prokka:
::
    cd ~/
    wget http://www.vicbioinformatics.com/prokka-1.12.tar.gz
    tar -xvzf prokka-1.12.tar.gz

We also will need some dependencies such as bioperl:
::
    sudo apt-get -y install bioperl libdatetime-perl libxml-simple-perl \
        libdigest-md5-perl python ncbi-blast+ fastqc

and then we have to install some Perl libraries, too::

    sudo bash
    export PERL_MM_USE_DEFAULT=1
    export PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"
    perl -MCPAN -e 'install "XML::Simple"'
    exit

Now, you should be able to add Prokka to your ``$PATH`` and set up the index for the sequence database:
::
    export PATH=$PATH:$HOME/prokka-1.12/bin
    echo 'export PATH=$PATH:$HOME/prokka-1.12/bin' >> ~/.bashrc
    prokka --setupdb

Prokka should be good to go now-- you can check to make sure that all is well by typing ``prokka``. This should print the help screen with all available options.

Running Prokka
==============

Make a new directory for the annotation:
::
    cd ~/
    mkdir annotation
    cd annotation

Link the metagenome assembly file into this directory:
::
    ln -fs ~/work/ecoli-assembly.fa

Now it is time to run Prokka! There are tons of different ways to specialize the running of Prokka. We are going to keep it simple for now, though. It will take a little bit to run.
::
    prokka ecoli-assembly.fa --outdir prokka_annotation --prefix myecoli

This will generate a new folder called ``prokka_annotation`` in which will be a series of files, which are detailed `here <https://github.com/tseemann/prokka/blob/master/README.md#output-files>`__.

The output of prokka is pretty interesting to watch scroll by - it does a
lot of different things, including finding tRNAs.

One of the main files output by prokka for annotation is a GFF file,
which stands for "gene feature format".  Take a quick look at the
annotations --
::
   grep -v ^## prokka_annotation/myecoli.gff | head

and you can see what kinds of things Prokka outputs.  (Note the
``putative deoxyribonuclease Rhsc`` halfway down the page.)

Interestingly, we have now "closed the loop" on variant calling --
if you recall `from the variant calling lesson <variant-calling.html#look-at-the-vcf-file-with-bedtools>`__, in order to determine where the variants were
in the genome, we needed a GFF file.  And, before that, in order to map
the reads to a reference, we need a reference!  Well, we assembled a reference
in Genome Assembly, and we just constructed a GFF annotation file using
Prokka, so now if you had a *different* set of reads from a different
(mutated or evolved) genome, you could map those reads to this set and
call variants.

-----

There are several other things you can do with these annotation files.  It's
close to being ready to upload to NCBI as a provisional genome
anonotation for this microbe, believe it or not!  You can also
download the GFF file and use it in a genome browser like
[Artemis](http://www.sanger.ac.uk/science/tools/artemis), and we'd be
happy to show you how to do that (but it involves installing software on
your laptop, so we avoid doing that in class).

Searching the annotated genes
=============================

Another thing you can do is BLAST.  If you recall from `the very first
thing we did <running-command-line-blast.html>`__, BLAST lets you
search gene sequences against gene sequences.  Helpfully, prokka
outputs all of the predicted protein coding sequences as
``prokka_annotation/myecoli.faa``::

   head prokka_annotation/myecoli.faa

Let's go grab ALL of
the E. coli genes from NCBI's `assembly database <https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2>`__::

  curl -L -o ncbi-ecoli.faa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_protein.faa.gz

and make a BLAST database of it::

  gunzip ncbi-ecoli.faa.gz
  makeblastdb -in ncbi-ecoli.faa -dbtype prot

and now search it with the genes predicted by prokka::

  blastp -query prokka_annotation/myecoli.faa -db ncbi-ecoli.faa | less

page down through the matches using spacebar ('q' to quit) and look to
see if you find some good ones.  (The first few seem spurious which is
odd, since we're searching E. coli against E. coli! But the matches are not
significant - see the e-values.)
  
References
===========

* http://www.vicbioinformatics.com/software.prokka.shtml
* https://www.ncbi.nlm.nih.gov/pubmed/24642063
* https://github.com/tseemann/prokka/blob/master/README.md
