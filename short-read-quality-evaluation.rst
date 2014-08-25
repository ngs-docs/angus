==========================
Short Read Quality Control
==========================

As useful as BLAST is, we really want to get into sequencing data,
right?  One of the first steps you must do with your data is
evaluate its quality and try to improve it.

.. highlights::
	**Summary**: a sequencing run may contain data of
	low reliability. It may also contain various contaminants and artificial
	sequence fragments. Some (but not all) of these problem can be corrected.

	**Caution**: Don't apply any type of correction without evaluating the results
	it produces.
	In general it is best to be conservative with QC. We are altering
	the data based on our expectations of what it should be like! The process may
	also introduce its own biases into the dataset.

Biostar QoD (questions of the day)
----------------------------------

QC is one of the most *mis-under-estimated* tasks of NGS data analysis. People
assume there is very little to it once they know how to run the tool.
The reality is a more complicated than that.

QC also happens to be
a pet peeve of mine (Istvan) as demonstrated below in the following
Biostar threads (and others):

 1. `FastQC quality control shootout <https://www.biostars.org/p/53528/>`__
 2. `So What Does The Sequence Duplication Rate Really Mean In A Fastqc Report <https://www.biostars.org/p/83842/>`__
 3. `Revisiting the FastQC read duplication report <https://www.biostars.org/p/107402/>`__

Quick Start
===========

The first part of this tutorial will run on your own computer. It assumes that you have
Java installed. Download both FastQC and two smaller datasets onto your system

  1. http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  2. http://apollo.huck.psu.edu/data/SRR519926_1.fastq.zip
  3. http://apollo.huck.psu.edu/data/SRR447649_1.fastq.zip

Run FastQC on each the files from the graphical interface. Let's discuss
the output in more detail.

FastQC at the command line
==========================

Before you can do that, though, you need to install a bunch o' software.

We will use a so called ``Makefile``, a simple text file
that can contains a series of commands that you could otherwise
type in yourself. There will be more information on
shell programming and automation later. For now think of a ``Makefile``
as a simple way to pick which commands you can execute yourself.
Let's get started. Install `make`::

	sudo apt-get install make -y

You can also investigate the `Makefile` yourself: https://github.com/ngs-docs/angus/blob/2014/files/Makefile-short-read-quality

This tutorial will download datasets. You may want to create a
directory (folder) that stores this data::

	mkdir qc
	cd qc

We assume that you are running the subsequent scripts from this folder.

.. important::

	Obtain the `Makefile` and save it onto your cloud system::

		# This is where you get the Makefile
		wget https://raw.githubusercontent.com/ngs-docs/angus/2014/files/Makefile-short-read-quality -O Makefile

You can investigate the file::

		# Look at the Makefile
		# more Makefile or pico Makefile

So we now have a `Makefile` and our system can
execute this `Makefile` via the `make` command.::

	make

Setup
-----

In our case you have to always specify which section of the
`Makefile` do you wish to execute. For example you can type::

	make setup

This will execute the parts of the Makefile that is listed below:

.. include:: files/Makefile-short-read-quality
	:code:
	:start-after: setup:
	:end-before: install:

Note that you could also just type in these commands yourself for
the same results. The Makefile just automates this.

Software Install
----------------

The next step is installing FastQC and Trimmomatic on your instance::

	make install

command will execute the following lines.

.. include:: files/Makefile-short-read-quality
	:code:
	:start-after: install:
	:end-before: download:

Where did stuff go? The downloaded code went into `~/src` the binaries are linked into `~/bin`
To test that everything works well type::

	fastqc -v

This will print the version number of `fastqc`

Download Data
-------------

Start at your EC2 prompt, then type ::

   make download

This will execute the following lines. Remember that you could also type
these in yourself.

.. include:: files/Makefile-short-read-quality
	:code:
	:start-after: download:
	:end-before: fastqc:

The FASTQ Format
----------------

In class explanation of the format. See a good description at
http://en.wikipedia.org/wiki/FASTQ_format

If you don't understand the format, you don't understand the basic premise of
the data generation!

Run a FastQC analysis on each dataset::

	make fastqc

would run the commands:

.. include:: files/Makefile-short-read-quality
	:code:
	:start-after: fastqc:
	:end-before: trim:

This command will generate an HTML file for each file. Copy these files to your dropbox and
look at them (a short walkthrough on what each plot means).

Alternatively you can generate the fastqc output directly to your Dropbox like so::

	fastqc *.fastq -o /mnt/Dropbox

Pattern Matching
----------------

We can also investigate what the files contain by matching::

	# Find start codons
	grep ATG SRR519926_1.fastq  --color=always | head

	# Find a subsequence
	grep AGATCGGAAG SRR519926_1.fastq  --color=always | head

Pattern matching via expressions is an extremely powerful concept. We'll revisit them later.

Trimming
--------

Note that there are vary large number of tools that perform quality/adapter trimming.

Now, run `Trimmomatic <http://www.usadellab.org/cms/index.php?page=trimmomatic>`__
to eliminate Illumina adapters from your sequences. First we need to find the adapter sequences::

	ln -s ~/src/Trimmomatic-0.32/adapters/TruSeq3-SE.fa

.. tip::
	You may also want to shorten the command line like so::

		alias trim='java -jar ~/src/Trimmomatic-0.32/trimmomatic-0.32.jar'

	You can now invoke the tool just by typing::

		trim

Among the (many) agonizing decisions that you will have to make is what
parameters to pick: how big should be my window be, how long should the reads be, what
should be the average quality be? What kinds of contaminants do I have. Run, rerun
and evaluate. Err on the side of caution.

Trim by quality alone:

.. include:: files/Makefile-short-read-quality
	:code:
	:start-after: trim:
	:end-before: trimPE:

Quality and clipping:

.. include:: files/Makefile-short-read-quality
	:code:
	:start-after: trimPE:
	:end-before: end:

Now a homework:

.. note::

	Read the manual for `Trimmomatic <http://www.usadellab.org/cms/index.php?page=trimmomatic>`__.
	Trim the reads in parallel for both readfiles in a sample.

.. Note::
	BTW: cool kids have pretty prompts, but you too can be cool, all you need to do is::

		echo "export PS1='\[\e]0;\w\a\]\n\[\e[32m\]\u@\h \[\e[33m\]\w\[\e[0m\]\n\$ '" >> ~/.bashrc

	Then relog. Don't ask why this works, it is one of those things that is best left undisturbed.



