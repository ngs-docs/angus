Transrate
###############

`Transrate <http://hibberdlab.com/transrate/index.html>` is a tool for assessing the quality of a de
novo transcriptome assembly.

We will use the assembly produced from the earlier Trinity tutorial.

Booting an Amazon AMI
~~~~~~~~~~~~~~~~~~~~~

Start up an Amazon computer (m1.large or m1.xlarge) using AMI
ami-7607d01e (see :doc:`amazon/start-up-an-ec2-instance` and
:doc:`amazon/starting-up-a-custom-ami`).

Log in `with Windows <amazon/log-in-with-ssh-win.html>`__ or
`from Mac OS X <amazon/log-in-with-ssh-mac.html>`__.

Updating the operating system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy and paste the following two commands.::
   apt-get update
   
   apt-get -y install screen git curl gcc make g++ python-dev unzip \
           default-jre pkg-config libncurses5-dev r-base-core \
           r-cran-gplots python-matplotlib sysstat

to update the computer with all the bundled software you'll need.

Install software
~~~~~~~~~~~~~~~~

First, we need to install `Transrate
<http://hibberdlab.com/transrate/installation.html>`__::

   cd /root
   mkdir transrate
   cd transrate
   wget -O transrate-1.0.1-linux-x86_64.tar.gz https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz
   tar xvzf transrate-1.0.1-linux-x86_64.tar.gz
   cp bwa /usr/local/bin

Overview
~~~~~~~~~~~~~~

Transrate analyses a transcriptome assembly in three key ways:

- by inspecting the contig sequences
- by mapping reads to the contigs and inspecting the alignments
- by aligning the contigs against proteins or transcripts from a related species and inspecting the alignments

We'll try out each one by one.
 
Level 1
-------

First level of inspection - generate basic statistics.::

	/root/transrate/transrate-1.0.1-linux-x86_64/transrate \
 	--assembly ../trinity/Trinity.fasta \
	--threads 1
 
The main output file is transrate_assemblies.csv (CSV = comma separated value)

Transfer via dropbox or scp and look at it in Excel.

Lets see what each of these mean by looking at the `documentation page <http://hibberdlab.com/transrate/metrics.html#contig-metrics>`_.

Level 2
-------
Lets map the reads and see how that looks. Transrate will
- map the provided reads to the assembly using `SNAP <http://snap.cs.berkeley.edu/>`_
- infer the most likely contig of origin for any multi-mapping reads with `Salmon <https://github.com/kingsfordgroup/sailfish/releases/tag/v0.3.0>`_
- inspect the resulting alignments with transrate-tools and use them to evaluate each contig in the assembly

Run transrate again, this time pointing to the input fastq files as well. Use the trimmed reads. ::

 /root/transrate/transrate-1.0.1-linux-x86_64/transrate \
 --assembly ../trinity/Trinity.fasta \
 --threads 1 \
 --left ../reads/allR1.fastq \
 --right ../reads/allR2.fastq

Lots of output files. Lets focus on the main metrics that were printed to the terminal 
as output. `Explanations are here <http://hibberdlab.com/transrate/metrics.html#read-mapping-metrics>`_

Level 3
-------
Lets try out the third capability of transrate - to compare our de novo assembly to a known proteome.

How to choose a reference species?

The ideal reference is one from a very closely related species that has a well 
annotated genome. If a well annotated genome is not available from a closely 
related species, then a set of proteins from a distantly related but well 
annotated genome is preferable to a closely related but poorly annotated one.

If your reference is from a species that is not very closely related, it is 
greatly preferable to use a set of proteins with only one protein 
per protein-coding gene. This is because most annotated genomes will have 
multiple isoforms for many genes, each producing a protein. The similar 
isoforms lead to confusing BLAST alignments and lower the quality of the 
results. For plant species, http://phytozome.jgi.doe.gov/pz/portal.html provides a 'single representative transcript per locus' set 
of proteins for every genome.

Get a proteome::

	curl -o 
	
Run transrate::

	/root/transrate/transrate-1.0.1-linux-x86_64/transrate \
	--assembly ../trinity/Trinity.fasta \
	--threads 1 \
	--reference /root/ref/genes.fa

