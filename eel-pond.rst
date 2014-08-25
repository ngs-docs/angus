A complete de novo assembly and annotation protocol for mRNASeq
===============================================================

The goal of this tutorial is to run you through (part of) a real
mRNAseq analysis protocol, using a small data set that will complete
quickly.

Prepare for this tutorial by working through
:doc:`amazon/start-up-an-ec2-instance`, but follow the instructions to
start up :doc:`amazon/starting-up-a-custom-ami` instead; use AMI
ami-7607d01e.

Switching to root
~~~~~~~~~~~~~~~~~

Start by making sure you're the superuser, root::

   sudo bash

Updating the software on the machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy and paste the following two commands
::

   apt-get update
   apt-get -y install screen git curl gcc make g++ python-dev unzip \
           default-jre pkg-config libncurses5-dev r-base-core \
           r-cran-gplots python-matplotlib sysstat samtools python-pip

If you started up a custom operating system, then this should finish
quickly; if instead you started up Ubuntu 14.04 blank, then this will
take a minute or two.

Downloading the sample data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The mRNAseq protocol works with the data set that you put in '/data'.
Here, we will download a small data set (a subset of the data from
`this paper <http://www.ncbi.nlm.nih.gov/pubmed/23731568>`__, data
from embryonic Nematostella>`__), and put it in /data  ::

   mkdir /mnt/data
   ln -fs /mnt/data /data
   cd /data
   curl -O http://athyra.idyll.org/~t/mrnaseq-subset.tar
   tar xvf mrnaseq-subset.tar

Check it out::

   ls

You'll see a bunch of different files -- these are the kinds of files
you'll get from your sequencing facility.

Starting on the protocols
~~~~~~~~~~~~~~~~~~~~~~~~~

We're going to work with a special version of the protocols today, one
that we adapted specifically for this course.

**In general**, you should use the latest version, which will be at
https://khmer-protocols.readthedocs.org/.

For today, we'll be using http://khmer-protocols.readthedocs.org/en/ngs2014/
instead.

Work through the following:

#. `Quality trimming <http://khmer-protocols.readthedocs.org/en/ngs2014/mrnaseq/1-quality.html>`__

#. `Applying digital normalization <http://khmer-protocols.readthedocs.org/en/ngs2014/mrnaseq/2-diginorm.html>`__

#. `Running the actual assembly <http://khmer-protocols.readthedocs.org/en/ngs2014/mrnaseq/3-big-assembly.html>`__

#. `BLASTing your assembly <http://khmer-protocols.readthedocs.org/en/ngs2014/mrnaseq/installing-blastkit.html>`__

Actually using the BLAST Web server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To connect to your BLAST Web server, you need to enable inbound traffic
on your computer.  Briefly:

* go to your instance and look at what security group you're using
(should be 'launch-wizard-' something).  On the left panel, under
Network and Security, go into Security Groups. Select your
security group, and select Inbound, and Edit.   Click "Add rule", and
change "Custom TCP rule" to "http".  Then click "save".  Done!

You can try pasting this into your BLAST server::

   MDRSVNVIQCAAAPTRIQCEEINAKLMLGVGVFGLCMNIVLAVIMSFGAAHPHSHGMLSSVEFDHDVDYH
   SRDNHHGHSHLHHEHQHRDGCSHSHGNGGADMQRLECASPESEMMEEVVVETTSSNAESICSHERGSQSM
   NLRAAVLHVFGDCLQSLGVVLAACVIWAGNNSSVGVPSSAHSYYNLADPLLSVLFGVITVYTTLNLFKEV
   IVILLEQVPPAVEYTVARDALLSVEKVQAVDDLHIWAVGPGFSVLSAHLCTNGCATTSEANAVVEDAECR
   CRQLGIVHTTIQLKHAADVRNTGA
