RNASeq Transcript Mapping and Counting (BWA and HtSeq Flavor)
###############

The goal of this tutorial is to show you one of the ways to map RNASeq
reads to a transcriptome and to produce a file with counts of mapped reads for
each gene. 

We will be using `BWA <http://bio-bwa.sourceforge.net/>`__ for the mapping (previously used in the variant
calling example) and `HtSeq <http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html>`__ for the counting.

Booting an Amazon AMI
~~~~~~~~~~~~~~~~~~~~~
Start up an Amazon computer (m1.large or m1.xlarge) using AMI
ami-7607d01e (see :doc:`amazon/start-up-an-ec2-instance` and
:doc:`amazon/starting-up-a-custom-ami`).

Go back to the Amazon Console. Now select "snapshots" from the left had column.
Changed "Owned by me" drop down to "All Snapshots".  Search for "snap-028418ad" - 
(This is a snapshot with our test RNASeq Drosophila data from Chris)  The
description should be "Drosophila RNA-seq data".  Under "Actions" select "Create Volume", then ok.

Now on the left select "Volumes".  You should see an "in-use" volume - this is for your running instance, as well as an "available" volume - this is the one you just created from the snapshot from Chris and should have the snap-028418ad label. Select the available volume and from the drop down select
"Attach Volume". The white box pop up will appear - select in the empty instance box, your running 
instance should appear as an option.  Select it.  For the device, enter /dev/sdf.  Now attach.

Log in `with Windows <amazon/log-in-with-ssh-win.html>`__ or
`from Mac OS X <amazon/log-in-with-ssh-mac.html>`__.


Updating the operating system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Become root
::
   sudo bash

Copy and paste the following two commands
::

   apt-get update
   apt-get -y install screen git curl gcc make g++ python-dev unzip \
           default-jre pkg-config libncurses5-dev r-base-core \
           r-cran-gplots python-matplotlib sysstat

to update the computer with all the bundled software you'll need.

Mount the data volume.  (This is for Chris's data that we added as a snapshot).
::
   cd /root
   mkdir /mnt/ebs
   mount /dev/xvdf /mnt/ebs


Install software
~~~~~~~~~~~~~~~~

First, we need to install the `BWA aligner
<http://bio-bwa.sourceforge.net/>`__::

   cd /root
   wget -O bwa-0.7.10.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download

   tar xvfj bwa-0.7.10.tar.bz2
   cd bwa-0.7.10
   make

   cp bwa /usr/local/bin

We also need a new version of `samtools <http://samtools.sourceforge.net/>`__::

   cd /root
   curl -O -L http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
   tar xvfj samtools-0.1.19.tar.bz2
   cd samtools-0.1.19
   make
   cp samtools /usr/local/bin
   cp bcftools/bcftools /usr/local/bin
   cd misc/
   cp *.pl maq2sam-long maq2sam-short md5fa md5sum-lite wgsim /usr/local/bin/


