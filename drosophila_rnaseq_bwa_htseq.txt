RNA-seq: mapping to a reference genome with BWA and counting with HTSeq
-----------------------------------------------------------------------

The goal of this tutorial is to show you one of the ways to map RNASeq
reads to a transcriptome and to produce a file with counts of mapped
reads for each gene. This is an alternative approach to mapping to the
reference genome, and by using the same dataset as the previous lesson
(see `drosophila\_rnaseq1 <drosophila_rnaseq1.txt>`__, we can see
the differences between the two approaches.

We will again be using `BWA <http://bio-bwa.sourceforge.net/>`__ for the
mapping (previously used in the variant calling example) and
`HTSeq <http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html>`__
for the counting.

Booting an Amazon AMI
---------------------

Start up an Amazon computer (m1.large or m1.xlarge) using AMI
ami-7607d01e (see :doc:`amazon/start-up-an-ec2-instance` and
:doc:`amazon/starting-up-a-custom-ami`.

Go back to the Amazon Console.

-  Select "snapshots" from the left side column.
-  Changed "Owned by me" drop down at the top to "All Snapshots"
-  Search for "snap-028418ad" - (This is a snapshot with our test RNASeq
   Drosophila data from Chris) The description should be "Drosophila
   RNA-seq data".
-  Under "Actions" select "Create Volume", then ok.

*Make sure to create your EC2 instance and your EBS volume in the same
availability zone, for this course we are using N. Virginia.*

-  Select "Volumes" from the left side column
-  You should see an "in-use" volume - this is for your running
   instance. You will also see an "available" volume - this is the one
   you just created from the snapshot from Chris and should have the
   snap-028418ad label. Select the available volume
-  From the drop down select "Attach Volume".
-  A white box pop up will appear - click in the empty instance box,
   your running instance should appear as an option. Select it. 
-  For the device, enter /dev/sdf.
-  Attach.

Log in `with Windows <amazon/log-in-with-ssh-win.html>`__ or
`from Mac OS X <amazon/log-in-with-ssh-mac.html>`__.

Updating the operating system
-----------------------------

Become root:
::
   sudo bash

Copy and paste the following two commands to update the computer with
all the bundled software you'll need.
::
      apt-get update
      apt-get -y install screen git curl gcc make g++ python-dev unzip \
              pkg-config libncurses5-dev r-base-core \
              r-cran-gplots python-matplotlib sysstat git python-pip
      pip install pysam

Mount the data volume. (This is the volume we created earlier from
Chris's snapshot - this is where our data will be found):
::
      cd /root
      mkdir /mnt/ebs
      mount /dev/xvdf /mnt/ebs

Install software
----------------

First, we need to install the `BWA aligner <http://bio-bwa.sourceforge.net/>`__
::
      cd /root
      wget -O bwa-0.7.10.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download
      tar xvfj bwa-0.7.10.tar.bz2
      cd bwa-0.7.10
      make
      cp bwa /usr/local/bin

We also need a new version of `samtools <http://samtools.sourceforge.net/>`__:
::
      cd /root
      curl -O -L http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
      tar xvfj samtools-0.1.19.tar.bz2
      cd samtools-0.1.19
      make
      cp samtools /usr/local/bin
      cp bcftools/bcftools /usr/local/bin
      cd misc/
      cp *.pl maq2sam-long maq2sam-short md5fa md5sum-lite wgsim /usr/local/bin/

Create a working directory to hold some more software that we're going
to install
::
       cd /mnt/ebs
       mkdir tools
       cd tools
       
Download and install HTSeq
::
       curl -O https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1.tar.gz
       tar -xzvf HTSeq-0.6.1.tar.gz
       cd HTSeq-0.6.1/
       python setup.py build
       python setup.py install
       chmod u+x ./scripts/htseq-count

We are also going to get a project, chado-test, from Scott Cain's git
hub account that will allow us to use a convenient file format
conversion script.
::
      cd /mnt/ebs/tools
      git clone https://github.com/scottcain/chado_test.git

For that we will need bioperl installed
::
       cpan

Answer yes until you get a prompt that looks like
::
       cpan[1]>

And type
::
       install Bio::Perl

When it asks "Do you want to run tests that require connection to
servers across the internet", answer no. The final line when finished
should be:
::
      ./Build install  -- OK

Now exit the CPAN shell
::
      exit

Preparing the reference
-----------------------

Next, we are going to work with our reference transcriptome. Drosophila
has a reference genome, but for this adventure, we are going to pretend
that it doesn't. Instead we are going to use the Trinity assembly as our
reference - Chris has provided this file, named Trinity\_all\_X.fasta.
Notice the fasta format; each line beginning with a > is a new sequence,
followed by another line (or multiple lines) containing the sequence
itself. If we want to count how many transcripts are in the file, we can
just count the number of lines that begin with >
::
      cd /mnt/ebs/trinity_output
      grep '>' Trinity_all_X.fasta | wc -l

You should see 8260. Now lets use bwa to index the file, this enables
the file to be used a reference for mapping:
::
      bwa index Trinity_all_X.fasta

To generate count files, we will use HTSeq. But HTSeq is expecting a
genome annotation file, which we don't have (since we're using the
transcriptome). So we have to do some data massaging. We will will
create an annotation file that says that the entire length of each
"scaffold" is in fact a coding region.
::
     cd /mnt/ebs/rnaseq_mapping2
     /mnt/ebs/tools/chado_test/chado/bin/gmod_fasta2gff3.pl \
     --fasta_dir /mnt/ebs/trinity_output/Trinity_all_X.fasta \
     --gfffilename Trinity_all_X.gff3 \
     --type CDS \
     --nosequence

Now you should have a file named Trinity\_all\_X.gff3 in your current
directory.

Mapping
-------

Lets check out the reads to be mapped
::
       cd /mnt/ebs/drosophila_reads
       ls -lh

Don't forget that with your reads, you'll want to take care of the usual
QC steps before you actually begin your mapping. The drosophila\_reads
directory contains raw reads; the trimmed\_x directory contains reads
that have already been cleaned using Trimmomatic. We'll use these for
the remainder of the tutorial, but you may want to try running it with
the raw reads for comparison.

We've got 12 sets of data, each with two files (R1 and R2). Let's run
bwa on the first pair to map our paired-end sequence reads to the
transcriptome. To make our code a little more readable and flexible,
we’ll use shell variables in place of the actual file names. In this
case, let’s first specify what the values of those variables should be:
::
      reference=/mnt/ebs/trinity_output/Trinity_all_X.fasta
      sample=HYB_sdE3_rep1

Now we can use these variable names in our mapping commands. The
advantage here is that we can just change the variables later on if we
want to apply the same pipeline to a new set of samples:
::
      cd /mnt/ebs
      mkdir rnaseq_mapping2
      cd rnaseq_mapping2
      bwa mem ${reference} /mnt/ebs/trimmed_x/${sample}_1_pe /mnt/ebs/trimmed_x/${sample}_2_pe > ${sample}.sam

The output is a file named HYB\_sdE3\_rep1\_2.sam in the current working
directory. This file contains all of the information about where each
read hits on the reference. Next, we want to use SAMTools to convert it
to a BAM, and then sort and index it:
::
      samtools view -Sb ${sample}.sam > ${sample}.unsorted.bam
      samtools sort ${sample}.unsorted.bam ${sample}
      samtools index ${sample}.bam

Now we can generate a counts file with the HTSeq-count script:
::
      htseq-count --format=bam --stranded=no --type=CDS --order=pos --idattr=Name ${sample}.bam Trinity_all_X.gff3 > ${sample}_htseq_counts.txt 

Optional - Script these steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since we have a lot of files to map, it would take a long time to
re-write the mapping commands for each one. And with so many parameters,
we might make a mistake or typo. It's usually safer to use a simple
shell script with shell variables to be sure that we do the exact same
thing to each file. Using well-named shell variables also makes our code
a little bit more readable. Open a file named map\_and\_count.sh and
paste in the following code:
::
       #Create an array to hold the names of all our samples
       #Later, we can then cycle through each sample using a simple foor loop
       samples[1]=ORE_wt_rep1
       samples[2]=ORE_wt_rep2
       samples[3]=ORE_sdE3_rep1
       samples[4]=ORE_sdE3_rep2
       samples[5]=SAM_wt_rep1
       samples[6]=SAM_wt_rep2
       samples[7]=SAM_sdE3_rep1
       samples[8]=SAM_sdE3_rep2
       samples[9]=HYB_wt_rep1
       samples[10]=HYB_wt_rep2
       samples[11]=HYB_sdE3_rep1
       samples[12]=HYB_sdE3_rep2
       
       #Create a shell variable to store the location of our reference genome 
       reference=/mnt/ebs/trinity_output/Trinity_all_X.fasta
       
       #Make sure we are in the right directory
       #Let's store all of our mapping results in /mnt/ebs/rnaseq_mapping2/ to make sure we stay organized
       #If this directory already exists, thats ok, but files might get overwritten
       cd /mnt/ebs
       mkdir rnaseq_mapping2
       cd rnaseq_mapping2
       
       #Now we can actually do the mapping and counting
       for i in 1 2 3 4 5 6 7 8 9 10 11 12
       do
           sample=${samples[${i}]}
           #Map the reads
           bwa mem ${reference} /mnt/ebs/trimmed_x/${sample}_1_pe /mnt/ebs/trimmed_x/${sample}_2_pe  > ${sample}.sam
           samtools view -Sb ${sample}.sam > ${sample}.unsorted.bam
           samtools sort ${sample}.unsorted.bam ${sample}
           samtools index ${sample}.bam
           htseq-count --format=bam --stranded=no --type=CDS --order=pos --idattr=Name ${sample}.bam Trinity_all_X.gff3 > ${sample}_htseq_counts.txt
       done

To run this script, change the permissions and run:
::
      chmod u+x ./map_and_count.sh
      ./map_and_count.sh

We now have count files for each sample. Take a look at one of the count
files using less. You'll notice there are a lot of zeros, but that's
partially because we've already filtered the dataset for you to include
only reads that map to the X chromosome.

You can also visualize these read mapping using tview
:doc:`variant`.
