==========================
Running command-line BLAST
==========================

The goal of this tutorial is to run you through a demonstration of the
command line, which you may not have seen or used much before.

Start up an m1.xlarge Amazon EC2 instance.

All of the commands below can copy/pasted.

Updating the software on the machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy and paste the following commands
::

   sudo apt-get update && sudo apt-get -y install python ncbi-blast+

(make sure to hit enter after the paste -- sometimes the last line doesn't
paste completely.)

This updates the software list and installs the Python programming
language and NCBI BLAST+.

Running BLAST
~~~~~~~~~~~~~

First! We need some data.  Let's grab the mouse and zebrafish RefSeq
protein data sets from NCBI, and put them in /mnt, which is the
scratch disk space for Amazon machines
::

   sudo chmod a+rwxt /mnt
   cd /mnt

   curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.1.protein.faa.gz
   curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.2.protein.faa.gz
   curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.3.protein.faa.gz

   curl -O ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.1.protein.faa.gz


If you look at the files in the current directory, you should see four
files, along with a directory called lost+found which is for system
information::

   ls -l

should show you::

   total 28904
   drwx------ 2 root   root      16384 Nov  2  2014 lost+found
   -rw-rw-r-- 1 ubuntu ubuntu  8132346 Aug 10 21:44 mouse.1.protein.faa.gz
   -rw-rw-r-- 1 ubuntu ubuntu  8091255 Aug 10 21:44 mouse.2.protein.faa.gz
   -rw-rw-r-- 1 ubuntu ubuntu   565224 Aug 10 21:44 mouse.3.protein.faa.gz
   -rw-rw-r-- 1 ubuntu ubuntu 12735506 Aug 10 21:44 zebrafish.1.protein.faa.gz

All four of the files are FASTA protein files (that's what the .faa
suggests) that are compressed by gzip (that's what the .gz suggests).

Uncompress them
::

   gunzip *.faa.gz

and let's look at the first few sequences in the file::

   head mouse.1.protein.faa 

These are protein sequences in FASTA format.  FASTA format is something
many of you have probably seen in one form or another -- it's pretty
ubiquitous.  It's just a text file, containing records; each record
starts with a line beginning with a '>', and then contains one or more
lines of sequence text.

Let's take those first two sequences and save them to a file.  We'll
do this using output redirection with '>', which says "take
all the output and put it into this file here."
::

   head -11 mouse.1.protein.faa > mm-first.fa

So now, for example, you can do 'cat mm-first.fa' to see the contents of
that file (or 'less mm-first.fa').

Now let's BLAST these two sequences against the entire zebrafish
protein data set. First, we need to tell BLAST that the zebrafish
sequences are (a) a database, and (b) a protein database.  That's done
by calling 'makeblastdb'
::

   makeblastdb -in zebrafish.1.protein.faa -dbtype prot

Next, we call BLAST to do the search
::

   blastp -query mm-first.fa -db zebrafish.1.protein.faa

This should run pretty quickly, but you're going to get of output!!
To save it to a file instead of watching it go past on the screen,
do::

   blastp -query mm-first.fa -db zebrafish.1.protein.faa -out mm-first.x.zebrafish.txt

and then you can 'page' through this file at your leisure by typing::

   less mm-first.x.zebrafish.txt

(Type spacebar to move down, and 'q' to get out of paging mode.)

-----

Let's do some more sequences::

   head -500 mouse.1.protein.faa > mm-second.fa
   blastp -query mm-second.fa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.txt

will compare the first 83 sequences.  You can look at the output file with::

   less mm-second.x.zebrafish.txt

(and again, type 'q' to get out of paging mode.)

Note:

* you can execute multiple commands at a time;

* You should see a warning - ::

    Selenocysteine (U) at position 310 replaced by X

  what does this mean?

* why did it take longer to BLAST ``mm-second.fa`` than ``mm-first.fa``?

Things to mention and discuss:

* blastp options and -help.
* command line options, more generally - why??
* automation rocks!
