==========================
Running command-line BLAST
==========================

The goal of this tutorial is to run you through a demonstration of the
command line, which you may not have seen or used much before.

Prepare for this tutorial by working through
:doc:`amazon/start-up-an-ec2-instance`, but follow the instructions
to start up :doc:`amazon/starting-up-a-custom-ami` instead; use
AMI ami-7606d01e.

All of the commands below can and should be copy/pasted rather than
re-typed.

Note: on Windows using TeraTerm, you can select the commands in
the Web browser, then go to TeraTerm and click your right mouse
button to paste.  On Mac OS X using Terminal, you can select the
commands in the Web browser, use Command-C to copy, and then go
the terminal and use Command-V to paste.

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
           r-cran-gplots python-matplotlib sysstat

(make sure to hit enter after the paste -- sometimes the last line doesn't
paste completely.)

If you started up a custom operating system, then this should finish
quickly; if instead you started up Ubuntu 14.04 blank, then this will
take a minute or two.

Install BLAST
~~~~~~~~~~~~~

Here, we're using curl to download the BLAST distribution from NCBI;
then we're using 'tar' to unpack it into the current directory; and
then we're copying the program files into the directory
/usr/local/bin, where we can run them from anywhere.
::

   cd /root

   curl -O ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.26/blast-2.2.26-x64-linux.tar.gz
   tar xzf blast-2.2.26-x64-linux.tar.gz
   cp blast-2.2.26/bin/* /usr/local/bin
   cp -r blast-2.2.26/data /usr/local/blast-data

OK -- now you can run BLAST from anywhere!

Again, this is basically what "installing software" means -- it just
means copying around files so that they can be run, and (in some cases)
setting up resources so that the software knows where specific data
files are.

Running BLAST
~~~~~~~~~~~~~

Try typing::

   blastall

You'll get a long laundry list of output, with all sorts of options and
arguments.  Let's play with some of them.

First! We need some data.  Let's grab the mouse and zebrafish RefSeq
protein data sets from NCBI, and put them in /mnt, which is the
scratch disk space for Amazon machines
::

   cd /mnt

   curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.protein.faa.gz
   curl -O ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.protein.faa.gz

If you look at the files in the current directory, you should see both
files, along with a directory called lost+found which is for system
information::

   ls -l

should show you::

   drwx------ 2 root root   16384 2013-01-08 00:14 lost+found
   -rw-r--r-- 1 root root 9454271 2013-06-11 02:29 mouse.protein.faa.gz
   -rw-r--r-- 1 root root 8958096 2013-06-11 02:29 zebrafish.protein.faa.gz

Both of these files are FASTA protein files (that's what the .faa suggests)
that are compressed by gzip (that's what the .gz suggests).

Uncompress them
::

   gunzip *.faa.gz

and let's look at the first few sequences::

   head -11 mouse.protein.faa

These are protein sequences in FASTA format.  FASTA format is something
many of you have probably seen in one form or another -- it's pretty
ubiquitous.  It's just a text file, containing records; each record
starts with a line beginning with a '>', and then contains one or more
lines of sequence text.

Let's take those first two sequences and save them to a file.  We'll
do this using output redirection with '>', which says "take
all the output and put it into this file here."
::

   head -11 mouse.protein.faa > mm-first.fa

So now, for example, you can do 'cat mm-first.fa' to see the contents of
that file (or 'less mm-first.fa').

Now let's BLAST these two sequences against the entire zebrafish
protein data set. First, we need to tell BLAST that the zebrafish
sequences are (a) a database, and (b) a protein database.  That's done
by calling 'formatdb'
::

   formatdb -i zebrafish.protein.faa -o T -p T

Next, we call BLAST to do the search
::

   blastall -i mm-first.fa -d zebrafish.protein.faa -p blastp

This should run pretty quickly, but you're going to get a LOT of output!!
What's going on?  A few things --

 - if you BLAST a sequence against a large database, odds are it will turn
   up a lot of spurious matches.  By default, blastall uses an e-value cutoff
   of 10, which is very relaxed.

 - blastall also reports the first 100 matches, which is usually more than
   you want.

 - a lot of proteins also have trace similarity to other proteins!

For all of these reasons, generally you only want the first few BLAST
matches, and/or the ones with a "good" e-value.   We do that by adding
'-b 2 -v 2' (which says, report only two matches and alignments); and
by adding '-e 1e-6', which says, report only matches with an e-value
of 1e-6 or better
::

   blastall -i mm-first.fa -d zebrafish.protein.faa -p blastp -b 2 -v 2 -e 1e-6

Now you should get a lot less text!  (And indeed you do...) Let's put it an
output file, 'out.txt'
::

   blastall -i mm-first.fa -d zebrafish.protein.faa -p blastp -b 2 -v 2 -o out.txt

The contents of the output file should look exactly like the output before
you saved it into the file -- check it out::

   cat out.txt

Converting BLAST output into CSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose we wanted to do something with all this BLAST output.  Generally,
that's the case - you want to retrieve all matches, or do a reciprocal
BLAST, or something.

As with most programs that run on UNIX, the text output is in some
specific format.  If the program is popular enough, there will be one
or more parsers written for that format -- these are just utilities
written to help you retrieve whatever information you are interested
in from the output.

Let's conclude this tutorial by converting the BLAST output in out.txt
into a spreadsheet format, using a Python script.  (We're not doing this
just to confuse you; this is really how we do things around here.)

First, we need to get the script.  We'll do that using the 'git' program
::

   git clone https://github.com/ngs-docs/ngs-scripts.git /root/ngs-scripts

We'll discuss 'git' more later; for now, just think of it as a way
to get ahold of a particular set of files.  In this case, we've placed
the files in /root/ngs-scripts/, and you're looking to run the 
script blast/blast-to-csv.py using Python
::

   python /root/ngs-scripts/blast/blast-to-csv.py out.txt

This outputs a spread-sheet like list of names and e-values.  To save this
to a file, do::

   python /root/ngs-scripts/blast/blast-to-csv.py out.txt > /root/Dropbox/out.csv

The end file, 'out.csv', should soon be in your Dropbox on your local
computer.  If you have Excel installed, try double clicking on it.

----

And that's the kind of basic workflow we'll be teaching you:

1. Download program
2. Download data
3. Run program on data
4. Look at results

...but in many cases more complicated :).

----

Note that there's no limit on the number of sequences you BLAST, etc.
It's just sheer compute speed and disk space that you need to worry
about, and if you look at the files, it turns out they're not that big --
so it's mostly your time and energy.

This will also maybe help you understand why UNIX programs are so
powerful -- each program comes with several, or several dozen, little
command line "flags" (parameters), that help control how it does its
work; then the output is fed into another such program, etc.  The possibilities
are literally combinatorial.

----

We're running a Python program 'blast-to-csv.py' above -- if you're
interested in what the Python program does, take a look at the source
code:

   https://github.com/ngs-docs/ngs-scripts/blob/master/blast/blast-to-csv.py

Summing up
~~~~~~~~~~

Command-line BLAST lets you do BLAST searches of any sequences you have,
quickly and easily.  It's probably the single most useful skill a
biologist can learn if they're doing anything genomics-y ;).

Its main computational drawback is that it's not fast enough to deal
with some of the truly massive databases we now have, but that's
generally not a problem for individual users. That's because they just
run it and "walk away" until it's done!

The main practical issues you will confront in making use of BLAST:

 - getting your sequence(s) into the right place.
 - formatting the database.
 - configuring the BLAST parameters properly.
 - doing what you want after BLAST!

----

Other questions to ponder:

 - if we're using a pre-configured operating system, why did we have to
   install BLAST?
