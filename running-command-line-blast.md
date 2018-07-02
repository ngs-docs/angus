# Running command-line BLAST

Learning objectives:

* gain hands-on exposure to the UNIX command line or shell,
  
* understand how data is turned into results by programs run at the command line.

* bridge BLAST data output to R/RStudio visualizations.

## Getting started

Boot an m1.medium instance on Jetstream and connect to your shell prompt.

Make sure you are starting in your home directory:

```
cd ~/
```

and let's make a new subdirectory to work in:

```
mkdir -p ~/blast
cd ~/blast
```

Now, install some software. We'll need NCBI BLAST for the below tutorial:

```
conda install -y blast
```

## Running BLAST

We need some data!  Let's grab the mouse and zebrafish RefSeq
protein data sets from NCBI, and put them in our home directory.

Now, we'll use `curl` to download the files from a Web site onto our
computer; note, these files originally came from the
[NCBI FTP site](ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot)

```
curl -o mouse.1.protein.faa.gz -L https://osf.io/v6j9x/download
curl -o mouse.2.protein.faa.gz -L https://osf.io/j2qxk/download

curl -o zebrafish.1.protein.faa.gz -L https://osf.io/68mgf/download
```

If you look at the files in the current directory, you should see four
files, along with a directory called lost+found which is for system
information. Running the command:

```
ls -l
```

should show you:

```
total 29908
-rw-rw-r-- 1 titus titus 12553742 Jun 29 08:41 mouse.1.protein.faa.gz
-rw-rw-r-- 1 titus titus  4074490 Jun 29 08:41 mouse.2.protein.faa.gz
-rw-rw-r-- 1 titus titus 13963093 Jun 29 08:42 zebrafish.1.protein.faa.gz
```

The three files you just downloaded are the last three on the list - the
`.faa.gz` files.

All three of the files are FASTA protein files (that's what the .faa
suggests) that are compressed with `gzip` (that's what the .gz means).

Uncompress them:

```
gunzip *.faa.gz
```

and let's look at the first few sequences in the file:

```
head mouse.1.protein.faa 
```

These are protein sequences in FASTA format.  FASTA format is something
many of you have probably seen in one form or another -- it's pretty
ubiquitous.  It's a text file, containing records; each record
starts with a line beginning with a '>', and then contains one or more
lines of sequence text.

Let's take those first two sequences and save them to a file.  We'll
do this using output redirection with '>', which says "take
all the output and put it into this file here."

```
head -11 mouse.1.protein.faa > mm-first.fa
```

So now, for example, you can do `cat mm-first.fa` to see the contents of
that file (or `less mm-first.fa`).

Now let's BLAST these two sequences against the entire zebrafish
protein data set. First, we need to tell BLAST that the zebrafish
sequences are (a) a database, and (b) a protein database.  That's done
by calling 'makeblastdb':

```
makeblastdb -in zebrafish.1.protein.faa -dbtype prot
```

Next, we call BLAST to do the search:

```
blastp -query mm-first.fa -db zebrafish.1.protein.faa
```

This should run pretty quickly, but you're going to get a lot of output!!
To save it to a file instead of watching it go past on the screen,
ask BLAST to save the output to a file that we'll name `mm-first.x.zebrafish.txt`:

```
blastp -query mm-first.fa -db zebrafish.1.protein.faa -out mm-first.x.zebrafish.txt
```

and then you can 'page' through this file at your leisure by typing:

```
less mm-first.x.zebrafish.txt
```

(Type spacebar to move down, and 'q' to get out of paging mode.)

-----

Let's do some more sequences (this one will take a little longer to run):

```
head -500 mouse.1.protein.faa > mm-second.fa
blastp -query mm-second.fa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.txt
```

will compare the first 83 sequences.  You can look at the output file with:

```
less mm-second.x.zebrafish.txt
```

(and again, type 'q' to get out of paging mode.)

Notes:

* you can copy/paste multiple commands at a time, and they will execute in order;

* why did it take longer to BLAST ``mm-second.fa`` than ``mm-first.fa``?

Things to mention and discuss:

* `blastp` options and -help.
* command line options, more generally - why so many?
* automation rocks!

----

Last, but not least, let's generate a more machine-readable version of that
last file --

```
blastp -query mm-second.fa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.tsv -outfmt 6
```

You can open the file with `less mm-second.x.zebrafish.tsv` to see how the file looks like.

See [this link](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) for a description of the possible BLAST output formats.

