================================================
So you want to get some sequencing data out of NCBI?
================================================

Requirements:  You'll need to start an Ubuntu EC2 instance and have root access.  The first part of this tutorial will be on your local computer and then we'll move onto the EC2 instance.  Also, this tutorial assumes that someone has talked to you about paths and you know how to change directories and execute a program on a file.  If you get an error that a program or file does not exist, make sure you are in the right path.

First, let's think about how these databases are structured.  I am going to cre\
ate a database for folks to deposit whole genome sequences.  What kind of infor\
mation am I going to store in this?  Many of you may be familiar with such a da\
tabase, hosted by the `NCBI <http://www.ncbi.nlm.nih.gov/>`_.  The scripts that complement this tutorial can be downloaded with the following::

    git clone https://github.com/adina/scripts-for-ngs.git

Let's come up with a list of things we'd like stored in	this database and discuss some	of the challenges involved in database creation, management, and access.

The challenge
-------------
So, you've been	given a	list of	genomes	and been asked to create a phylogenetic tree of these genomes.	 How big would this list be before you thought about hiring an undergraduate to download sequences?

Say the list is only 3 genomes::

   CP000962
   CP000967
   CP000975
   
The following are some ways with which I've used to grab genome sequences:

#. Use the web portal and look up each FASTA
#. Use the `FTP site <ftp://ftp.ncbi.nlm.nih.gov/refseq/>`_, find each genome, and download manually
#. Use the NCBI Web Services API to download the data

Among these, I'm going to assume many of you are familiar with the first two.  This tutorial then is going to go through an example of the first approach and then focus on using APIs.  

What is an API and how does it relate to NCBI?
----------------------------------------------

Here's some `answers <http://stackoverflow.com/questions/7440379/what-exactly-is-the-meaning-of-an-api>`_, among which my favorite is "an interface through which you access someone else's code or through which someone else's code accesses yours -- in effect the public methods and properties."

The NCBI has a whole toolkit which they call *Entrez Programming Utilities* or *eutils* for short.  You can read all about it in the `documentation <http://www.ncbi.nlm.nih.gov/books/NBK25501/>`_.  There are a lot of things you can do to interface with all things NCBI, including publications, etc., but I am going to focus today on downloading sequencing data.

To do this, you're going to be using one tool in *eutils*, called *efetch*.  There is a whole chapter devoted to `efetch <http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch>`_ -- when I first started doing this kind of work, this documentation always broke my heart.  Its easier for me to just show you how to use it.

You can use the NCBI efetch utility on your web browser (as is true of many APIs).  Open a web browser, and try the following URL to download the nucleotide genome for CB00962::

    http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000962&rettype=fasta&retmode=text

You'll note that a file downloaded on your computer.  Take a look at it.  You'll note that it looks a lot like what you would see if you had searched for CP00962 on the NCBI search page.  Check it out `here <http://www.ncbi.nlm.nih.gov/nuccore/CP000962>`_.

If I want to access other kinds of data associated with this genome.  For example, if I want the Genbank file as an output rather than a FASTA file, I would try the following command::

   http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000962&rettype=gb&retmode=text

Do you notice the difference in these two commands?  Let's breakdown the command here:

#.  <http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?>  This is command telling your computer program (or your browser) to talk to the NCBI API tool efetch.
#.  <db=nuccore>  This command tells the NCBI API that you'd like it to look in this particular database for some data.  Other databases that the NCBI has available can be found `here <http://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi>`_.
#.  <id=CP000962>  This command tells the NCBI API efetch the ID of the genome you want to find.
#.  <rettype=gb&retmode=text>  These two commands tells the NCBI how the data is returned.  You'll note that in the two examples above this command varied slightly.  In the first, we asked for only the FASTA sequence, while in the second, we asked for the Genbank file.  Here's some elusive documentation on where to find these `"return" objects <http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly>`_.  

Also, a useful command is also <version=1>.  There are different versions of sequences and some times that is useful.  For reproducibility, I try to specify versions in my queries, see these `comments <http://www.ncbi.nlm.nih.gov/Class/MLACourse/Modules/Format/exercises/qa_accession_vs_gi.html>`_.

.. Note:: 

   Notice the "&" that comes between each of these little commands, it is necessary and important.   

Automating with an API
----------------------

Ok, let's think of automating this sort of query.  

In the shell, you could run the same commands above with the addition of *curl* (a program to get information from remote sources) on your EC2 instance::

    curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000962&rettype=fasta&retmode=text"

You'll see it fly on to your screen.  Don't panic - you can save it to a file with the redirect command ">" and make it more useful.::

    curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000962&rettype=fasta&retmode=text" > CP000962.fa

You've now saved a query you've sent to the NCBI API to a file.  Now, you could now imagine writing a program where you made a list of IDs you want to download and put it in a for loop, *curling* each genome and saving it to a file.  The following is a `script <https://github.com/adina/scripts-for-ngs/blob/master/fetch-genomes.py>`_.  Thanks to Jordan Fish who gave me the original version of this script before I even knew how and made it easy to use.

This script has some nifty documentation that you can see by trying to execute the script.  To see the documentation for this script::

    python fetch-genomes.py

You'll see that you need to provide a list of IDs and a directory where you want to save the downloaded files.  What do you need to provide to this script?  The first thing is a file that contains a list of IDs (note that this is a required format, each ID on a new line) to fetch the data from NCBI.  Second, you need to name a directory where you want the program to put the files you fetch from the NCBI API.   Note that the lazy programmer who wrote this script requires you to identify a directory that does not currently exist.

So to use this script, two things are needed.  What are they?

To help out, I have provided a list of 50 IDs in a file called "interesting-genomes.txt."  How can you tell there are 50 genomes in this file?  But I don't want to just go wild and download all 50 at once.  
 
.. Note::
    
    You may want to run this on just a few of these IDs to begin with.  You can create a smaller list using the *head* command with the -n parameter in the shell.  For example, head -n 3 interesting-genomes.txt > 3genomes.txt. 
 
To run the script::
    python fetch-genomes.py 3genomes.txt 3genomes
    
Take a look at what happened, and when you're ready to try it for more files you could try the following.  Note that the directory you name in this command is arbitrary and defined by you::
    python fetch-genomes.py interesting-genomes.txt genbank-files

Let's take a look inside this script.  The meat of this script uses the following code::

    url_template = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=gb&retmode=text"

You'll see that the *id* here is a string character which is obtained from list of IDs contained in a separate file.  The rest of the script manages where the files are being placed and what they are named.  It also prints some output to the screen so you know its running.

Exercise - Downloading data
---------------------------

Try modifying the fetch_genomes.py script to download just the FASTA sequences of the genes.  

Running this script should allow you to download genomes to your heart's content.  But how do you grab specific genes from this data then?  Specifically, the challenge was to make a phylogenetic tree of sequences, so let's target the conserved bacterial gene, *16S ribosomal RNA gene*.

Comment on Genbank files
------------------------

Genbank files have a special structure to them.  You can look at it and figure it out for the most part, or read about it in detail `here <http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>`_.  To find out if your downloaded Genbank files contain 16S rRNA genes, I like to run the following command::

    grep 16S *gbk

This should look somewhat familiar from your shell lesson, but basically we're looking for anylines that contain the character "16S" in any Genbank file we've downloaded.  Note that you'll have to run this in the directory where you downloaded these files.

The structure of the Genbank file allows you to identify 16S genes.  For example, ::

         rRNA        9258..10759
                     /gene="rrs"
                     /locus_tag="CLK_3816"
                     /product="16S ribosomal RNA"
                     /db_xref="Pathema:CLK_3816"

You could write code to find text like 'rRNA' and '/product="16S ribosomal RNA"', grab the location of the gene, and then go to the FASTA file and grab these sequences.  To make things easy, there are existing packages to parse Genbank files.  I have the most experience with BioPython.  To begin with, let's just use BioPython to help us with our program.  

First, we'll have to install BioPython on your instance and they've made that pretty easy::

    sudo apt-get update
    sudo apt-get install python-biopython

Fan Yang (Iowa State University) and I wrote a script to extract 16S rRNA sequences from Genbank files, `here <https://github.com/adina/scripts-for-ngs/blob/master/parse-genbank.py>`_.  It basically searches for text strings in the Genbank structure that is appropriate for these particular genes.  You can read more about BioPython `here <http://biopython.org/DIST/docs/tutorial/Tutorial.html>`_ and its Genbank parser `here <http://biopython.org/DIST/docs/api/Bio.GenBank-module.html>`_.  In this script, we are looking for an "rRNA" feature and looking for specific text in its "/product" line.  If this is true, we go through the genome sequence and extract the coordinates of these genes, providing the specific gene sequence.

To run this script on the Genbank file for CP000962.  Note make sure you are in the right directory for both the program and the files::

    python parse-genbank.py genbank-files/CP000962.gbk > genbank-files/CP000962.gbk.16S.fa

The resulting output file contains all 16S rRNA genes from the given Genbank file.

To run this for multiple files, I use a shell for loop::

    for x in genbank-files/*; do python parse-genbank.py $x > $x.16S.fa; done

There are multiple ways to get this done -- but this is how I like to do it.  Now, you can figure out how you like to do it.

And there you have it, you can now pretty much automatically grab 16S rRNA genes from any number of genomes in NCBI databases.

Challenge:  
----------
Find your favorite gene, download a database of it from NCBI, and find matching sequences from a sequencing dataset.





