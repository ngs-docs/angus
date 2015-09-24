.. _aligncount

Align and count gene features
=============================

Several aligners exist for mapping reads back to a reference genome (Bowtie, Bowtie2, BWA, SOAP, etc.).
Each aligner was developed for a specific use case in mind and one should be aware of these
differences. For example, Bowtie was optimized for identification of splice junctions. We
don't have to worry about that for microbes, but the aligner itself works just fine even in
gene dense genomes like bacteria. Another caveat is there are *two* flavors of Bowtie: Bowtie1
and Bowtie2. Bowtie1 seems to perform better for 50 bp reads while 100+ bp reads tend to
be faster using Bowtie2.

But what about the other aligners? Well, to start, look and see what other people are using
in the literature, talk to your colleagues, or your friendly neighborhood bioinformatician.

The take home messages are: 
    
	* Many of these tools have great documentation, so it's worth having a look to see if they explicitly state the use case in mind and whether your experiment/data fits.
 
 	* Check the literature
 	
 	* Try them out and compare the results
 	
 The next component in all of this is whether you are doing reference or reference free assembly.
 More on reference free assembly later. For our purposes, we have a reference.
 
 My favorite place to download a reference genome and annotation files are from `Ensembl <http://bacteria.ensembl.org/info/website/ftp/index.html>`__.

Align the trimmed reads with Bowtie
-----------------------------------

Download the reference genome::

	cd ../../Alignment
	wget http://s3-us-west-1.amazonaws.com/microgenomicstranscriptomics/MtbCDC1551.fa
	
Build the index of the genome for bowtie to use to perform the mapping/alignment process::

	bowtie-build MtbCDC1551.fa MtbCDC1551

The counting transcript abundance step won't read in compressed files (files that end with .gz),
so we will unzip them now and hand them to bowtie (the aligner) unzipped.

We could do this in individual lines by typing something like::

    gunzip -c gly7a.fq.gz > gly7a.fq
    
But this would be rather tedious to do this for every file. Then we have to move the files
to our Alignment folder by typing something like::

    mv gly7a.fq ~/TranscriptomicsWorkshop/Alignment
    
But again, this would be rather tedious to do... So what we will do is script it and combine
the decompressing with moving the file into a single command with &&. This is *really* useful
in a Linux type system to string two commands together. What it does is say "if the command
on the left of the && works, then do the command on the right of the &&. Finally, the . in
the move (mv) command means "the folder you are in right now".

So, let's do it!

Unzip the files and move the files to the Alignment folder (which we should already be in). Copy and paste the lines below into the terminal::

	#decompress the files
	for compfiles in ../QC/Trimmomatic/*.fq.gz
	do
		FILENAME=${compfiles%.*.*}
		gunzip -c $FILENAME.fq.gz > $FILENAME.fq && mv $FILENAME.fq .
	done
	
Now, we have been adding prefixes to files so we have an idea based on the name, how the file
has been modified during the analysis. We want to strip off this prefix before moving on so that
we can add a new one. To do this, we can use the move (mv) command again to rename a file.

The syntax looks something like this as an example::

    mv trimmedgly7a.fq gly7a.fq
    
So essentially, you have your command, mv and then the file you want to rename and then the new
name of the file.

But to do this for every file is a pain and error prone. Script it!::

	#rename the files by stripping off the trimmed prefix
	for renametrim in *.fq
	do
		LESSPREFIX=`basename ${renametrim//trimmed/}`
		mv $renametrim $LESSPREFIX
	done

Whew, that was a lot of work to get to the alignment step. We are going to use Bowtie to do
the alignment of the reads. Bowtie has **a lot** of options to tailor the alignment procedure.
Here we only specifiy the -S option to let Bowtie know that we want our output in .sam format. The
general syntax for running bowtie is as follows::

    bowtie -S [index file prefix from the bowtie-build command] [your input file] > [your output file that ends in .sam]

We have another new thing that looks like the greater than (>) sign. This says "take the output
from the command on the left and put it in a new file on the right".

So let's align the reads::

	#align the reads
	for alignment in *.fq
	do	
		FILENAME=`basename ${alignment%.*}`
		PREFIX=align
		NEWALIGNFILE=${PREFIX}${FILENAME}
		bowtie -S MtbCDC1551 $FILENAME.fq > $NEWALIGNFILE.sam
	done
	
	#clean up the FASTQ files
	rm *.fq

Count gene features/quantify transcript abundance with HTSeq
------------------------------------------------------------

So now that we have aligned the reads to the genome, we want to count transcript abundance
per gene using an annotation file (MtbCDC1551.gtf) that was generated when the genome is 
annotated. This file tells the software the gene coordinates and strandedness. 

Download the .gtf file::

    cd ../TranscriptAbund
    wget http://s3-us-west-1.amazonaws.com/microgenomicstranscriptomics/MtbCDC1551.gtf
    
Let's have a look at what the first few lines of this file look like::

    head MtbCDC1551.gtf
    
It should look something like this::

    Chromosome	protein_coding	exon	1	1524	.	+	.	 gene_id "MT0001"; transcript_id "AAK44224"; exon_number "1"; gene_name "dnaA"; transcript_name "dnaA/AAK44224"; seqedit "false";
    Chromosome	protein_coding	CDS	1	1521	.	+	0	 gene_id "MT0001"; transcript_id "AAK44224"; exon_number "1"; gene_name "dnaA"; transcript_name "dnaA/AAK44224"; protein_id "AAK44224";
	.
	.
	.

We can see the entire exon, it's coordinates, the strand, gene_id, gene_name, etc. The software
we will use (HTSeq) requires this file to "know" how to interpret the alignment file to count
whether a transcript was observed.

We need to strip off prefixes again and we encounter another new command called copy (cp).
The general syntax is similar to move (mv)::

    cp [the file you want to copy] [where you want to copy the file to]
    
So let's take off the prefixes and copy the stripped files to the TranscriptAbund folder::

    #strip off the align prefix and move the files into the TranscriptAbund folder
	for renamealign in ../Alignment/*.sam
	do
		LESSPREFIX=`basename ${renamealign//align/}`
		cp $renamealign ../TranscriptAbund/$LESSPREFIX
	done
	
HTSeq has a spectacular function called htseq-count. This function will quantify transcript
abundances for us. It also has several options we can specify, with two particularly important
ones that tell it how to call whether a read is within a coding region (-m) and whether our
data is stranded (--stranded). Have a look at the documentation on the three choices for the
-m option. Further, if you are doing a single-end RNA-seq experiment with the Illumina TruSeq
library preparation kit, your data will be reverse stranded. You can experiment with this to
see by specifying --stranded=forward or --stranded=no.

The syntax is as follows::

    htseq-count -m [how you want HTSeq to call a read in a gene region] --stranded=[reverse, forward, no] [alignment file that ends in .sam] [genome annotation file that ends with .gtf] > [your output file that ends in .counts]

Let's count!

To run the software::

	for counts in *.sam
	do
		FILENAME=`basename ${counts%.*}`
		PREFIX=abundcount_
		NEWMAPFILE=${PREFIX}${FILENAME}
		htseq-count -m intersection-nonempty --stranded=reverse $FILENAME.sam MtbCDC1551.gtf > $NEWMAPFILE.counts
	done
	
	#clean up the .sam files
	rm *.sam
	
Congratulations! We are now ready to do differential gene expression. 