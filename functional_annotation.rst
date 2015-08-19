===========================================================================
Functional Annotation of transcripts
===========================================================================

In this tutorial, we'll use some sample data from green ash to demonstrate functional annotation of transcript sequences. Green ash is from a relatively under-sampled lineage of plants, the order Lamiales. We have an assembled transcriptome, but no information on the function of these transcripts.

Amazon log in

1. Launch an EC2 instance and log in (Ubuntu 14.04, c3.4xlarge). Note the availability zone.
2. Click "Snapshots" item in the menu.  Click on dropdown that says "Owned by me" and select "Public Snapshots". In the search box, enter snap-153de366, and search. A snapshot with the name "KBS NGS 2015 functional Annotation" should appear. Select it.
3. Once you have the correct snapshot selected, choose "Create Volume" from the dropdown Actions menu.  Select the same availability zone and create the volume.
4. Click on the "volumes" link in the menu. Your volume should be available. If you have more than one, click on each and look for the one that has the correct snapshot.
5. Under actions, click "Attach Volume". Select the running instance. Click Attach.

Now log in as usual. We need to mount the volume:

::

	sudo bash 
	mkdir /mnt/ebs 
	mount /dev/xvdf /mnt/ebs 
	cd /mnt/ebs/

Next, we'll need to install a bunch of software. Some of these tools can be installed using apt-get. Note that apt-get does not necessarily always install the most up-to-date versions of this software! You should always double check versions when you do this. 

::

    #Housekeeping stuff for a new instance
    sudo bash
    apt-get update
    apt-get -y install screen git curl gcc make g++ python-dev unzip default-jre pkg-config libncurses5-dev r-base-core r-cran-gplots python-matplotlib sysstat python ncbi-blast+ sqlite 
    apt-get hmmer

I have preloaded a set of 1000 green ash transcripts to annotate. They are located at

::

	/mnt/ebs/green_ash_transcripts.fasta

I have also preloaded software and data in 

::

	/mnt/ebs/functional_annotation

Lets make a directory to hold all of our analysis.

::

	cd /mnt/ebs
	mkdir green_ash
	cd green_ash
	mv ../green_ash_transcripts.fasta .

First, lets find open reading frames (orfs). we will do this with `transdecoder <https://transdecoder.github.io/>`_. 

TransDecoder identifies likely coding sequences based on the following criteria:

- a minimum length open reading frame (ORF) is found in a transcript sequence
- a log-likelihood score similar to what is computed by the GeneID software is > 0.
- the above coding score is greatest when the ORF is scored in the 1st reading frame as compared to scores in the other 5 reading frames.
- if a candidate ORF is found fully encapsulated by the coordinates of another candidate ORF, the longer one is reported. However, a single transcript can report multiple ORFs (allowing for operons, chimeras, etc).
- optional the putative peptide has a match to a Pfam domain above the noise cutoff score. (We will see this later!)


Install software:

::

	cd /mnt/ebs/functional_annotation
	curl -L https://github.com/transdecoder/transdecoder/archive/2.0.1.tar.gz > transdecoder.tar.gz
	tar xzf transdecoder.tar.gz
	cd TransDecoder-2.0.1
	make

Add paths:

::
	cd /mnt/ebs/functional_annotation/TransDecoder-2.0.1
	PATH=$PATH:$(pwd)
	cd /mnt/ebs/functional_annotation/interproscan-5.14-53.0
	PATH=$PATH:$(pwd)

and run. This is the simplest use of transdecoder, not using homology searches.

::

	cd /mnt/ebs/green_ash
	TransDecoder.LongOrfs -t green_ash_transcripts.fasta -m 30

your initial predictions are now available in the file:

::

	/mnt/ebs/green_ash/green_ash_transcripts.fasta.transdecoder_dir/longest_orfs.pep

lets check the number of sequences:

:: 

	grep -c '^>' longest_orfs.pep

Lots of ORFs!

PFam
~~~~

Next, lets use `Pfam <http://pfam.xfam.org/>`_. Pfam is a collection of protein families, represented by multiple sequence alignments and hidden markov models (HMMs). we can use the tool `HMMER <http://hmmer.janelia.org/>`_  to compare our sequences to the PFam database. 

(HMMER is to Pfam as blast is to refseq).

The PFam database is available from the PFam website. It is a bit large (1.3Gb), so its already been placed on your snapshot. the location of the PFam database is:

::

	/mnt/ebs/functional_annotation/Pfam-A.hmm

However, HMMER would like this file to be compressed further, for speed. lets do that next

::

	cd /mnt/ebs/functional_annotation
	hmmpress Pfam-A.hmm

While we wait, lets learn about E-values.

------------

We need a way to assess the strength of evidence for homology. I.E. How likely would it be to see this alignment or match by chance?

Expect value (E) 
	a parameter that describes the number of hits one can "expect" to see 
	by chance when searching a database of a particular size. It decreases 
	exponentially as the Score (S) of the match increases. Essentially, the 
	E value describes the random background noise. For example, an E value 
	of 1 assigned to a hit can be interpreted as meaning that in a database 
	of the current size one might expect to see 1 match with a similar score 
	simply by chance.

Find lots more info: `Statistics of Sequence Similarity Scores <http://www.ncbi.nlm.nih.gov/BLAST/tutorial/>`_

------------

Lets run the hmmer tool to compare our green ash transcripts to the pfam database.

:: 

	cd /mnt/ebs/green_ash
	## make a copy of the peptide sequence file for convenience
	cp green_ash_transcripts.fasta.transdecoder_dir/longest_orfs.pep ./green_ash_peptides.fasta
	## and run
	hmmscan --cpu 14 --domtblout green_ash_peptides.pfam.out /mnt/ebs/functional_annotation/Pfam-A.hmm green_ash_peptides.fasta > pfam.log
	
Check out the output file, green_ash_peptides.pfam.out

:: 

	less green_ash_peptides.pfam.out

------------

less command

	space to scroll forward

	q to quit

------------

Output formats (easier to look at in Excel):

- [1] target name: The name of the target sequence or profile.
- [4] query name: The name of the query sequence or profile.
- [7] E-value: E-value of the overall sequence/profile comparison (including all domains).
- [11] of: The total number of domains reported in the sequence, ndom.
- [23] Description of target

The rest of the column descriptions can be found in the `Hmmer User Guide <ftp://selab.janelia.org/pub/software/hmmer/CURRENT/Userguide.pdf>`_.


Transdecoder (again and better)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run again, this time using the Pfam results to guide ORF identification. 

::

	cd /mnt/ebs/green_ash
	TransDecoder.Predict \
	-t green_ash_transcripts.fasta \
	--retain_pfam_hits green_ash_peptides.pfam.out

New predictions are now available in the file:

::

	green_ash_transcripts.fasta.transdecoder.pep

lets check the number of sequences:

:: 

	grep -c '^>' green_ash_transcripts.fasta.transdecoder.pep

We have much more reasonable ORF predictions now.


InterProScan
~~~~~~~~~~~~

`InterProScan <https://www.ebi.ac.uk/interpro/interproscan.html>`_ is a bit complex - this tool is a wrapper for 11 different databases. They all provide protein domain/function information.

The download for interproscan includes most of the databases, so it is quite large (3.3gb). It is on the snapshot.

::

	/mnt/ebs/functional_annotation/interproscan-5.14-53.0

Lets run without parameters to see what is available.

::

	interproscan.sh | more

We need to remove the astericks from the peptide file - ips does not like these. This is a common problem - astericks are often used to denote a stop codon.

::

	sed -i 's_*__g' /mnt/ebs/green_ash/green_ash_peptides.fasta

And we can make the software faster. It does not accept a parameter on the command line to increase the number of processors used, but it does have a properties file. Lets edit it.

::
	nano /mnt/ebs/functional_annotation/interproscan-5.14-53.0/interproscan.properties

Change 

::

	number.of.embedded.workers=1
	maxnumber.of.embedded.workers=2

To

::

	number.of.embedded.workers=14
	maxnumber.of.embedded.workers=15

Save with (Control-O, enter to save, Control-X to exit).
	
And we will now make a results directory and run the software

::

	mkdir /mnt/ebs/green_ash/ips_results
	interproscan-5.14-53.0/interproscan.sh \
	-d /mnt/ebs/green_ash/ips_results \
	-dp \
	-goterms \
	-i /mnt/ebs/green_ash/green_ash_peptides.fasta \
	-iprlookup \
	-pa

----------

Parameters

-dp         disable precalculation
-goterms    lookup the GO terms associated with the database entry
-iprlookup  lookup the global InterPro accession number
-pa         lookup pathway annotation

----------

Check out results

::

	cd ips_results/
	less green_ash_peptides.fasta.tsv
	wc -l green_ash_peptides.fasta.tsv


blast to swiss-prot
~~~~~~~~~~~~~~~~~~~

`Swiss-prot <http://web.expasy.org/docs/swiss-prot_guideline.html>`_ is the manually annotated and reviewed section of the uniprot knowledgebase (uniprotkb).

I've already downloaded it on our snapshot to:

::

	/mnt/ebs/functional_annotation/uniprot_sprot.fasta

Lets check out the number of reads with our handy one liner for fasta files

::

	grep -c '^>' uniprot_sprot.fasta

And we need to format this for blast searching

::

	makeblastdb -in uniprot_sprot.fasta -dbtype prot

Now lets get back into our functional_annotation directory and run blast. we can run blastp to search our predicted transcript proteins (orfs) to known proteins, or we can run blastx to search all transcripts against known proteins. lets do the latter. This is an important step, because we don't know if transdecoder found the correct open reading frame.

::

	blastx \
	-query /mnt/ebs/green_ash/green_ash_transcripts.fasta \
	-db /mnt/ebs/functional_annotation/uniprot_sprot.fasta \
	-out green_ash_transcripts-vs-swissprot.blastx \
	-num_threads 15

----------

How do you know which blast program you need?

======= ===================================	===================================
Program Database	                        Query 
======= ===================================	===================================
blastp  protein                             protein
blastn	nucleotide                          nucleotide
blastx  protein                             nucleotide (translated in 6 frames)
tblastn nucleotide (translated in 6 frames) protein
tblastx nucleotide (translated in 6 frames) nucleotide (translated in 6 frames)
======= ===================================	===================================

----------

This returns a lot of results in a relatively unusable format. lets try again with some new parameters

::

	blastx \
	-query /mnt/ebs/green_ash/green_ash_transcripts.fasta \
	-db /mnt/ebs/functional_annotation/uniprot_sprot.fasta \
	-out green_ash_transcripts-vs-swissprot.blastx.tsv \
	-num_threads 15 \
	-outfmt "6 std stitle" \
	-max_target_seqs 1


Lets figure out how many transcipts have a match with a bash one liner.

::

	cut -f1 green_ash_transcripts-vs-swissprot.blastx.tsv  | sort | uniq | wc -l

This is good information to have, especially if the open reading frame was incorrectly identified.


We can also search the peptide sequences.

::

	blastp \
	-query /mnt/ebs/green_ash/green_ash_peptides.fasta \
	-db /mnt/ebs/functional_annotation/uniprot_sprot.fasta \
	-out green_ash_peptides-vs-swissprot.blastx.tsv \
	-num_threads 15 \
	-outfmt "6 std stitle" \
	-max_target_seqs 1

---------


There are a lot of ways to do functional annotation. Here are a few more:

- `KEGG Automatic Annotation Server <http://www.genome.jp/tools/kaas/>`_
- `Trinotate <https://trinotate.github.io/>`_ (This includes licensed tools like `RNAMMER <http://www.cbs.dtu.dk/cgi-bin/sw_request?rnammer>`_)


