Interval Analysis and Visualization
===================================

The results generate below are based on a question posed by a participant in the course.
She wanted to know how well contigs of an unfinished genomic build of and ecoli strain
match the common (K-12 strain MG1655) genomic build.

Download the results from:

http://apollo.huck.psu.edu/data/ms115.zip

How did we get the results in the file above? A short description follows:

Data collection
---------------

The partial genomic build is located at:

http://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161608

From this we downloaded the summary file ``code/ADTL01.txt``
that happens to be a tab delimited file that lists accession numbers.
We then wrote a very simple program ``code/getdata.py`` to parse
the list of accessions and download the data like so ::

    # requires BioPython
    from Bio import Entrez
    Entrez.email = "A.N.Other@example.com"
    stream = file("ADTL01.txt")
    stream.next()

    for line in stream:
        elems = line.strip().split()
        val = elems[1]
        handle = Entrez.efetch(db="nucleotide", id=val, rettype="fasta", retmode="text")
        fp = file("data/%s.fa" % val, 'wt')
        fp.write(handle.read())
        fp.close()

Finally we merged all data with::

	cat *.fa > MS115.fa

Then we went hunting for the EColi genome, we found it here:

http://www.genome.wisc.edu/sequencing/updating.htm

Turns out that this site only distributes a GBK (Genbank file).
We now need to extract the information from the
GBK file to FASTA (genome) and GFF (genomic feature) file. For this we need to
install the ReadSeq program:

http://iubio.bio.indiana.edu/soft/molbio/readseq/java/

Once we have this we typed::

	# GBK to GFF format
	java -jar readseq.jar -inform=2 -f 24 U00096.gbk

	# GBK to FASTA
	java -jar readseq.jar -inform=2 -f 24 U00096.gbk

This will create the files U00096.gbk.fasta and U00096.gbk.gff

Now lets map the ms115.fa contigs to the U00096.fa reference::

	bwa index U00096.fa
	bwa mem U00096.fa ms115.fa | samtools view -bS - | samtools sort - U00096

will produce the U00096.bam file. We have converted the U00096.bam to BED format
via the::

	bedtools bamtobed -i  U00096.bam  > U00096.bed

Visualizing the data

Download and run IGV

http://www.broadinstitute.org/igv/

Create a custom genome via `Genomes -> Create .genome` options

We will now  visualize the BAM, GFF and BED files and discuss the various aspects of it.

Running bedtools
================

Install bedtools::

	sudo apt-get bedtools

This works best if you store your files in Dropbox, that way you can
edit the file on your computer then load them up on your IGV instance.


