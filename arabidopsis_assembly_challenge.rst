===========================================================
Arabidopsis - De novo assembly vs Reference Guided Assembly
===========================================================



Goal:
====

Split into teams and choose one of these strategies:

1. De novo assembly with trinity, following instructions from earlier
2. Reference-guided assembly with trinity, following `these instructions from Trinity <https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly>`__.

   This will follow a slightly altered workflow:

   1. correct reads (same as before)
   2. run skewer (same as before)
   3. map with STAR (new step)

      Note: In addition to the aforementioned options, for GFF3 formatted annotations you need to use --sjdbGTFtagExonParentTranscript Parent.

   4. run trinity (altered command with bam input)
   5. run busco (same as before)
   6. run transrate (same as before)

I have notes for all of these commands if you need help.


Log your results here:
https://docs.google.com/spreadsheets/d/12X06LqGM8j4a4oV3_IsM91Hvop6B7L84xSWR71dXK3E/edit?usp=sharing


Get Data
========

Arabidopsis lyrata RNASeq Flower reads::

    wget http://www.hardwoodgenomics.org/sites/default/files/kbs_temp/SRR3993765_1.subsample.fastq
    wget http://www.hardwoodgenomics.org/sites/default/files/kbs_temp/SRR3993765_2.subsample.fastq

Arabidopsis thaliana Genome::

	wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

Arabidopsis thaliana Genome Annotation::

	wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff

We still have a slight problem - the chromosome names in the fasta file don't match the chromosome names in the annotation file. This is a surprisingly common problem and breaks any tool that needs both files. So lets fix the names.::

	sed -i 's/>\([1-5]\)/>Chr\1/' TAIR10_chr_all.fas
	sed -i 's/>mitochondria/>ChrM/' TAIR10_chr_all.fas
	sed -i 's/>chloroplast/>ChrC/' TAIR10_chr_all.fas

BUSCO plant:
===========

To get the plant database file::

	wget http://buscos.ezlab.org/files/plant_early_release.tar.gz
	gunzip plant_early_release.tar.gz


