# Finding Resources online

## Objectives

+ Discover common biological databases
+ Understand common bioinformatics file formats

## Biological Databases

There are many sources available for acquiring biological data. A partial, and by no means exhaustive, list includes:

+ NCBI
+ Genbank
+ RefSeq
+ Ensembl
+ Microbial Genome Database
+ DDBJ
+ FungiDB

Note that many of these databases cross-reference one-another, with the same depositions found in multiple sources.

## NCBI

+ [NCBI](https://www.ncbi.nlm.nih.gov/search/) - NCBI houses many major databases, including GenBank, which houses genome assemblies, and the SRA, which houses sequencing data. It's also home to databases like PubMed. 

## NCBI Taxonomy

NCBI Taxonomy is one particularly useful database. The Taxonomy Browser displays a hierarchical view of the taxonomy, as well as more detailed taxon-specific pages that highlight the internal links to other Entrez databases, the LinkOut links to external resources and other taxon-specific data. It is particularly useful for finding genomes and assemblies. 

1. Go to [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) page and search for "Saccharomyces cerevisiae"

![](/static/taxonomy1.png)

2. Click on the `1` in the `Genome` column to launch a new page. This page contains links to download the genome sequence and annotation files. 

![](/static/taxonomy2.png)

## NCBI SRA

NCBI Sequence Read Archive (SRA) stores sequence and quality data (fastq files) in aligned or unaligned formats from NextGen sequencing platforms.

1. Open [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) page and search for "Saccharomyces cerevisiae"
![](/static/sra1.png)

2. click on any of the listed results to show project information
![](/static/sra2.png)

3. Under "Project Data" click on links to "SRA Experiments"
![](/static/sra3.png)

4. In, this page, you can find:
    + Experimental Design
    + Sample Information
    + Library Preparation Information
    + Sequencing "Run" information with number usually starting with "SRR"
![](/static/sra4.png)

### SRA ID's from publications

Most journals require researches to publish their sequencing data in a public repository upon manuscript publication.
There's often a section at the end of a publication that indicates the sequence accession number. For example, the 
dataset we're using in these lessons comes from this [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/). The text says the following:

```
DATA DEPOSITION

The data sets supporting the results of this article are available in the European Nucleotide Archive repository (ENA) (PRJEB5348, http://www.ebi.ac.uk/ena/data/view/ERX425102). All the code for this work is publicly available (https://github.com/bartongroup/profDGE48).
```

The ENA and SRA mirror each other, so we can view this accession from either location. It is often easier to download 
single or a handful of accessions from the ENA. In our R lesson, we will work with the ENA table from this dataset. 

## File formats

There are a lot of different file formats in bioinformatics! We've listed a few common ones below:

+ FASTA – plain sequences (.fa, .fasta, .faa, .fna, .fnn)
+ FASTQ – sequencing reads
+ GFF – general feature format
+ GTF - variation of GFF
+ VCF – sequence variants
+ SAM – sequence alignments
+ BAM – alignments in binary (compressed)


## Optional: Ensemble

Ensembl is another useful database if you're working with a model organism. It provides a browser for vertebrate (and sometimes other) genomes and supports research in comparative genomics, evolution, sequence variation and transcriptional regulation. Ensembl annotates genes, computes multiple alignments, predicts regulatory function and collects disease data. Ensembl tools include BLAST, BLAT, BioMart and the Variant Effect Predictor (VEP) for all supported species. Below we demonstrate how to find annotated genes using Ensemble.


1. Open [Ensemble](https://www.ensembl.org/index.html?redirect=no) webpage and search for "Saccharomyces cerevisiae"
![](/static/ensemble1.png)

2. Under the "Gene Annotation" section you can download FASTA files for genes, cDNA, ncRNA, proteins and annotations for the same
![](/static/ensemble2.png)

3. This is the ftp site of ensemble 
![](/static/ensemble3.png)
![](/static/ensemble4.png)

