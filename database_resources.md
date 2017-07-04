# Publicly available databases
There are many, many databases around for sequence data and for downstream analysis of sequence data. Below we have listed some of hte most commone ones and their function. This is not an exhaustive list, but it can help you get started with finding relevant data to help with your analysis.

## Finding data of interest
* [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed): citations, abstract and links to more than 27 million scientific papers.
* [Google scholar](https://scholar.google.se)
* [Dryad](https://datadryad.org/): a "curated general-purpose repository that makes the data underlying scientific publications discoverable, freely reusable, and citable." Integrated with many journals.
* [FigShare](https://figshare.com/): is a repository where users can make all of their research outputs available in a citable, shareable and discoverable manner.

## NCBI
[NCBI](https://www.ncbi.nlm.nih.gov) has a lot of really wonderful resources. These all have different interfaces, and some are better organized than others, but the data housed within the various databases is gold. 
+ [GEO](https://www.ncbi.nlm.nih.gov/geo/) (Gene expression omnibus): gene expression data; array- and sequence-based data are catalogued within. Experimental design is also reported, although some experiments give more details than others. 
+ [Assembly](https://www.ncbi.nlm.nih.gov/assembly/): organisms with genomic assemblies.
+ [Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy): names of organisms, taxonomic ID. It can be a bit of a mess, but also really useful when you go to do anything with phylogeny.
+ [Sequence read archive](https://www.ncbi.nlm.nih.gov/sra): raw data files. Can filter based on DNA, RNA, whole genome sequencing, organism, etc. See below for accessing & downloading this data.
+  [WGS](https://www.ncbi.nlm.nih.gov/genbank/wgs/): Whole Genome Shotgun projects (complete or incomplete assemblies)
+  Many others!


## Downloading data from NCBI
+ [European nucleotide archive](http://www.ebi.ac.uk/ena): links to fastq files. You can search for SRA project data here to download fastq files & avoid SRA format (below).
+ [SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software): command-line interface -- recommended only if you have many samples to download

## Other Protein data
+ [UniProt](http://www.uniprot.org)
    + Uniprot is composed of 2 resources: Swissprot and TrEMBL. Swissprot is a databse of manually curated protein sequences (very high quality!) while trEMBL is automatically annotated (but contains a lot more sequences)
+ [NCBI protein](https://www.ncbi.nlm.nih.gov/protein/)
    + A database that includes protein sequence records from a variety of sources, including GenPept, RefSeq, Swiss-Prot, PIR, PRF, and PDB.

## Genomes & Genome Browsers
+ [Ensembl](http://www.ensembl.org/index.html): a genome browser that supports research in comparative genomics, evolution, sequence variation and transcriptional regulation. Ensembl annotate genes, computes multiple alignments, predicts regulatory function and collects disease data. Ensembl tools include BLAST, BLAT, BioMart and the Variant Effect Predictor (VEP) for all supported species. Also contains reference genomes and annotation files that can be downloaded.
    + [EnsemblMetazoa](http://metazoa.ensembl.org/index.html): The same as Ensembl but addition organisms that are not available on the primary Ensembl page.
+ [UCSC Genome browser](https://genome.ucsc.edu): contains many reference genomes, and many tools to search these including BLAT, in silico PCR, LiftOver (lift sequence from organism onto another) and many other cool and useful tools.
+ [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/): a non-redundant and well-annotated set of reference sequences including genomic, transcript, and protein.
+ [GENCODE](http://www.gencodegenes.org): high quality reference gene annotation and experimental validation for human and mouse genomes.
+ Joint Genome Institute
    + [Mycocosm](http://genome.jgi.doe.gov/programs/fungi/index.jsf): 1000 fungal genomes project. Well-annotated and diverse fungal reference genomes with many supporting tools.
    + [GOLD](https://gold.jgi.doe.gov/): genomes online database
    + [Integrated Microbial Genomes & Microbiomes](https://img.jgi.doe.gov/)

## Other databases full of many things
+ [GenBank](https://www.ncbi.nlm.nih.gov/genbank/): NIH genetic sequence database, an annotated collection of all publicly available DNA sequences
+ [EBML-EBI](http://www.ebi.ac.uk/): European bioinformatics institute. 
+ [DDBJ](http://www.ddbj.nig.ac.jp/): DNA databank of Japan



## BLAST
+ [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
+ [UniProt](http://www.uniprot.org/blast/)
+ Ensembl: organism-specific blast supported. 

## Pathways
+ [KEGG](http://www.genome.jp/kegg/pathway.html)
+ [MetaCyc](https://metacyc.org/)
+ [BioCyc](https://biocyc.org/)
+ [Pathway tools](http://brg.ai.sri.com/ptools/)
+ [GO](http://www.geneontology.org)
+ [MsigDB](http://software.broadinstitute.org/gsea/msigdb)

## Metagenomes
+ [MG-RAST](http://metagenomics.anl.gov/): full of great data, some reports of odd uploading and downloading so be careful when using this resource!
+ [EBI metagenomics](https://www.ebi.ac.uk/metagenomics/)
+ [UniMES](http://www.uniprot.org/help/unimes): UniProt Metagenomic and Environmental Sequences

## Marine organism resources 
- [Aniseed](www.aniseed.cnrs.fr): Ascidian Network for in situ Expression and Embryological Data
- [Echinobase](http://www.echinobase.org/Echinobase/): Several echinoderm genomes, expression data also available 
- [OIST Marine Genomic Unit](http://marinegenomics.oist.jp/cots/gallery): Several genomes and omics data for various marine organisms

## Human Genomes
+ [1000 genomes project](http://www.internationalgenome.org/data/)
+ [The Cancer Genome Atlas](https://cancergenome.nih.gov/)
+ [The Human Microbiome Project](http://hmpdacc.org/)

## Tool Aggregators
* [OMICtools](https://omictools.com/)

## Blogs and other useful links
* [RNA-Seq Blog](http://www.rna-seqblog.com/)
* [Titus's blog](http://ivory.idyll.org/blog/)
* [Staying Current in Bioinformatics](http://www.gettinggeneticsdone.com/2017/02/staying-current-in-bioinformatics-genomics-2017.html) - Stephen Turner
