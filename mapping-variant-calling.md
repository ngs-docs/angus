# Mapping and variant calling on yeast transcriptome

Learning objectives:

* define and explore the concepts and implications of shotgun
  sequencing;
  
* explore coverage;

* understand the basics of mapping-based variant calling;

* learn basics of actually calling variants & visualizing.

## Boot up a Jetstream

[Boot an m1.medium Jetstream instance](jetstream/boot.md) and log in.

## Install software

```
conda install -y bwa samtools bcftools
```

## Change to a new working directory and map data
After installing the necessary software, we will create the working directory for the mapping as follows:
```
cd ~/
mkdir mapping
cd mapping
```
Next, we will create links from the previously downloaded and quality-trimmed yeast dataset:
```
ln -fs ~/quality/*.qc.fq.gz .
ls
```

## Variant Calling Workflow

<center><img src="_static/variant_calling_workflow.png" width="50%"></center>
<br>

## Map data

Goal: perform read alignment or mapping to determine where in the genome our reads originated from.

### Download and gunzip the reference:

Here we are using **open coding regions** to do variant calling because we are working with mRNA sequences.
It's important to think about what reference is appropriate for your experiment. Many biologically important
variants exist in non-coding regions, so if we were looking at genomic sequences, it would be important to
use a different reference such as the whole genome.

```
curl -O https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
gunzip orf_coding.fasta.gz
```

Let's take a look at our reference:

```
head orf_coding.fasta
```
        
### Indexing: Prepare reference for mapping:

Our first step is to index the reference genome for use by BWA. Indexing allows the aligner to quickly find potential
alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to 
be run once. The only reason you would want to create a new index is if you are working with a different reference genome
or you are using a different tool for alignment.

```
bwa index orf_coding.fasta
```
        
### Map!

We use an algorithm called `bwa mem` to perform mapping.

```
bwa mem -t 4 orf_coding.fasta ERR458493.qc.fq.gz  > ERR458493.sam
```
Have a look at the [bwa options](http://bio-bwa.sourceforge.net/bwa.shtml) page. While we are running bwa with the default parameters here, your use case might require a change of parameters. NOTE: Always read the manual page for any tool before using and make sure the options you use are appropriate for your data.
        
**What is the difference between Salmon and bwa mem?** In an earlier tutorial, we used Salmon "quasi-mapping" to generate
counts for transcripts. Salmon uses exact matching of k-mers (sub-strings in reads) to approximate which read a transcipt
originated from. This is enough information for read quantification, and is really fast. BWA `mem` produces an alignment,
where an entire read is mapped exactly against a reference sequence. This produces more information that is important
for things like variant calling. 

### Observe!

```
head ERR458493.sam
tail ERR458493.sam
```

The SAM file is a tab-delimited text file that contains information for each individual read and its 
alignment to the genome. While we do not have time to go in detail of the features of the SAM format, 
the paper by [Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

The compressed binary version of SAM is called a BAM file. We use this version to reduce size and to 
allow for indexing, which enables efficient random access of the data contained within the file.

The file begins with a header, which is optional. The header is used to describe source of data, 
reference sequence, method of alignment, etc., this will change depending on the aligner being used. 
Following the header is the alignment section. Each line that follows corresponds to alignment information 
for a single read. Each alignment line has 11 mandatory fields for essential mapping information and a 
variable number of other fields for aligner specific information. An example entry from a SAM file is 
displayed below with the different fields highlighted.

<center><img src="_static/sam_file.png" width="60%"></center>
<br>

<center><img src="_static/bam_file.png" width="60%"></center>
<br>

**Challenge**

Construct a `for` loop to to generate `.sam` alignment files for all of our quality controlled samples.

## Visualize mapping

Goal: make it possible to go look at a specific bit of the genome.

### Index the reference genome:

Before we indexed the reference for BWA, now we reference the index for samtools. Although both
tools use different indexing methods, they both allow the tools to find specific sequences within
the reference quickly.

```
samtools faidx orf_coding.fasta
```
        
### Convert the SAM file into a BAM file:

```
samtools import orf_coding.fasta.fai ERR458493.sam ERR458493.bam
```
        
### Sort the BAM file by position in genome:

You can sort on many different fields within a sam or bam file. After mapping, our
files are sorted by read number. Most programs require mapping files to be sorted by
position in the reference. You can sort a file using the `samtools sort` command.

```
samtools sort ERR458493.bam -o ERR458493.sorted.bam
```
        
### Index the BAM file so that we can randomly access it quickly:

```
samtools index ERR458493.sorted.bam
```
        
### Visualize with `tview`:

Samtools implements a very simple text alignment viewer based on the GNU ncurses library, called tview. 
This alignment viewer works with short indels and shows MAQ consensus. It uses different colors to display 
mapping quality or base quality, subjected to users' choice. Samtools viewer is known to work with an 
130 GB alignment swiftly. Due to its text interface, displaying alignments over network is also very fast.

In order to visualize our mapped reads we use tview, giving it the sorted bam file and the reference file:

```
samtools tview ERR458493.sorted.bam orf_coding.fasta
```
        
   `tview` commands of relevance:
   
   * left and right arrows scroll
   * `q` to quit
   * CTRL-h and CTRL-l do "big" scrolls
   * Typing `g` allows you to go to a specific location, in this format chromosome:location. Here are some locations you can try out:
   
     - `YLR162W:293` (impressive pileup, shows two clear variants and three other less clear)
     - `YDR034C-A:98` (impressive pileup, shows two clear variants)
     - `YDR366C:310` (impressive pileup, less clear variants)
     - `YLR256W:4420` (impressive pileup, less clear variants)
     - `YBL105C:2179` (less depth, shows two clear variants)
     - `YDR471W:152` (impressive pileup, shows one clear variant)
   
Get some summary statistics as well:

```
samtools flagstat ERR458493.sorted.bam
```
   
## Call variants!

Goal: find places where the reads are systematically different from the
genome.

A variant call is a conclusion that there is a nucleotide difference vs. some reference at a given 
position in an individual genome or transcriptome, often referred to as a Single Nucleotide Polymorphism (SNP). 
The call is usually accompanied by an estimate of variant frequency and some measure of confidence. Similar 
to other steps in this workflow, there are number of tools available for variant calling. In this workshop 
we will be using `bcftools`, but there are a few things we need to do before actually calling the variants.


```
bcftools mpileup -O b -f orf_coding.fasta ERR458493.sorted.bam | \
    bcftools call -m -v -o variants.vcf
```

```
vcfutils.pl varFilter variants.vcf  > variants_filtered.vcf
```

**Challenge** 

How many fewer lines are there in the `variants.vcf` vs. in the `variants_filtered.vcf` file?

To look at the entire `variants.vcf` file you can do `cat
variants.vcf`; all of the lines starting with `#` are comments.  You
can use `tail variants.vcf` to see the last ~10 lines, which should
be all of the called variants.

The first few columns of the VCF file represent the information we have about a predicted variation:

| column | info |
| ------ | ------ | 
| CHROM |	contig location where the variation occurs |
| POS |	position within the contig where the variation occurs |
| ID |	a `.` until we add annotation information |
| REF |	reference genotype (forward strand) |
| ALT |	sample genotype (forward strand) |
| QUAL |	Phred-scaled probablity that the observed variant exists at this site (higher is better) |
| FILTER |	a `.` if no quality filters have been applied, PASS if a filter is passed, or the name of the filters this variant failed |

The Broad Institute’s [VCF guide](https://www.broadinstitute.org/gatk/guide/article?id=1268) is an excellent place to learn more about VCF file format.
