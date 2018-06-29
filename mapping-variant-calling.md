# Mapping and variant calling on yeast transcriptome

TODO CTB:

* add instructions on downloading/viewing output HTML with RStudio
* put fastqc/multiqc output files in this repo

Learning objectives:

* define and explore the concepts and implications of shotgun
  sequencing;
  
* explore coverage;

* understand the basics of mapping-based variant calling;

* learn basics of actually calling variants & visualizing.

## Boot up a Jetstream

[Boot an m1.medium Jetstream instance](../lab1-jetstream/boot.md) and log in.

## Install software

```
conda install -y bwa samtools bcftools
```

## Download data

Goal: get the sequence data!

1. Run:

```
cd ~/
mkdir mapping
cd ~/mapping
```

Download some data from [Schurch et al, 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/):

```
curl -L https://osf.io/5daup/download -o ERR458493.fastq.gz
curl -L https://osf.io/8rvh5/download -o ERR458494.fastq.gz
curl -L https://osf.io/2wvn3/download -o ERR458495.fastq.gz
curl -L https://osf.io/xju4a/download -o ERR458500.fastq.gz
curl -L https://osf.io/nmqe6/download -o ERR458501.fastq.gz
curl -L https://osf.io/qfsze/download -o ERR458502.fastq.gz
```
        
...and take a quick look at it:

```
gunzip -c ERR458493.fastq.gz | head
```

## Map data

Goal: execute a basic mapping



3. Download and gunzip the reference:

```
curl -O https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
gunzip orf_coding.fasta.gz
```

and look at it:

```
head orf_coding.fasta
```
        
4. Prepare it for mapping:

```
bwa index orf_coding.fasta
```
        
5. Map!

```
bwa mem -t 4 orf_coding.fasta ERR458493.fastq.gz  > ERR458493.sam
```
        
6. Observe!

```
head ERR458493.sam
```

what does all this mean??
        
## Visualize mapping

Goal: make it possible to go look at a specific bit of the genome.

        
2. Index the reference genome:

```
samtools faidx orf_coding.fasta
```
        
3. Convert the SAM file into a BAM file:

```
samtools import orf_coding.fasta.fai ERR458493.sam ERR458493.bam
```
        
4. Sort the BAM file by position in genome:

```
samtools sort ERR458493.bam -o ERR458493.sorted.bam
```
        
5. Index the BAM file so that we can randomly access it quickly:

```
samtools index ERR458493.sorted.bam
```
        
6. Visualize with `tview`:

```
samtools tview ERR458493.sorted.bam orf_coding.fasta
```
        
   `tview` commands of relevance:
   
   * left and right arrows scroll
   * `q` to quit
   * CTRL-h and CTRL-l do "big" scrolls
   * `g NYLR162W:220` will take you to a specific location with two variants
   
7. Get some summary statistics as well:

```
samtools flagstat ERR458493.sorted.bam
```
   
## Call variants!

Goal: find places where the reads are systematically different from the
genome.
   
Now we can call variants using
[samtools mpileup](http://samtools.sourceforge.net/mpileup.shtml):

```
samtools mpileup -u -t DP -f orf_coding.fasta ERR458493.sorted.bam | \
    bcftools call -mv -Ov > variants.vcf
```

To look at the entire `variants.vcf` file you can do `cat
variants.vcf`; all of the lines starting with `#` are comments.  You
can use `tail variants.vcf` to see the last ~10 lines, which should
be all of the called variants.

## Discussion points / extra things to cover

* What are the drawbacks to mapping-based variant calling? What are
  the positives?

* Where do reference genomes come from?
