# Mapping and variant calling on yeast transcriptome

TODO:

* find a better variant with variants.vcf

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
mkdir -p mapping
cd mapping
```
Next, we will create links from the previous download yeast dataset:
```
ln -fs ~/data/*.fastq.gz .
ls
```

## Map data

Goal: execute a basic mapping

### Download and gunzip the reference:

```
curl -O https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
gunzip orf_coding.fasta.gz
```

and look at it:

```
head orf_coding.fasta
```
        
### Prepare it for mapping:

```
bwa index orf_coding.fasta
```
        
### Map!

```
bwa mem -t 4 orf_coding.fasta ERR458493.fastq.gz  > ERR458493.sam
```
        
### Observe!

```
head ERR458493.sam
```

what does all this mean??
        
## Visualize mapping

Goal: make it possible to go look at a specific bit of the genome.

### Index the reference genome:

```
samtools faidx orf_coding.fasta
```
        
### Convert the SAM file into a BAM file:

```
samtools import orf_coding.fasta.fai ERR458493.sam ERR458493.bam
```
        
### Sort the BAM file by position in genome:

```
samtools sort ERR458493.bam -o ERR458493.sorted.bam
```
        
### Index the BAM file so that we can randomly access it quickly:

```
samtools index ERR458493.sorted.bam
```
        
### Visualize with `tview`:

```
samtools tview ERR458493.sorted.bam orf_coding.fasta
```
        
   `tview` commands of relevance:
   
   * left and right arrows scroll
   * `q` to quit
   * CTRL-h and CTRL-l do "big" scrolls
   * Typing `g` allows you to go to a specific location, in this format chromosome:location. Here are some locations you can try out:
   -- `YLR162W:293` (impressive pileup, shows two clear variants and three other less clear)
   -- `YDR034C-A:98` (impressive pileup, shows two clear variants)
   -- `YDR366C:310` (impressive pileup, less clear variants)
   -- `YLR256W:4420` (impressive pileup, less clear variants)
   -- `YBL105C:2179` (less depth, shows two clear variants)
   -- `YDR471W:152` (impressive pileup, shows one clear variant)
   
Get some summary statistics as well:

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
