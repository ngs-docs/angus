# Genome assembly - some basics

Learning objectives:

	* Learn what is genome assembly?
	* Different types of assemblies
	* How do assemblers work?
	* Checking the quality of assembly
	* Understanding genome assembly
	
In [variant calling](http://angus.readthedocs.io/en/2018/mapping-variant-calling.html), we mapped reads to a reference and looked systematically for differences.

But what if you don't have a reference? How do you construct one?

The answer is *de novo* assembly, and the basic idea is you feed in your reads and you get out a bunch of *contigs*, that represent stretches of DNA present in the reads that don't have any long repeats or much significant polymorphism.  Like everything else, the basic idea is that you run a program, feed in the reads, and get out a pile of assembled DNA.

MEGAHIT, used below, works well for assembly short-read data sets from genomes and metagenomes.  For transcriptomes, you might use Trinity - see
[the eel-pond protocol](https://eel-pond.readthedocs.io/en/latest/) for our guide to doing RNA-seq assembly.

## Boot up a Jetstream

[Boot an m1.medium Jetstream instance](jetstream/boot.md) and log in.

## Install the MEGAHIT assembler

The [MEGAHIT assembler](https://www.ncbi.nlm.nih.gov/pubmed/27012178) is one of the first NGS metagenome assembler that can assemble genome sequences from metagenomic datasets of hundreds of Giga base-pairs (bp) in a time- and memory-efficient manner on a single server. We can also install it through `conda`:

```
conda install -y megahit quast 
```

## Change to a new working directory and link the original data

We will be using the same data as before ([Schurch et al, 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/)), so the following commands will create a new folder `assembly` and link the data in:

```
cd ~/
mkdir -p assembly
cd assembly

ln -fs ~/data/*.fastq.gz .
ls
```
    
### Run the assembler


Let's assemble the yeast data with MEGAHIT as follows:

```
megahit -r ERR458493.fastq.gz,ERR458494.fastq.gz,ERR458495.fastq.gz,ERR458500.fastq.gz,ERR458501.fastq.gz,ERR458502.fastq.gz -o yeast
```

(This will take about 3 minutes.)  You should see something like:

```
	--- [STAT] 8920 contigs, total 7438631 bp, min 200 bp, max 6752 bp, avg 834 bp, N50 1293 bp
	--- [Fri Jun 29 18:29:27 2018] ALL DONE. Time elapsed: 113.842486 seconds ---
```

at the end.

Questions while we're waiting:

* how many reads are there?

* how long are they?

* are they paired end or single-ended?

* are they trimmed?

...and how would we find out?

Also, what expectation do we have for this genome in terms of size, content, etc?


### Looking at the assembly

First, save the assembly:

```
cp yeast/final.contigs.fa yeast-assembly.fa
``` 
 
Now, look at the beginning:

```
head yeast-assembly.fa
```
    
It's DNA! Yay!

So this is the curse and the benefit of assembly - you go through some amount of work to get your data, QC it, clean it up, and assemble it, but then you're faced with a pile of assembled but unannotated results! (We'll talk about annotation in another lesson)

But before you put effort into annotating the assembly, you should think about whether it's any good...

### Measuring the assembly

Run QUAST on your assembly:

```
quast yeast-assembly.fa -o yeast_report
```

and now take a look at the report:

```
less yeast_report/report.txt
```

You should see something like:

```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                    yeast-assembly
# contigs (>= 0 bp)         8920
# contigs (>= 1000 bp)      2706
# contigs (>= 5000 bp)      6
# contigs (>= 10000 bp)     0
# contigs (>= 25000 bp)     0
# contigs (>= 50000 bp)     0
Total length (>= 0 bp)      7438631
Total length (>= 1000 bp)   4625556
Total length (>= 5000 bp)   35293
Total length (>= 10000 bp)  0
Total length (>= 25000 bp)  0
Total length (>= 50000 bp)  0
# contigs                   4861
Largest contig              6752
Total length                6184445
GC (%)                      39.33
N50                         1481
N75                         997
L50                         1463
L75                         2719
# N's per 100 kbp           0.00
```

This is a set of summary stats about your assembly. Are they good? Bad? How would you know?

## What are other metrics you could use to evaluate your assembly?

This is a good opportunity for brainstorming and group thinking :)

## End of day

Question: why so many contigs?!

And what do you do with a bunch of assembled DNA sequence anyway?
