# De novo transcriptome assembly

Learning objectives:

* What is transcriptome assembly?
* How do assemblers work?
* Checking the quality of assembly
* Understanding transcriptome assembly

What if you are working with a species with no existing reference genome or transcriptome. How do you construct one?

The answer is *de novo* assembly. The basic idea with *de novo* transcriptome assembly is you feed in your reads and you get out a bunch of *contigs* that represent transcripts, or stretches of RNA present in the reads that don't have any long repeats or much significant polymorphism. You run a  *de novo* transcriptome assembly program using the trimmed reads as input and get out a pile of assembled RNA.

Trinity, one of the leading *de novo* transcriptome assemblers, was developed at the [Broad Institute](http://www.broadinstitute.org/) and the [Hebrew University of Jerusalem](http://www.cs.huji.ac.il/). We will be losely following steps from the [Eel pond protocol](https://eel-pond.readthedocs.io/en/latest) for our guide to doing RNA-seq assembly.

## Download and trim the data

We will be using a set of *Nematostella vectensis* mRNAseq reads from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16).

To avoid needing to go though the trimming steps as before, we'll download a snakefile to take care of these steps for us. If you look at this file, you may notice it is very similar to the one you generated during the snakemake challenge.


Download the snakefile and a corresponding conda environment file:
```
cd ~ # cd to home directory
curl -L https://osf.io/nqh6p/download -o nema_trim.snakefile
curl -L https://osf.io/xs6k7/download -o trim-env.yml
```

Run the snakefile to download and trim the Nematostella reads:
```
snakemake -s nema_trim.snakefile --use-conda --cores 6
```

The trimmed data should now be within the `nema_trimmed` folder.


## Generate a _de novo_ assembly

The following commands will create a new folder `assembly` and link the trimmed data we prepared earlier in the newly created folder:

```
cd 
mkdir assembly
cd assembly

ln -fs ../nema_trimmed/*.qc.fq.gz .
ls
```

Combine all fq into 2 files (left.fq and right.fq)
```
cat *_1.pe.qc.fq.gz *se.qc.fq.gz > left.fq.gz
cat *_2.pe.qc.fq.gz > right.fq.gz
```

_Note: this step just makes it easier for us to type out the trinity command. Trinity can accept comma-separated lists of files as inputs - check the trinity documentation for details._


## Run the assembler

Trinity works both with paired-end reads as well as single-end reads (including simultaneously both types at the same time). In the general case, the paired-end files are defined as `--left left.fq` and `--right right.fq` respectively. Our single-end reads (orphans) have been concatenated onto the `left.fq` file. 


So let's run the assembler as follows:

```
time Trinity --seqType fq --max_memory 16G --CPU 6 --left left.fq.gz --right right.fq.gz --output nema_trinity
```

(This will take about 5 minutes)

You should see something like:

```
** Harvesting all assembled transcripts into a single multi-fasta file...

Thursday, October 25, 2018: 21:55:15	CMD: find /home/dibbears/work/assembly/nema_trinity/read_partitions/ -name '*inity.fasta'  | /opt/miniconda3/opt/trinity-2.8.4/util/support_scripts/partitioned_trinity_aggregator.pl --token_prefix TRINITY_DN --output_prefix /home/dibbears/work/assembly/nema_trinity/Trinity.tmp
-relocating Trinity.tmp.fasta to /home/dibbears/work/assembly/nema_trinity/Trinity.fasta
Thursday, October 25, 2018: 21:55:15	CMD: mv Trinity.tmp.fasta /home/dibbears/work/assembly/nema_trinity/Trinity.fasta


###################################################################
Trinity assemblies are written to /home/dibbears/work/assembly/nema_trinity/Trinity.fasta
###################################################################


Thursday, October 25, 2018: 21:55:15	CMD: /opt/miniconda3/opt/trinity-2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/dibbears/work/assembly/nema_trinity/Trinity.fasta > /home/dibbears/work/assembly/nema_trinity/Trinity.fasta.gene_trans_map

real	7m7.692s
user	23m59.929s
sys	13m32.485s
```

at the end.



## Looking at the assembly

First, save the assembly:

```
cp nema_trinity/Trinity.fasta nema-transcriptome-assembly.fa
``` 
 
Now, look at the beginning:

```
head nema-transcriptome-assembly.fa
```
    
These are the transcripts! Yay!

Let's capture also some statistics of the Trinity assembly. Trinity provides a handy tool to do exactly that:

```
TrinityStats.pl nema-transcriptome-assembly.fa
```

The output should look something like the following:

```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	217
Total trinity transcripts:	220
Percent GC: 48.24

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 1763
	Contig N20: 819
	Contig N30: 548
	Contig N40: 407
	Contig N50: 320

	Median contig length: 245.5
	Average contig: 351.60
	Total assembled bases: 77353


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 1034
	Contig N20: 605
	Contig N30: 454
	Contig N40: 357
	Contig N50: 303

	Median contig length: 245
	Average contig: 328.43
	Total assembled bases: 71270
```

This is a set of summary stats about your assembly. Are they good? Bad? How would you know?

## Suggestions for next steps 

After generating a *de novo* transcriptome assembly:
* [annotation](https://angus.readthedocs.io/en/2018/dammit_annotation.html)
* [evaluation](https://dibsi-rnaseq.readthedocs.io/en/latest/evaluation.html)
