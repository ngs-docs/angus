# Transcriptome assembly - some basics

Learning objectives:

	* Learn what is Transcriptome assembly?
	* Different types of assemblies
	* How do assemblers work?
	* Checking the quality of assembly
	* Understanding Transcriptome assembly
	
In [variant calling](http://angus.readthedocs.io/en/2018/mapping-variant-calling.html), we mapped reads to a reference and looked systematically for differences.

But what if you don't have a reference? How do you construct one?

The answer is *de novo* assembly, and the basic idea is you feed in your reads and you get out a bunch of *contigs*, that represent stretches of DNA present in the reads that don't have any long repeats or much significant polymorphism.  Like everything else, the basic idea is that you run a program, feed in the reads, and get out a pile of assembled DNA.

Trinity, developed at the [Broad Institute](http://www.broadinstitute.org/) and the [Hebrew University of Jerusalem](http://www.cs.huji.ac.il/, represents a novel method for the efficient and robust de novo reconstruction of transcriptomes from RNA-seq data. We will be using the [eel-pond protocol](https://eel-pond.readthedocs.io/en/latest) for our guide to doing RNA-seq assembly.

## Boot up a Jetstream

[Boot an m1.medium Jetstream instance](jetstream/boot.md) and log in.

## Install the MEGAHIT assembler

The [Trinity assembler](https://www.ncbi.nlm.nih.gov/pubmed/21572440) can also install it through `conda`:

```
conda install -y trinity 
```

## Change to a new working directory and link the original data

We will be using the same data as before ([Schurch et al, 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/)), so the following commands will create a new folder `assembly` and link the trimmed data we prepared earlier in the newly created folder:

```
cd ~/
mkdir -p assembly
cd assembly

ln -fs ~/quality/*.qc.fq.gz .
ls
```

## Applying Digital Normalization

In this section, we’ll apply [digital normalization](http://arxiv.org/abs/1203.4802) and [variable-coverage k-mer abundance trimming](https://peerj.com/preprints/890/) to the reads prior to assembly. This has the effect of reducing the computational cost of assembly [without negatively affecting the quality of the assembly](https://peerj.com/preprints/505/). Although the appropriate approach would be to use all 6 samples, for time consideration we will be using just the first one, i.e. ERR458493.qc.fq.gz

```
normalize-by-median.py --ksize 20 --cutoff 20 -M 10G --savegraph normC20k20.ct --force_single ERR458493.qc.fq.gz
```

This tools works mainly for paired-end reads, combined with a one file containing single-end reads. Given that all our samples are single end, we'll use the `--force_single` flag to force all reads to be considered as single-end. The parameter `--cutoff` indicates that when the median k-mer coverage level is above this number the read is not kept. Also note the `-M` parameter. This specifies how much memory diginorm should use, and should be less than the total memory on the computer you’re using. (See [choosing hash sizes for khmer for more information](http://khmer.readthedocs.io/en/v2.1.1/user/choosing-table-sizes.html).

(This step should take about 2-3 minutes to complete)


## Trim off likely erroneous k-mers

Now, run through all the reads and trim off low-abundance parts of high-coverage reads

```
filter-abund.py --threads 4 --variable-coverage --normalize-to 18 normC20k20.ct *.keep
```

The parameter `--variable-coverage` requests that only trim low-abundance k-mers from sequences that have high coverage. The parameter `--normalize-to` bases the variable-coverage cutoff on this median k-mer abundance.

(This step should take about 2-3 minutes to complete)

### Run the assembler


Trinity works both with paired-end reads as well as single-end reads (including simultaneously both types at the same time). In the general case, the paired-end files are defined as `--left left.fq` and `--right right.fq` respectively. The single-end reads (a.k.a _orphans_) are defined by the flag `--single`. 

First of all though, we need to make sure that there are no whitespaces in the header of the input fastq file. This is done using the following command:

```
cat ERR458493.qc.fq.gz.keep.abundfilt | tr -d ' ' > ERR458493.qc.fq.gz.keep.abundfilt.clean
```

So let's run the assembler as follows:

```
time Trinity --seqType fq --max_memory 10G --CPU 4 --single ERR458493.qc.fq.gz.keep.abundfilt.clean --output yeast_trinity
```

(This will take about 20 minutes)

You should see something like:

```
** Harvesting all assembled transcripts into a single multi-fasta file...

Saturday, June 30, 2018: 16:42:08       CMD: find /home/tx160085/assembly/yeast_trinity/read_partitions/ -name '*inity.fasta'  | /opt/miniconda/opt/trinity-2.6.6/util/support_scripts /partitioned_trinity_aggregator.pl TRINITY_DN > /home/tx160085/assembly/yeast_trinity/Trinity.fasta.tmp 
-relocating /home/tx160085/assembly/yeast_trinity/Trinity.fasta.tmp to /home/tx160085/assembly/yeast_trinity/Trinity.fasta
Saturday, June 30, 2018: 16:42:08       CMD: mv /home/tx160085/assembly/yeast_trinity/Trinity.fasta.tmp /home/tx160085/assembly/yeast_trinity/Trinity.fasta

###################################################################
Trinity assemblies are written to /home/tx160085/assembly/yeast_trinity/Trinity.fasta
###################################################################

Saturday, June 30, 2018: 16:42:08       CMD: /opt/miniconda/opt/trinity-2.6.6/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/tx160085/assembly/yeast_trinity/Trinity.fasta > /home/tx160085/assembly/yeast_trinity/Trinity.fasta.gene_trans_map
```

at the end.



### Looking at the assembly

First, save the assembly:

```
cp yeast_trinity/Trinity.fasta yeast-transcriptome-assembly.fa
``` 
 
Now, look at the beginning:

```
head yeast-transcriptome-assembly.fa
```
    
It's DNA! Yay!