# Transcriptome assembly - some basics

Learning objectives:

	* Learn what transcriptome assembly is
	* Learn to differentiate different types of assemblies
	* Discuss how do assemblers work
	* Learn to check the quality of a transcriptome assembly
	
In [variant calling](http://angus.readthedocs.io/en/2019/mapping-variant-calling.html), we mapped reads to a reference and looked systematically for differences.

But what if you don't have a reference? How do you construct one?

The answer is *de novo* assembly, and the basic idea is you feed in your reads and you get out a bunch of *contigs*, that represent stretches of RNA present in the reads that don't have any long repeats or much significant polymorphism.  Like everything else, the basic idea is that you run a program, feed in the reads, and get out a pile of assembled RNA.

Trinity, developed at the [Broad Institute](http://www.broadinstitute.org/) and the [Hebrew University of Jerusalem](http://www.cs.huji.ac.il/, represents a novel method for the efficient and robust de novo reconstruction of transcriptomes from RNA-seq data. We will be using the [eel-pond protocol](https://eel-pond.readthedocs.io/en/latest), which uses trinity, as our guide to doing RNA-seq assembly.

## Boot up a Jetstream

[Boot an m1.medium Jetstream instance](jetstream/boot.md) and log in.

## Install the TRINITY assembler in a conda environment

Create and environment and install the [Trinity assembler](https://www.ncbi.nlm.nih.gov/pubmed/21572440) through `conda`:

```
conda create -y -n trinity-env trinity khmer
```
where: 
   - `-y` says to install the pkgs without double checking with me
   - `-n` names the environment "trinity-env"

Now, let's activate that environment:

```
conda activate trinity-env
```

## Make sure you've got the data

We will be using the same data as before ([Schurch et al, 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/)), so the following commands will create a new folder `assembly` and link the trimmed data we prepared earlier in the newly created folder:

```
cd ~
ls -l quality/*qc* .
```

You should see the following:

```
-rw-rw-r-- 1 dibada dibada  59465832 Jul 11 15:45 quality/ERR458493.qc.fq.gz
-rw-rw-r-- 1 dibada dibada  58503534 Jul 11 15:46 quality/ERR458494.qc.fq.gz
-rw-rw-r-- 1 dibada dibada  58054460 Jul 11 15:46 quality/ERR458495.qc.fq.gz
-rw-rw-r-- 1 dibada dibada 102164315 Jul 11 15:44 quality/ERR458500.qc.fq.gz
-rw-rw-r-- 1 dibada dibada 101187090 Jul 11 15:44 quality/ERR458501.qc.fq.gz
-rw-rw-r-- 1 dibada dibada 100550095 Jul 11 15:45 quality/ERR458502.qc.fq.gz
```

If you don't see the data, you can download it using the following commands:

```
cd ~
mkdir -p quality
cd quality
curl -L https://osf.io/wfz34/download -o ERR458493.qc.fq.gz
curl -L https://osf.io/jxh4d/download -o ERR458494.qc.fq.gz
curl -L https://osf.io/zx7n3/download -o ERR458495.qc.fq.gz
curl -L https://osf.io/96mrj/download -o ERR458500.qc.fq.gz
curl -L https://osf.io/wc8yn/download -o ERR458501.qc.fq.gz
curl -L https://osf.io/sdtz3/download -o ERR458502.qc.fq.gz
```


## Change to a new working directory and link the original data

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
filter-abund.py --threads 6 --variable-coverage --normalize-to 18 normC20k20.ct *.keep
```

The parameter `--variable-coverage` requests that only trim low-abundance k-mers from sequences that have high coverage. The parameter `--normalize-to` bases the variable-coverage cutoff on this median k-mer abundance.

(This step should take about 2-3 minutes to complete)

### Run the assembler


Trinity works both with paired-end reads as well as single-end reads (including simultaneously both types at the same time). In the general case, the paired-end files are defined as `--left left.fq` and `--right right.fq` respectively. If assembling from single-end reads (as in this case), we use the flag `--single`. 

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
All commands completed successfully. :-)



** Harvesting all assembled transcripts into a single multi-fasta file...

Thursday, July 11, 2019: 17:15:47	CMD: find /home/dibada/assembly/yeast_trinity/read_partitions/ -name '*inity.fasta'  | /opt/miniconda/envs/trinity-env/opt/trinity-2.8.5/util/support_scripts/partitioned_trinity_aggregator.pl --token_prefix TRINITY_DN --output_prefix /home/dibada/assembly/yeast_trinity/Trinity.tmp
-relocating Trinity.tmp.fasta to /home/dibada/assembly/yeast_trinity/Trinity.fasta
Thursday, July 11, 2019: 17:15:47	CMD: mv Trinity.tmp.fasta /home/dibada/assembly/yeast_trinity/Trinity.fasta


###################################################################
Trinity assemblies are written to /home/dibada/assembly/yeast_trinity/Trinity.fasta
###################################################################


Thursday, July 11, 2019: 17:15:47	CMD: /opt/miniconda/envs/trinity-env/opt/trinity-2.8.5/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/dibada/assembly/yeast_trinity/Trinity.fasta > /home/dibada/assembly/yeast_trinity/Trinity.fasta.gene_trans_map
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
    
It's RNA! Yay!

Let's capture also some statistics of the Trinity assembly. Trinity provides a handy tool to do exactly that:

```
TrinityStats.pl yeast-transcriptome-assembly.fa
```

The output should look something like the following:

```

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	3296
Total trinity transcripts:	3320
Percent GC: 42.02

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 1420
	Contig N20: 1069
	Contig N30: 802
	Contig N40: 634
	Contig N50: 514

	Median contig length: 321
	Average contig: 448.98
	Total assembled bases: 1490618


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 1385
	Contig N20: 1021
	Contig N30: 787
	Contig N40: 618
	Contig N50: 502

	Median contig length: 319
	Average contig: 442.56
	Total assembled bases: 1458676

```

This is a set of summary stats about your assembly. Are they good? Bad? How would you know?
