# De novo transcriptome assembly

Learning objectives:

	* Learn what transcriptome assembly is
	* Learn to differentiate different types of assemblies
	* Discuss how do assemblers work
	* Learn to check the quality of a transcriptome assembly

What if you are working with a species with no existing reference genome or transcriptome...how do you construct a reference?

The answer is *de novo* assembly. The basic idea with *de novo* transcriptome assembly is you feed in your reads and you get out a bunch of *contigs* that represent transcripts, or stretches of RNA present in the reads that don't have any long repeats or much significant polymorphism. You run a  *de novo* transcriptome assembly program using the trimmed reads as input and get out a pile of assembled RNA.

Trinity, one of the leading *de novo* transcriptome assemblers, was developed at the [Broad Institute](http://www.broadinstitute.org/) and the [Hebrew University of Jerusalem](http://www.cs.huji.ac.il/). We will be losely following steps from the [Eel pond protocol](https://eel-pond.readthedocs.io/en/latest) for our guide to doing RNA-seq assembly.

## Download and trim the data

We will be using a set of *Nematostella vectensis* mRNAseq reads from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16).

To avoid needing to go though the trimming steps as before, we'll download a snakefile to take care of these steps for us. If you look at this file, you may notice it is very similar to the one you generated during the [snakemake challenge](./snakemake_for_qc.md).


Make sure you have snakemake installed in your base environment, where we will be running the snakefile:
```
conda install -c conda-forge snakemake-minimal
```

Download the snakefile and a corresponding conda environment file:
```
cd ~ # cd to home directory
curl -L https://osf.io/nqh6p/download -o nema_trim.snakefile
curl -L https://osf.io/xs6k7/download -o trim-env.yml
```

Let's take a look at the environment file to see what software snakemake will be putting in the environment it creates:

```
less trim-env.yml
```

Run the snakefile to download and trim the *Nematostella* reads:
```
snakemake -s nema_trim.snakefile --use-conda --cores 6
```

Here, we run snakemake. We use the `-s` command to tell snakemake where it can find our snakefile, 
tell it to build and use the environment above using the `--use-conda` flag, and have it execute 
the workflow over 6 cores using `--cores 6`.

The trimmed data should now be within the `nema_trimmed` folder.


## Generate a _de novo_ assembly

The following commands will create a new folder `assembly` and link the trimmed data we prepared in our 
snakemake workflow into the `assembly` folder:

```
cd 
mkdir assembly
cd assembly

ln -fs ../nema_trimmed/*.qc.fq.gz .
ls
```

Next we will combine our all forward reads into a single file and all reverse reads into a single file.
Usually, you want to have many samples from a single individual that are combined, but that minimize 
polymorphism and improve assembly. See this 
[preprint](https://www.biorxiv.org/content/10.1101/035642v2) for best practices for care and 
feeding of your transcriptome.

Use the following command to combine all fq into 2 files (left.fq and right.fq). 

```
cat *_1.pe.qc.fq.gz *se.qc.fq.gz > left.fq.gz
cat *_2.pe.qc.fq.gz > right.fq.gz
```

_Note: this step just makes it easier for us to type out the trinity command. Trinity can accept comma-separated lists of files as inputs - check the trinity documentation for details._


## Run the assembler

First, we'll create and activate an environment where Trinity is installed:

```
conda create -y -n trinity-env trinity
conda activate trinity-env
```

Trinity works both with paired-end reads as well as single-end reads (including with both types of reads at the same time). In the general case, the paired-end files are defined as `--left left.fq` and `--right right.fq` respectively. Our single-end reads (orphans) have been concatenated onto the `left.fq` file. 


So let's run the assembler as follows:

```
time Trinity --seqType fq --max_memory 16G --CPU 6 --left left.fq.gz --right right.fq.gz --output nema_trinity
```

(This will take a few minutes)

The `time` command allows us to see how long Trinity takes to run, but is not a part of the Trinity command.

You should see something like:

```
** Harvesting all assembled transcripts into a single multi-fasta file...

Thursday, July 11, 2019: 18:21:18       CMD: find /home/dibada/assembly/nema_trinity/read_partitions/ -name '*inity.fasta'  | /opt/miniconda/envs/trinity-env/opt/trinity-2.8.5/util/support_scripts/partitioned_trinity_aggregator.pl --token_prefix TRINITY_DN --output_prefix /home/dibada/assembly/nema_trinity/Trinity.tmp
-relocating Trinity.tmp.fasta to /home/dibada/assembly/nema_trinity/Trinity.fasta
Thursday, July 11, 2019: 18:21:18       CMD: mv Trinity.tmp.fasta /home/dibada/assembly/nema_trinity/Trinity.fasta


###################################################################
Trinity assemblies are written to /home/dibada/assembly/nema_trinity/Trinity.fasta
###################################################################


Thursday, July 11, 2019: 18:21:18       CMD: /opt/miniconda/envs/trinity-env/opt/trinity-2.8.5/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/dibada/assembly/nema_trinity/Trinity.fasta > /home/dibada/assembly/nema_trinity/Trinity.fasta.gene_trans_map
```

at the end.



## Looking at the assembly

First, save the assembly:

```
cp nema_trinity/Trinity.fasta nema-trinity.fa
``` 
 
Now, look at the beginning:

```
head nema-trinity.fa
```
    
These are the transcripts! Yay!

Let's capture also some statistics of the Trinity assembly. Trinity provides a handy tool to do exactly that:

```
TrinityStats.pl nema-trinity.fa
```

The output should look something like the following:

```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  18
Total trinity transcripts:      46
Percent GC: 46.68

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 3631
        Contig N20: 2497
        Contig N30: 2145
        Contig N40: 2097
        Contig N50: 1977

        Median contig length: 1593
        Average contig: 1459.98
        Total assembled bases: 67159


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 4447
        Contig N20: 3631
        Contig N30: 3631
        Contig N40: 2497
        Contig N50: 2151

        Median contig length: 1107
        Average contig: 1422.83
        Total assembled bases: 25611
```

This is a set of summary stats about your assembly. Are they good? Bad? How would you know?

Lastly, we'll deactivate the environment where Trinity is installed.

```
conda deactivate
```

## Suggestions for next steps 

After generating a *de novo* transcriptome assembly:
* [annotation](https://angus.readthedocs.io/en/2019/dammit_annotation.html)
* [evaluation](https://dibsi-rnaseq.readthedocs.io/en/latest/evaluation.html)
