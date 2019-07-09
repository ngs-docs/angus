# De novo genome assembly

<h1> TEMPORARY NOTICE: UNDER CONSTRUCTION!! (July-9-2019)</h1>

Here we're going to run through some of the typical steps for taking a newly sequenced bacterial isolate genome from raw fastq files through to an assembled, curated genome we can then begin to explore 🙂  

Before we get started, a public service announcement:

> **ATTENTION!**
> 
> This is not an authoritative or exhaustive workflow for working with a newly sequenced genome. No such thing exists! All genomes, datasets, and goals are different, and new tools are constantly being developed. The point of this page is just to give examples of some of the things we <i>can</i> do. <b>Don't let anything here, or anywhere, constrain your science to doing only what others have done!</b></div>

Now that that's out of the way, let's get to it!  

## Setting up our working environment
Throughout this process we'll be using a variety of tools. They are all going to be installed with the specific versions used at the time this was put together, using [conda](https://conda.io/docs/). First, we will create an environment, activate it, and then install all the needed programs.

```bash
conda create -y -n de_novo_example
conda activate de_novo_example

conda install -y -c bioconda -c conda-forge fastqc=0.11.5 \
              trimmomatic=0.36 spades=3.11.1 megahit=1.1.1 \
              quast=5.0.2 bowtie2=2.2.5 anvio=5.5.0 \
              centrifuge=1.0.4
```

> **CODE BREAKDOWN**
> 
> - **`conda create`** - creating our conda environment
>   - **`-y`** - say yes to all questions automatically (just do it and leave me alone)
>   - **`-n de_novo_example`** - names the environment we are creating "de_novo_example"
> 
> - **`conda activate de_novo_example`** - activates our new environment (puts us in it)
> 
> - **`conda install`** - installs the programs we need (in our current conda environment, which is "de_novo_example")
>   - **`-y`** - same as above, just do it, no questions
>   - **`-c bioconda`** and **`-c conda-forge`** - specifies to search these channels for the following programs
>   - **`fastqc=0.11.5`** and rest of positional arguments - names the program wanted, followed by an equals sign preceding the version we want

## The data
The practice data we're going to use here was provided by colleagues at the [J. Craig Venter Institute](http://www.jcvi.org/). In working out the details for a rather large-scale project, which in part involves sequencing a bunch of *Burkholderia* isolates from the ISS and performing de novo genome assemblies, [Aubrie O'Rourke](https://www.jcvi.org/about/aorourke) and her team put an already sequenced isolate – [*Burkholderia cepacia* (ATCC 25416)](https://www.atcc.org/products/all/25416.aspx) – through their pipeline in order to test things out and to have something to benchmark their expectations against. The sequencing was done on Illumina's Nextseq platform as paired-end 2x150 bps, with about a 350-bp insert size.  

The data files are rather large for this example, at 250 MB compressed, and 1 GB uncompressed. The download below comes with starting and output files, so that we can skip some of the more time-consuming steps and just copy over results. 

Copying and pasting this code block will get us the data:

```bash
cd ~
curl -L https://ndownloader.figshare.com/files/15458522 -o genomics_de_novo_temp.tar.gz
tar -xzvf genomics_de_novo_temp.tar.gz
rm genomics_de_novo_temp.tar.gz

cd genomics_de_novo_temp/
```

> **CODE BREAKDOWN**
> 
> - **`cd ~`** - this moves us into our home directory in case we aren't already there, so we are all starting in the same place
> 
> - **`curl`** - this let's us download the data
>   - **`-L`** - this specifies to follow "soft links" if there are any (these are like redirects when we click to download a file and are sent somewhere else first)
>   - **`https:...`** - this is the link to download the data
>   - **`-o`** - let's us specify the name of the output file, which we are naming "genomics_de_novo_temp.tar.gz"
> 
> - **`tar`** - program that unpacks "tarballs" (.tar) files, and can also optionally decompress them if they are gzipped (.gz)
>   - **`-x`** - exctract files from a .tar
>   - **`-z`** - decompress if gzipped (.gz)
>   - **`-v`** - print verbose output
>   - **`-f genomics_de_novo_tmp.tar.gz`** - specifies the file we want to act on
> 
> Then we delete the downloaded file with **`rm`** and change into our new directory with the last **`cd`** command

There are two subdirectories here, we are going to change into the "working" directory:

```bash
cd working/
```

## Quality trimming/filtering
Assessing the quality of our sequence data and filtering appropriately should pretty much always be one of the first things we do with our data. A great tool we've already used for this is [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Here we'll start with that on the raw reads.

## FastQC
Remember that [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) scans the fastq files we give it to generate a broad overview of some summary statistics. And it has several screening modules that test for some commonly occurring problems. But as we talked about before, its modules are expecting random sequence data, and any warning or failure notices the program generates should be interpreted within the context of our experiment (e.g. high duplication levels are to be expected with RNAseq or amplicon data, but FastQC would still flag them as something possibly wrong). Here we'll provide our forward and reverse reads as 2 positional arguments:

```bash
fastqc B_cepacia_raw_R1.fastq.gz B_cepacia_raw_R2.fastq.gz -t 4
```

Looking at our output from the forward reads (B_cepacia_raw_R1_fastqc.html), not too much stands out other than the quality scores are pretty mediocre from about 2/3 of the way through the read on: 

<center><img src="_static/fastqc_R1_before.png" width="90%"></center>  
<br>

> Here: the read length is stretched across the x-axis; the blue line is the mean quality score of all reads at the corresponding positions; red line is the median; the yellow boxplots represent the interquartile range; and the whiskers represent the 10th and 90th percentiles. 

The reverse reads look very similar: 

<center><img src="_static/fastqc_R2_before.png" width="90%"></center>  
<br>

## Trimmomatic
Also as we've seen, [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is a pretty flexible tool that enables us to trim up our sequences based on several quality thresholds and some other metrics (like minimum length or removing adapters and such). The summary from [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) didn't look all that terrible other than semi-low quality scores towards the end of each read, so for a first pass let's run Trimmomatic with pretty generic, but stringent settings and see what we get:

```bash
trimmomatic PE B_cepacia_raw_R1.fastq.gz B_cepacia_raw_R2.fastq.gz \
            BCep_R1_paired.fastq.gz BCep_R1_unpaired.fastq.gz \
            BCep_R2_paired.fastq.gz BCep_R2_unpaired.fastq.gz \
            LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:151 \
            -threads 6
```

> **CODE BREAKDOWN**
> 
> - **`PE`** - specifies are data are paired-end
> - **`B_cepacia_raw_R1.fastq.gz`** and **`B_cepacia_raw_R1.fastq.gz`** - first and second positional arguments are our input forward and reverse reads
> - **`BCep_R1_paired.fastq.gz `** and **`BCep_R1_unpaired.fastq.gz `** - third and fourth positional arguments are our output files for our forward reads: the third is for those that remain paired (i.e. their partner in the reverse files also passed the filtering thresholds); the fourth is for those that passed filtering but their partner did not 🥺 (gah, that emoji really tugs at the heart strings!)
> - **`BCep_R2_paired.fastq.gz `** and **`BCep_R2_unpaired.fastq.gz `** - fifth and sixth positional arguments are the same as 3rd and 4th above, but now for the reverse reads
> - **`LEADING:20`** - this states, starting from the left, if a base has a quality score less than 20, cut it off (will stop at first base above 20)
> - **`TRAILING:20`** - same as above, except working from the end of the read in
> - **`SLIDINGWINDOW:5:20`** - starting at base 1, look at a window of 5 bps, and if the average quality score of those 5 cut the read at the first position of that window and only keep up to that point
> - **`MINLEN:151`** - after performing all steps above, if the read is now shorter than 151, throw it away. For our purposes this means if any trimming happened the read will be removed. - Note that not all things happen in the order in which we provide them at the command line, but Trimmomatic does work that way 
> - **`-threads 6`** - specifies how many threads will be used, here 6 (["threads"](https://en.wikipedia.org/wiki/Thread_(computing)) are complicated)

The syntax for how to run Trimmomatic can be found in their [manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf), but our filtering thresholds here start with "LEADING:10". This says cut the bases off the start of the read if their quality score is below 10, and we have the same set for the end with "TRAILING:10". Then the sliding window parameters are 5 followed by 20, which means starting at base 1, look at a window of 5 bps and if the average quality score drops before 20, truncate the read at that position and only keep up to that point. The stringent part comes in with the MINLEN:151 at the end. Since the reads are already only 151 bps long, this means if any part of the read is truncated due to those quality metrics set above the entire read will be thrown away.  

The output from that shows us that only about 14% of the read pairs (both forward and reverse from the same fragment) passed, leaving us with only ~600,000 (just over ~300,000) read pairs. That sounds low, but since we know what we're working with here (meaning it's an isolate of a known genus and not a metagenome or something completely unknown), we can pretty quickly estimate if this could even possibly be enough coverage for us to assemble the genome we're expecting. Assuming those reads were perfect quality and perfectly evenly distributed (which they're not), that would be (350,000 read pairs) (600,000 paired reads) * (302 bps per paired read) = 181.2 Mbps covered. Most *Burkholderia* are around 8.5 Mbps, meaning we'd have around 20X (12X) coverage right now, if all were perfect. This confirms that this is a little low and we should probably adjust our stringency on filtering – I don't think there are solid rules on this either that always hold true, but in my (albeit limited) experience ~50–100X coverage is more around where you want to be for de novo assembly of a typical prokaryotic genome.  

So, I went back and altered how I was filtering a bit. Since the worst part of the reads quality-wise is at the end, and this is likely getting chopped off in a lot of reads during the sliding window scan, I next decided to try allowing a shorter minimum length of 141 to retain some of those reads: 

```bash
trimmomatic PE B_cepacia_raw_R1.fastq.gz B_cepacia_raw_R2.fastq.gz \
            BCep_R1_paired.fastq.gz BCep_R1_unpaired.fastq.gz \
            BCep_R2_paired.fastq.gz BCep_R2_unpaired.fastq.gz \
            LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 \
            MINLEN:140 -threads 4
```

These settings allowed ~49% of paired reads through (~1.2 million pairs) (~1.5 million pairs), which by the same quick estimation we did above suggests we could possibly have around (40-50X coverage) 50X coverage. Throwing away half the data here may feel uncomfortable, but it's always important to think about what we're doing with it when making these decisions. We are trying to assemble a good genome here, if we can do it with less "good data", that's better than doing it with more "bad data" – especially when bad sequence data could severely inhibit our assembly efforts.  

> **Keep in mind**
> 
> We're going to move forward with this here, but keep in mind this process doesn't need to be linear. While we can't try *everything* in the world trying to get the best out of our data, sometimes we can try *a lot* of things (resources and time permitting, of course). So just like we're going to run a few different assemblers below with different settings and compare them, we could just as easily test different assemblers with different quality-filtered data going into them.  

And just for a peek at the FastQC output after our trimming:

```bash
fastqc BCep_R1_paired.fastq.gz BCep_R2_paired.fastq.gz -t 4
```

<center><img src="_static/fastqc_after.png" width="90%"></center>
<br>

These are the forward reads, the reverse looked basically the same. Things still don't look perfect, but they look much cleaner than before – now our interquartile boxes (yellow) are much more snuggly sitting up top telling us our distribution of higher qualities across the end of the reads is much better. And though they weren't much of a factor here, don't forget to keep an eye on all the other modules from FastQC with any data you throw into it. They are not designed to be perfect assessments, because different experimental conditions can lead to warnings or failures for expected reasons as mentioned above, but they are very useful in that they can point you toward potential problems you might not otherwise see. When something does catch your eye, open up the manual for that module [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) and look into what it can mean and what can cause it 🙂  

## Read-error correction
Even though we've done our best to quality trim and filter our reads, there can still be errors in there that may hurt our assembly attempts. This paper by [Heydari et al.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1784-8) goes a bit into trying to systematically evaluate Illumina read-error correction tools, which might be a good place to start if you'd like to look into it more. In my experience it seems to consistently improve assembly, but as with most of these things, that isn't *always* the case – which is discussed a bit in that paper. Here we'll be trying different 4 different assemblies all with error-corrected reads, but feel free to experiment with that too 🙂

The read-error correction I typically use is employed by the SPAdes assembler. So even when I end up using an assembly from another program, I usually run the error-correction step within SPAdes first, and take the error-corrected reads from that to be used with another assembler. If we run the SPAdes program with the `--error-correction-only` flag specified, it will stop after that point and not go into assembly. 

This was the most computationally intensive step, taking about 45 minutes with the resources employed in the following command (which is more than are available on our instances). So for the sake of this tutorial, instead of running this command (which is commented out below), we are just going to copy over the result files in the next code block:

```bash
# spades.py -1 BCep_R1_paired.fastq.gz -2 BCep_R2_paired.fastq.gz \
            -o spades_error_corrected_reads -t 50 -m 500 \
            --only-error-correction
```

> **CODE BREAKDOWN**
> 
> - **`-1`** and **`-1`** - specify the input forward and reverse reads
> - **`-o`** - specifies the output directory where the results will be stored
> - **`-t`** - specifies how many threads will be used (["threads"](https://en.wikipedia.org/wiki/Thread_(computing)) are complicated)
> - **`-m`** - specifies the maximum amount of memory to possibly be used in gigabytes; program would quit if it needed more than that (our instances have 16 GB available to them)
> - **`--only-error-correction`** - telling SPAdes to only runn the error-correction step, and not perform an assembly

And here we are going to copy over the result files:

```bash
cp ../downloaded_results/BCep_R?_QCd_err_cor.fastq.gz .
```

## Assembly
Now that we have our reads quality filtered and we have run an error-correction step we're ready to move on to assembling them! There are lots of assembly programs out there, and once again, there is no one-size-fits-all. The reason each new assembly paper shows it doing better than all the others is not because everyone is lying about things, it's just that the data have a lot to say about which assembler is going to work the "best" – and what's "best" is not really a straightforward criterion to shoot for either. 

Two of the most commonly used assemblers today are [SPAdes](http://cab.spbu.ru/software/spades/) and [MEGAHIT](https://github.com/voutcn/megahit). SPAdes uses much more memory than MEGAHIT, so it is often suitable for working with one or a few genomes (like from an isolate or enrichment culture). But if working with high-diversity metagenomic samples, sometimes the memory requirements for SPAdes get too high, and MEGAHIT (which uses much less memory) can handle things. That's *not* to say SPAdes is always better when both can be run, it's not necessarily. 

When working with a new dataset, it's not uncommon to generate a few assemblies testing different programs with different parameters and then compare the results to try to feel a little confident we are doing the best that can currently be done with the data.  

Here we'll run a couple with [SPAdes](http://cab.spbu.ru/software/spades/) and a couple with [MEGAHIT](https://github.com/voutcn/megahit).

> We're going to be comparing a total of 4 different assembly variations here, with each taking about 10-15 minutes to run. So for the sake of time while doing this tutorial, we will only run the first one together, and the rest we'll copy over the results like we did with the error-correction step above.

## SPAdes
As mentioned above, [SPAdes](http://cab.spbu.ru/software/spades/) tends to do well with isolate or enrichment cultures that aren't too diverse. Since this is an isolate culture that was sequenced, it should run just fine. We're going to try with default settings and then with an additional option set, for a total for 2 assemblies from SPAdes. 

### SPAdes with default settings
Let's first run a SPAdes assembly with default settings. Note that we are providing the `--only-assembler` flag because we already ran the read-error correction step above, and are providing the output from that here as the input: 

```bash
spades.py -1 BCep_R1_err_corr.fq.gz -2 BCep_R2_err_corr.fq.gz \
          -o spades-default/ -t 6 --only-assembler
```

While that's running, let's talk about some of the others 🙂

### SPAdes in careful mode
In the SPAdes documentation, they have [this note](http://cab.spbu.ru/files/release3.11.1/manual.html#sec3.4) suggesting specific assembly settings when using 2x150 paired-end Illumina data as we are here. Most of the suggestions there are defaults now, and were run with our command above, but one, running things in `--careful` mode, was not. This tries to find and fix mismatches after the assembly is finished. Here's how we would run one with careful mode on, though again we're skipping it for time and will copy over the results once our other assembly finishes:

```bash
     ## don't run this code block, we will copy over results ##
# spades.py -1 BCep_R1_paired.fastq.gz -2 BCep_R2_paired.fastq.gz \
            -o spades-careful-assembly -t 6 --only-assembler \
            --careful
```

Command to copy over results when our first assembly is done:

```bash
cp -r ../downloaded_results/assemblies/spades-careful-assembly/ spades-careful-assembly/
```

## MEGAHIT
Next we'll run some assemblies with [MEGAHIT](https://github.com/voutcn/megahit). 

### MEGAHIT with default settings
Here's how we could run one with default MEGAHIT settings:

```bash
     ## don't run this code block, we will copy over results ##
megahit -1 BCep_R1_paired.fastq.gz -2 BCep_R2_paired.fastq.gz \
        -o megahit-default-assembly -t 6
```

Command to copy over results when our first assembly is done:

```bash
cp -r ../downloaded_results/assemblies/megahit-default-assembly/ megahit-default-assembly/
```

### MEGAHIT adjusted min-count parameter
The [MEGAHIT documentation](https://github.com/voutcn/megahit/wiki) has an assembly tips page that notes [here](https://github.com/voutcn/megahit/wiki/Assembly-Tips#filtering-kmin1-mer) that the default settings are largely tuned for metagenomic assembly, and that for a generic assembly (like our case here with an isolate), it is suggested to set the `--min-count` parameter (which deals with the frequency of kmers and filtering out reads with unique kmers) to `3` when there is greater than 40X coverage. This is the case with what we estimate here, so let's also try a run with that set:

```bash
     ## don't run this code block, we will copy over results ##
# megahit -1 BCep_R1_paired.fastq.gz -2 BCep_R2_paired.fastq.gz \
          -o megahit-min-count-3-assembly/ -t 6 --min-count 3
```

Command to copy over results when our first assembly is done:

```bash
cp -r ../downloaded_results/assemblies/megahit-min-count-3-assembly/ megahit-min-count-3-assembly/
```

Now that we have a handful of assemblies done, let's see how they compare.  

## Comparing assemblies
Let's just get this out there right off the bat, **there is no individual metric that exists to determine if we have a good assembly or not**, especially if we have no reference, and *especially* especially if we're working with a metagenomic assembly. 

There are, however, some general things we can look at, like N50 (half of the total assembly size can be found in contigs >= this size) or largest contig, or fraction of reads that successfully recruit to your assembly, or how many genes can we identify, and more. But: 1) these don't really have any context to tell us if they're "good" or not unless we're comparing multiple assemblies of the same data; and 2) we can have metagenomic assemblies with "worse" overall summary statistics compared to some other assembly, but that may actually enable us to do whatever it is we want to do *better* than the assembly with more appealing summary metrics. So, as usual, we have to keep in mind what our goals are, and know that picking the "best" assembly we can get out of our data is not necessarily a trivial or straightforward task. Having a reference genome like we do in this case makes things a lot easier as we'll see, and we'll talk about what we could do if we didn't have a reference in a bit 🙂


### QUAST
[QUAST](https://github.com/ablab/quast) is a really nice tool for comparing multiple assemblies, and for metagenome assemblies there is a comparable [MetaQUAST](http://bioinf.spbau.ru/metaquast). We can provide QUAST with all of our assemblies, a fasta file of our reference genome, and a .gff (**g**eneral **f**eature **f**ormat) file of our reference genome that contains information about its genes. Then QUAST will compare compare our assemblies to the reference and generate several reference-independent summary statistics. Our two reference files for our *Bulkholderia cepacia* ATCC 25416 came from NCBI [here](https://www.ncbi.nlm.nih.gov/genome/10703?genome_assembly_id=255013). And here's how we can run QUAST:  

```bash
quast -o quast-B-cep-out -r ../reference_genome/BCep_ref.fna \
      -g ../reference_genome/BCep_ref.gff -t 6 -m 1000 \
      -l "SPAdes-default, SPAdes-careful, MEGAHIT-default, MEGAHIT-min-count-3" \
      spades-default-with-err-corr-reads/contigs.fasta \
      spades-careful-with-err-corr-reads/contigs.fasta \
      megahit-default-with-err-corr-reads-assembly/final.contigs.fa \
      megahit-min-count-3-with-err-corr-reads-assembly/final.contigs.fa
```

> **CODE BREAKDOWN**
> 
> Remember the backslashes (**`\`**) don't have anything to do with the actual command, they are just there so this is not one long line and we can copy and paste (otherwise the "enter" character that comes after them would try to execute the command on each line.
> 
> - **`-o quast-B-cep-out`** - setting the output directory
> - **`-r reference_genome/BCep_ref.fna`** - providing a fasta file of our reference genome
> - **`-g reference_genome/BCep_ref.gff`** - providing the gene coordinates of our reference (this is so QUAST can see if it finds them in our assemblies)
> - **`-l "..."`** - within these quotes are labels for each of the input assemblies, provided here in the same order in which their input fasta files are entered next as positional arguments
> - the end is 4 positional arguments for our 4 input assemblies

There is more information about QUAST [in its documentation here](http://quast.bioinf.spbau.ru/manual.html). The output directory contains lots of information, but it is also nicely summarized in an html file called "report.html", which we can open and look at. Here's a portion of it:

<center><img src="_static/quast_output.png" width="95%"></center>
<br>

The columns here hold information about each of our 4 assemblies and the rows are different metrics. The majority of the rows starting from the top are in relation to the reference we provided, then the last few starting with "# contigs" are reference-independent. In the interactive html page, you can highlight the row names to get some help on what they mean, and there is more info in the [QUAST manual](http://quast.bioinf.spbau.ru/manual.html). The cells are shaded across the assemblies for each row from red to blue, indicating "worst" to "best", but this is only a loose guide to help your eye, as differences can be negligible or up to interpretation, and some rows are more important than others.  

It seems all of them did pretty well. They each reconstructed over 98% of the reference genome, which is pretty stellar, but none of them got 100%. This is most likely because short Illumina reads alone aren't able to assemble repetitive regions that extend longer than the paired-read fragment length. We can also see that our assemblies aligned across about 8,180 genomic features (genes, transcripts, tRNAs, and such) out of the 8,418 that are annotated in the reference genome, and looking at mismatches per 100 kbp we can see we're down around 2-4 mismatches per 100 kbp. 


<blockquote>
<center><b>QUICK QUESTION!</b></center>

What might be causing the 2-4 mismatches per 100 kbps?

<div class="toggle-header closed">
    <strong>Solution</strong>
</div>

<div class="toggle-content docutils container" style="width:100%">

This could be a mix of sequencing and/or assembly error, but it is also very likely some actual biological variation too. There is no reason to expect the genome to be <i>exactly</i> identical to the reference even after just 1 culturing transfer.

</div>
</blockquote>

Moving down to the reference-independent section in the table we can see which assemblies cover the most bps with the fewest contigs. This doesn't mean everything either, but fewer contigs covering about the same amount of bps does provide more information on synteny, which can be very important for trying to understand functional annotations. 

### Read recruitment
Another metric to assess the quality of an assembly is to see how well the reads that went into the assembly recruit back to it. It's sort of a way of checking to see how much of the data that went in actually ended up getting used. This tends to be more informative with mixed community metagenomes, if one assembly allows recruitment of 85% of our input reads, and another allows only 10% of the reads to recruit, it is safe to say the first one does a better job of representing what was recovered with our metagenomic sequencing. But in this case –  with an isolate genome when the few assemblies tested performed pretty similarly – it turns out that they all pretty much recruited reads just as efficiently (with the default settings of [bowtie2 v2.2.5](https://github.com/BenLangmead/bowtie2) at least). So I'm leaving the mention of this in here because it can be helpful with some datasets.  

As far as selecting which of our assemblies to move forward with here, since they all performed reasonably well, mostly for the reason of maximizing insight into synteny as mentioned above, based on the QUAST output we'll move forward with the SPAdes assembly done with the `--careful` setting as that gave us the most genes on the fewest contigs.


## What about when we don't have a reference?
It's pretty nice here because we have a known reference genome. When that's not that case, it's much harder to be confident we are making the best decisions we can. With no reference, we might go further down the path of processing and analysis with multiple assemblies. Some might be better than others depending on our particular purpose, and that might not always be apparent based on overall summary statistics.  

For instance, if recovering representative genomes from a metagenomes is the goal, we might want to start that process and see which assembly is better suited for that. As mentioned above, in some cases it is possible to be able to better recover representative genomes from an assembly that has seemingly worse overall summary statistics. And if we just chose based on the assembly-level summary statistics, we wouldn't find that out. Or if our question is more about the functional potential of the whole community, or maybe looking for specific genes in a metagenomic assembly, then seeing which assembly gives more [open-reading frames](https://en.wikipedia.org/wiki/Open_reading_frame) or better annotation results might help us decide. The bottom line is it's difficult, and there is no one answer, but we do the best we can 🙂  

## Exploring our assembly with anvi'o
Now that we've selected the assembly we're going to move forward with, we can start to take a deeper look at it. To do that we're going to use a wonderful platform [anvi'o](http://merenlab.org/software/anvio/) - which stands for analysis and visualization of omics data.

anvi'o is very expansive, but to try to summarize it in a sentence, it is a powerful and user-friendly data visualization and exploration platform. There is a ton of excellent documentation at the [anvi'o site](http://merenlab.org/software/anvio/), including workflows, tutorials, and blog posts. 

Here we're going to put our isolate-genome assembly into the anvi'o framework and see just a few of the ways it can help us begin to look at our assembled genome - including a visual way of checking for possible contamination (which in this case would be anything in there that's *not* our target *Burkholderia cepacia*). Many of the major steps we're going to be performing here are laid out in the [metagenomic workflow presented here](http://merenlab.org/2016/06/22/anvio-tutorial-v2/), as much of the processing is the same even though we are not working with a metagenome here. 

### Putting our assembly into anvi'o and generating some information about it
For us to get our assembly into anvi'o, first we need to generate what it calls a [contigs database](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#creating-an-anvio-contigs-database). This contains the contigs from our assembly and some information intrisic to them (like GC content). The following command will organize our contigs in an anvi'o-friendly way, generate some basic stats about them, and use the program [Prodigal](https://github.com/hyattpd/Prodigal) to identify [open-reading frames](https://en.wikipedia.org/wiki/Open_reading_frame). This will just take a few minutes:

```bash
anvi-gen-contigs-database -f spades-careful-assembly/contigs.fasta \
                          -o contigs.db -n B_cepacia_assembly
```

> **CODE BREAKDOWN**
> 
> - **`-f`** - specifies our input fasta file of our chosen assembly
> - **`-o`** - specifies our output file name, which has the ".db" extension as preferred by anvi'o
> - **`-n`** - specifies a name for our anvi'o project

Now that we have our `contigs.db` that holds our sequences and some basic information about them, we can start adding more information. This is one of the places where the flexibility comes into play, but for now we'll just move forward with some parts of a general anvi'o workflow, including:

- using the program [HMMER](http://hmmer.org/) with profile hidden Markov models (HMMs) to scan for bacterial single-copy genes [(from Campbell et al. 2013)](http://www.pnas.org/content/110/14/5540.short)
  - if new to HMMs, see the bottom of page 7 [here](http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf) for a good explanation of what exactly a "hidden Markov model" is in the realm of sequence data
  - to summarize, instead of treating all columns of an alignment equally (like something like BLAST does), HMMs can weight individual columns in an alignment differently. So the highly-conserved stretches within, say, a conserved domain of a protein will matter more than a stretch in between two domains that may vary more.

- using [NCBI COGs](https://www.ncbi.nlm.nih.gov/COG/) to functionally annotate the open-reading frames [Prodigal](https://github.com/hyattpd/Prodigal) predicted

- using a tool called [centrifuge](https://ccb.jhu.edu/software/centrifuge/index.shtml) for taxonomic classification of the identified open-reading frames

- generate coverage information for our assembly so we can visualize it

The following block of code took ~45 minutes to run, so we are going to breakdown , mostly because of needing to download and setup the COG and centrifuge databases for the first time, the following code block doing the mapping took about 10 minutes, and then the `anvi-profile` step took about 10 minutes. So feel free to skip these next three code blocks and copy over the result files as shown after them to move forward if you'd like 🙂

```bash
  # HMM searching for bacterial single-copy genes
anvi-run-hmms -I Campbell_et_al -c contigs.db -T 6 # took < 1 min.

  # functional annotation with DIAMOND against NCBI's COGs
anvi-setup-ncbi-cogs -T 6 # only needed the first time, took ~5 min.
anvi-run-ncbi-cogs -c contigs.db --num-threads 6 # took ~ 1 min.

  # exporting Prodigal-identified open-reading frames from anvi'o
anvi-get-sequences-for-gene-calls -c contigs.db -o gene_calls.fa

  # setting up centrifuge for taxonomy (only needed the first time)
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz
tar -xzvf p_compressed+h+v.tar.gz && rm -rf p_compressed+h+v.tar.gz

  # running centrifuge taxonomy (this taxonomically classifies each identified coding sequence)
centrifuge -f -x p_compressed+h+v gene_calls.fa -S centrifuge_hits.tsv -p 6 # took < 1 min.

  # importing the taxonomy results into our anvi'o contigs database
anvi-import-taxonomy-for-genes -c contigs.db -i centrifuge_report.tsv centrifuge_hits.tsv -p centrifuge
```

The last thing we want to add right now is mapping information from recruiting our reads to the assembly. Here is generating that with [bowtie2](https://github.com/BenLangmead/bowtie2) and preparing for anvi'o:

```bash
  # building bowtie index from our selected assembly fasta file
bowtie2-build spades-careful-assembly/contigs.fasta \
              spades-careful-assembly.btindex

  # mapping our reads (takes ~1 min.)
bowtie2 -q -x spades-careful-assembly.btindex \
        -1 BCep_R1_err_corr.fq.gz -2 BCep_R2_err_corr.fq.gz \
        -p 6 -S spades-careful-assembly.sam

  # converting to a bam file (takes ~1 min.)
samtools view -bS spades-careful-assembly.sam > B-cep-assembly.bam

  # sorting and indexing our bam file (takes ~1 min.)
anvi-init-bam B-cep-assembly.bam -o B-cep.bam
```

We can then integrate this mapping information into anvi'o with the `anvi-profile` program, which generates another type of database anvi'o calls a "profile database". For our purposes right now, this is mostly about storing coverage information for each contig. We'll see what that looks like and why it's helpful soon!

```bash
  # integrating the read-mapping info into anvi'o (takes ~2 min.)
anvi-profile -i B-cep.bam -c contigs.db -M 1000 -T 6 --cluster-contigs -o B-cep-profiled/
```

If skipping the last 3 code blocks for time, you can copy over the appropriate files with the following:

```
cp ../downloaded_results/anvio_files/contigs.db .
cp -r ../downloaded_results/anvio_files/B-cep-profiled/ B-cep-profiled/
```

### Pulling out some sequences and summarizing our assembly
Ok, great, so we've just generated and put quite a bit of information about our assembly into our contigs and profile databases. And this being anvi'o, it's now very easy to access specifics of that information when we want it (you can see all of the available programs by typing `anvi-` and hitting tab twice).

For example, one of the commands we just ran, `anvi-run-hmms`, searched for bacterial single-copy genes. We can ask anvi'o to give us which genes were identified using `anvi-get-sequences-for-hmm-hits` (and adding the flag `-h` will tell us all about how to use the program). Including adding the `-l` flag to see which HMM sets we have :

```bash
anvi-get-sequences-for-hmm-hits -c contigs.db -l
```

Which shows us "Campbell_et_al" which we specified above. And adding the `-L` flag will list the names of the genes in that set:

```bash
anvi-get-sequences-for-hmm-hits -c contigs.db --hmm-sources Campbell_et_al -L
```

And say we want to pull out the "RecA" gene, which is a commonly targeted bacterial single-copy gene involved in DNA repair:

```bash
anvi-get-sequences-for-hmm-hits -c contigs.db --hmm-sources Campbell_et_al --gene-names RecA --get-aa-sequences --no-wrap -o RecA.faa
```

If we print out what's in that file with `cat`:

```bash
cat RecA.faa
```

Copy it and paste it into [NCBI's protein BLAST here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and BLAST it, we can see that the top hit is 100% identical across the entire query to *B. cepacia*, and all top hits are to the same species (which is a nice little sanity check). Here are just a few from the top:

<center><img src="_static/blastp-B-cep.png" width="95%"></center>  
<br>

We can also generate some summary statistics on our assembly, including estimated percent completion and redundancy based on the presence/absence of the single-copy marker genes we scanned for above. Here are the two steps needed to summarize our isolate-genome assembly: 

```bash
  # this is adding all contigs to a group called "DEFAULT"
anvi-script-add-default-collection -p B-cep-profiled/PROFILE.db
  # and here is our summary command
anvi-summarize -c contigs.db -p B-cep-profiled/PROFILE.db \
               -C DEFAULT -o B-cepacia-assembly-summary/
```

A lot was generated with that, now found in our new directory, "B_cepacia_assembly_summary/". If we glance at the "bins_summary.txt" file, we can see some summary statistics of our assembled genome, including the completion/redundancy estimates (`column -t` will keep the columns lined up for us):

```bash
column -t B-cepacia-assembly-summary/bins_summary.txt
```

```bash
## output from above command ##
# bins        taxon         total_length  num_contigs  N50     GC_content         percent_completion  percent_redundancy
# EVERYTHING  Burkholderia  8472061       92           194480  66.44548680306582  99.28057553956835   1.4388489208633093
```

Which shows us, in the second column, the majority of taxonomy calls by [Centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml#obtaining-centrifuge) were actually for the right genus, that's always nice. Also, based on the [Campbell et al. 2013](http://www.pnas.org/content/110/14/5540.short) bacterial single-copy marker genes, our assembly is estimated to be ~99.3% complete with ~1.4% redundancy. This approach doesn't actually add up to 100% completion and 0% redundancy in many, if any, organisms, so for comparison's sake, running the ATCC 25416 reference genome through the same pipeline gives the same results based on this single-copy gene set. Though we are still ~130 Mbps short, as we saw with the QUAST output above.

But great, things look pretty much as they should so far. We can also visually inspect our assembly, and how the reads that went into it recruit to it. In theory, if all the DNA in the assembly came from the same organisms (i.e. it's a clean assembly), there should be pretty even coverage across the whole thing. And if there is contamination in there (something other than our target *B. cepacia*, it will likely have very different coverage). So let's finally take a look with `anvi-interactive`. 

### Visualizing with anvi-interactive
First let's copy-and-paste this command to get the link we will need to view our anvi'o interactive interface:

```bash
echo "http://$(hostname -i):8080"
```

Now let's start up the anvi'o interactive interface:

```bash
anvi-interactive -c contigs.db -p B-cep-profiled/PROFILE.db \
                 --server-only -P 8080
```

Now if we go back to our browser, and refresh or load the address we copied above, it should open to the anvi'o interactive interface and look something like this: 

<center><img src="_static/fresh_anvi.png" width="90%"></center>
<br>

So there is a lot going on here at first glance, especially if you're not yet familiar with how anvi'o organizes things. So here's a quick crash course 🙂  

At the center of the figure is a hierarchical clustering of the contigs from our assembly (here clustered based on tetranucleotide frequency and coverage). So each tip (leaf) represents a contig (or a fragment of a contig as each is actually broken down into a max of ~20,000bps, but for right now we'll just be referring to them as contigs). Then radiating out from the center are layers of information ("Parent", "Length", "GC content", etc.), with each layer displaying information for each contig.

One of the first things that might jump out here is the outer layer (purple in this image, but may be different on yours), labeled "Taxonomy". There is actually a color for each contig for whatever taxonomy was assigned to the majority of genes in that particular contig. This solid bar all around tells us that the genes in almost the entire assembly were identified as *Burkholderia* – minus two contigs which were not classified as anything. The next thing that stands out is how stable the mean coverage is across all contigs, other than mostly just that same area on the right side. That likely holds things found in multiple copies in the genome. For instance, according to [IMG](https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=2509276048), *B. cepacia* ATCC 25416 has 7 16S copies and 9 23S copies, which would complicate assembly if they aren't all identical, and would inflate their coverage compared to the parts of the genome that exist in single copy. Overall this is great and shows the culture really seems to have been axenic.  

Just for a quick comparison, here is the same type of figure (taxonomy at the inner-part of the circle here), but from an enrichment culture, rather than axenic:  

<center><img src="_static/other_anvi.png" width="90%"></center>
<br>

Here the highlighted contigs, labeled "Bin 1", represent the targeted cultivar from this sequencing run, demonstrating a nice example of how anvi'o can help you manually curate bins you're trying derive from assemblies.  

While we didn't need much (any) manual curation in this case, it was still a good idea to visually inspect the coverage of our assembly to make sure nothing weird was going on. And if we wanted, we could further explore those parts with higher coverage to find out which parts of the genome seem to exist in greater than 1 copy.  


## Epilogue
This was all basically to get a high-quality draft of our isolate genome, that we could feel confident about investigating further. Once we feel comfortable with our assembled genome, we can go a lot of different ways. Going into individual approaches are beyond the scope of this particular page, but here are just a few examples.  


### Phylogenomics
Pull available reference genomes of close relatives and build a phylogenomic tree to get a robust estimate of where our newly acquired isolate(s) fit in evolutionarily with what is already known. 

<center><img src="_static/phylo_syn.png" width="90%"></center>
<br>

### Distributions
Pull available metagenomes from studies and recruit the reads to a reference library containing our isolate (and its close relatives if it has any) to begin assessing their distributions across samples. This example is done with ocean samples, but the same principle can be applied to any environments.

<center><img src="_static/dist_syn.png" width="90%"></center> 
<br>

### Pangenomics
Start investigating differences in the genetic complement of our new isolate as compared to its known close relatives. This example figure is combining pangenomics (the core of the figure showing the presence or absence of genes within each genome) with metagenomics (distributions of the genomes across samples in the top right corner) to try to associate genomic variability with ecological delineations:

<center><img src="_static/pan_syn.png" width="90%"></center>  
<br>

And so much more! This is where just pure data crunching slows down, and the actual science begins. The above are just some of the ways to get to the point where we can then consider our experimental design and our questions and let them guide where we go next 🙂
