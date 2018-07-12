# Assembling a metagenome and recovering "genomes" with Anvi'o

Though it is only one of many approaches you can take when working with a metagenomics dataset ([e.g. here's a generic overview pdf that tries to highlight some](https://ndownloader.figshare.com/files/12367187)), recovering genomes from metagenomes has become a powerful tool for microbial ecologists. Here we will assemble a metagenome, and go through the process of "binning" our assembled contigs into groups based on coverage and sequence composition using the analysis and visualization platform [anvi'o](http://merenlab.org/software/anvio/).

> **NOTE:** Even the highest quality genomes recovered from metagenomes are not the same as isolate genomes. It depends on the data and assembly, but in general they are more of an agglomeration of very closely related organisms from the sample due to the assembly process and fine-scale variation that exists in microbial populations. They are being more frequently referred to as "metagenome-assembled genomes", or MAGs, to better convey this. 

---

Learning objectives:

	* Think about assembling individual samples vs co-assembling multiple samples together
	* Assemble a metagenome
	* Visually identify clusters of contigs ("binning")

## Boot up a Jetstream
[Boot an m1.medium Jetstream instance](https://angus.readthedocs.io/en/2018/jetstream/boot.html) and log in.

## Software needed
We will be using several different tools here, all are installable with `conda`. This line should do the trick:

## What is a co-assembly?
"Co-assembly" refers to performing an assembly where the input files would be reads from multiple samples. This is in contrast to doing an independent assembly for each sample, where the input for each would be just the reads from that individual sample. Three major benefits of co-assembly include: higher read depth (this *can* allow you to have a more robust assembly that captures more of the diversity in your system, *but not always*); it facilitates the comparison across samples by giving you one reference assembly to use for all; and it can substantially improve your ability to recover genomes from metagenomes due to the awesome power of differential coverage (you can download a slide showing how coverage is used to do this from here -> [keynote](https://ndownloader.figshare.com/files/12367211), [powerpoint](https://ndownloader.figshare.com/files/12367226).

## To co-assembly or not to co-assemble
Though a co-assembly has its benefits, it will not be ideal in all circumstances. And as with most things in bioinformatics, there are no golden rules as for when it would be better to co-assemble multiple samples over when it would be better to run individual assemblies on each. It could depend on many factors (e.g. variation between the samples, diversity of the communities, the assembler(s) you're trying, you're overall questions, who knows how many others?, etc.), and there can be datasets where a co-assembly would give you a poorer output assembly than an individual assembly would. So it is something you should consider with each dataset you work with, and if you are unsure and have the time/resources, then trying and assessing the outcome is always good practice. 

## Our practice data
To work with a smaller dataset here that will let us do things in a reasonable amount of time, we're going to be working with a relatively simple microbial community here that comes from metagenomic sequencing of an enrichment culture of the nitrogen-fixing cyanobacterium *Trichodesmium*. Metagenomics still takes a lot of time, so we're going to skip over quality trimming here (already done), though assessing the quality and trimming/filtering needed as laid out in [this lesson](https://angus.readthedocs.io/en/2018/quality-and-trimming.html) should pretty much always be the first step. There are still some steps that would take a bit too long to just wait for, so in those cases there will be examples of how the code would be run, but we're also going to download a results directory that we can pull from to skip some of the more time-consuming steps 🙂

You can download the entire working directory using the following command:

```bash
cd
curl -O -L 
tar -xzvf 
rm tar.gz
cd metagen_tut
```

This main directory we just changed into holds 2 subdirectories: "data" which holds our 4 samples' forward (R1) and reverse (R2) reads; "results" which holds our result files we'll use from time to time to skip longer steps; and "working" where we are going to be running our commands from. So let's take a look at it with `ls`, and then change into the "working" directory:

```bash
ls
cd working/
ls
```

I typically make a "samples.txt" file that contains each of my sample names when I start with a new project as I'll often use it for some loops as we'll see below. One way to do that would be like this:

```bash
ls ../data/*.gz | cut -f3 -d "/" | cut -d "_" -f1,2 | uniq > samples.txt
```

> **Code breakdown:** build this code up one command at a time (press enter for each part, then add the `|` and the next part)
> * `ls` is the starting command, and we're telling it to list all the files that end in ".gz" from our "data" directory, then we pipe `|` the output from that into the next command
> * `cut` pulls out specific columns. It uses tabs as default to be the character that separates columns, but we are adding a flag to modify that
>   * `-d` is how you specify the delimiter (what separates columns) for `cut`, here we are saying to use the underscore character `_`, just because that would allow us to split the file names up where we wanted (working to get just the sample names with no extensions or read direction info)
>   * `-f` stands for "field" (column), and this is where we tell `cut` that we want the first and second column (which based on our delimiter grabs "Sample_A" for example). Right now we'd have double lines for each sample, which isn't what we want in our file of sample names, so we're going to pipe `|` this output also into one last command
> * `uniq` removes duplicate lines and leaves only singles of each (**note** things need to be sorted for `uniq` to work, here that was already done, but you should look into it further in the future if unfamiliar with it)
> * `>` finally we redirect the final output to a file names "samples.txt"

And now we're ready to get to work!

## Co-assembly
Above, we briefly touched on some plusses and minuses of co-assembly vs individual-sample assembly. Here we are working with 4 samples from an enrichment culture, and the communities *should be* pretty well constrained - certainly relative to complex environmental samples like, say, a transect across soil. That would be harder to decide, but for us, it's a pretty safe start to go with a co-assembly. 

There are many assemblers out there, and the reason each has a paper showing how it beats others is because every dataset is different. Some assemblers work better for some datasets, and others work better for others. That's not to say all are magically equally good in every sense, but most that gather a following will out-perform all others under certain conditions. I typically try several (SPAdes, Megahit, idba-ud), and compare them with QUAST (for individual genome assembly) or MetaQUAST (for metagenome assemblies). Choosing the "best" is also not straightforward and it can depend on what you're doing, there is more on that [here](https://astrobiomike.github.io/genomics/de_novo_assembly#comparing-assemblies) if you're interested, along with testing other assemblies/options and comparing them.

In this case, here is how the command would be run with the Megahit assembler. But even with the smaller dataset we're using, it takes about 40+ minutes on our cloud instances using 4 cpus. So rather than running this, we're just going to pull the assembly output file from our results directory:

```bash
  ## don't run this, would take about 40+ minutes ##
# megahit -1 ../data/SRR3880207_nonTricho_1.fastq.gz,../data/SRR3880208_nonTricho_1.fastq.gz,../data/SRR3880209_nonTricho_1.fastq.gz,../data/SRR3880210_nonTricho_1.fastq.gz -2 ../data/SRR3880207_nonTricho_2.fastq.gz,../data/SRR3880208_nonTricho_2.fastq.gz,../data/SRR3880209_nonTricho_2.fastq.gz,../data/SRR3880210_nonTricho_2.fastq.gz -o megahit_default -t 4
```

> **Code breakdown:** 
> * `megahit` is the command, it is specifying the program we are using
>   * `-1` is for specifying the forward reads. Since we are doing a co-assembly of multiple samples, we could have stuck them altogether first into one file, but `megahit` allows multiple samples to be added at once. Notice here we are specifying all the forward reads in a list, each separated by a **comma**, with **no** spaces. This is because spaces are important to the command-line's interpretation, whereas commas don't interfere with how it is interpreting the command and the information gets to megahit successfully. 
>   * `-2` contains all the reverse read files, separated by commas, *in the same order*.
>   * `-o` specifies the output directory (the program will produce information in addition to the final output assembly)
>   * `-t` specifies that we want to use 4 cpus (some programs can be run in parallel, which means splitting the work up and running portions of it simultaneously)

After waiting for that to finish, we would have the "megahit_default" directory that is currently within our "results" directory. So we're just going to copy over that final output fasta file (our co-assembly) into our working directory:

```bash
cp ../results/megahit_default/final.contigs.fa .
```

If you glance at this file, with `head final.contigs.fa` for example, you can see there are spaces and some other special characters in the headers (the lines that start with ">"). This can be problematic for some tools, and it's better to have simplified headers if you can. There are lots of ways to make modifications like that, and some scripts already exist. Anvi'o happens to have one, so we're just going to use that here, and also for processing speed moving forward we're going to filter out any contigs that are shorter than 1,000 bps.

```bash
anvi-script-reformat-fasta final.contigs.fa -o contigs.fa -l 1000 --simplify-names
```

> **Code breakdown:** 
> * `anvi-script-reformat-fasta` is the command, it is specifying the program we are using
>   * the first "positional" argument (no flag needed) is specifying the input file (our assembly)
>   * `-o` specifies the name of the output file
>   * `-l` specifies we wanted to filter out sequences with fewer than 1000 bps
>   * `--simplify-names` tells the program we want it to give all sequences a "clean" header (no special characters)

## Mapping our reads to the assembly they built
Among other things (like enabling variant detection), mapping our reads for each sample to the co-assembly they built gives us "coverage" information for each contig in each sample, which as discussed above will help us with our efforts to recover metagenome-assembled genomes (MAGs). 

Here we are going to use [bowtie2]() to do our mapping, and first need to create an index of our co-assembly:

```bash
bowtie2-build contigs.fa assembly
```

And here is where we would map our individual samples' reads to our co-assembly, but this also would take about 25 minutes, so we'll look at the commands here, but don't run them. This can be done one at a time like this:

```bash
bowtie2 -x assembly -q -1 ../data/Sample_A_1.fastq.gz -2 ../data/Sample_A_2.fastq.gz --no-unal -p 4 -S Sample_A.sam
```

> **Code breakdown:** 
> * `bowtie2` is the command, it is specifying the program we are using
>   * `-x` is where we specify the "basename" of the index we just built, which is "assembly", `bowtie2` will find all the appropriate files using that basename 
>   * `-1` specifies the forward reads of the sample we're doing
>   * `-2` specifies the reverse reads
>   * `--no-unal` tells the program we don't want to record the reads that don't align successfully (saves on storage space)
>   * `-p` specifies how many cpus we want to use
>   * `-S` specifies the name of the output "sam" file we want to create (**S**equence **A**lignment **M**ap)


Or with a loop, like this:

```bash
for sample in $(cat samples.txt); do bowtie2 -x assembly -q -1 ../data/"$sample"_1.fastq.gz -2 ../data/"$sample"_2.fastq.gz --no-unal -p 4 -S "$sample".sam; done
```

> **Code breakdown:** 
> The special words of the loop are `for`, `in`, `do`, and `done`.
> * `for` – here we are specifying the variable name we want to use ("sample" here), this could be anything
> * `in` – here we are specifying what we are going to loop through, in this case it is every line of the "samples.txt" file
>   * `$(cat ...)` – this is a special type of notation in shell. The operation within the parentheses here will be performed and the output of that replaces the whole thing. It would be the same as if we typed out each line in the file with a space in between them
> * `do` – indicates we are about to state the command(s) we want to be performed on each item we are looping through
>   * `bowtie2` here the explanation is the same as above, except the part that has to change for each sample, the sample name. The loop gets this from the "samples.txt" file for each iteration, and we tell it were to swap that in with *"$sample"*
> * `done` – tells the loop it is finished
> **Don't worry** if loops are confusing at first! They take a litting getting used to, and this not a lesson on loops, just a quick explanation. All exposure helps over time 🙂

This mapping 
