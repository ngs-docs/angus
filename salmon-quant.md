# RNA-seq read quantification with salmon

Learning objectives:

* Install read quantification program Salmon
* Learn quantification of RNA-seq data


## Boot up a Jetstream

[Boot an m1.medium Jetstream instance](jetstream/boot.md) and log in.

You should have the quality trimmed data  on your instance in the `~/quality/` 
folder. If not, you can get this data by running:

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

## Introduction to Salmon (adapted from salmon documentation)

[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) is a tool for fast
transcript quantification from RNA-seq data. It requires a set of target 
transcripts (either from a reference or de-novo assembly) to quantify and 
FASTA/FASTQ file(s) containing your reads. 

Salmon runs in two phases, indexing and quantification. The indexing step is 
independent of the reads, and only needs to be run one for a particular set of 
reference transcripts. The quantification step is specific to the set of RNA-seq
reads and is thus run more frequently. 

## Install software

Salmon is installed through conda.

Let's make sure our conda channels are loaded.

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

```
conda install -y -c bioconda salmon
```

## Make a new working directory and link the quality trimmed data

We will be using the same data as before 
([Schurch et al, 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/)),
so the following commands will create a new folder `quant` and link the data in:

```
cd ~
mkdir -p quant
cd quant

ln -fs ~/quality/*qc.fq.gz .
ls
```

If you don't have the data from the previous lesson, you can download it
## Download the yeast transcriptome:

```
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/146/045/GCA_000146045.2_R64/GCA_000146045.2_R64_rna_from_genomic.fna.gz
```

## Index the yeast transcriptome:

```
salmon index --index sc_index --type quasi --transcripts GCA_000146045.2_R64_rna_from_genomic.fna.gz
```

## Run salmon on all the samples:

```
for i in *.qc.fq.gz
do
   salmon quant -i sc_index --libType A -r ${i} -o ${i}_quant --seqBias --gcBias --validateMappings
done
```

What do all of these flags do?

| Flag | Meaning |
|--------------------|------------------------------------------------------------------|
| -i | path to index folder |
| --libType | The library type of the reads you are quantifying. `A` allows salmon to automatically detect the library type of the reads you are quantifying. |
| -r | Input file (for single-end reads) |
| -o | output folder |
| --seqBias | learn and correct for [sequence-specific biases](https://salmon.readthedocs.io/en/latest/salmon.html#seqbias) in the input data |
| --gcBias | learn and correct for [fragment-level GC biases](https://salmon.readthedocs.io/en/latest/salmon.html#seqbias) in the input data |
| --validateMappings | Enables [selective alignment](https://salmon.readthedocs.io/en/latest/salmon.html#validatemappings), which improves salmon's sensitivity |

As Salmon is running, a lot of information is printed to the screen. For example,
we can see the mapping rate as the sample finishes:

```
[2019-06-29 18:39:18.367] [jointLog] [info] Mapping rate = 86.2687%
```

When it finished, Salmon outputs a folder for each input RNA-seq sample. Let's
take a look at one of the output folders.

```
ls ERR458493.qc.fq.gz_quant
```

You should see output like this:

```
aux_info/               lib_format_counts.json
cmd_info.json           logs/
libParams/              quant.sf
```

The information we saw scroll through our screens is captured in a log file in
`aux_info/`. 

```
less ERR458493.qc.fq.gz_quant/aux_info/meta_info.json
```

We see information about our run parameters and performance. To exit out of the 
file, press `q`.

We can also look at the count file:
```
less -S ERR458493.qc.fq.gz_quant/quant.sf
```


We see our transcript names, as well as the number of reads that aligned to 
each transcript. 

<blockquote>
<center><b>PRACTICE!</b></center>

How can we use the `grep` command to find the percent of reads mapped in ALL our json files? Make sure to print out the name of each file before you print the percent of reads mapped!

<div class="toggle-header closed">
    <strong>Solution</strong>
</div>

<div class="toggle-content docutils container" style="width:100%">

<div class="highlight-bash notranslate">
<div class="highlight">
<pre>
<span class="nb">cd ~/quant/
for infile in *_quant/aux_info/meta_info.json
do
echo ${infile}
grep "percent_mapped" ${infile}
done</span>
</pre>
</div>
</div>

First we use `echo` to print the name of the file. Then, we use `grep` to find and print the line containing the percent of reads mapped.
<br>

</div>
</blockquote>

In our next lesson, we will be reading these quant files into 
R and performing differential expression with them.
"Salmon provides accurate, fast, and bias-aware transcript expression estimates using dual-phase inference" [Patro et al., 2016](http://biorxiv.org/content/early/2016/08/30/021592).

Also see [seqanswers](http://seqanswers.com/) and [biostars](https://www.biostars.org/).

