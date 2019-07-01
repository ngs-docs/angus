# RNA-seq read quantification with salmon

Learning objectives:

* Install read quantification data
* Learn quantification of RNA-seq data


## Boot up a Jetstream

[Boot an m1.medium Jetstream instance](jetstream/boot.md) and log in.

You should have the quality trimmed data  on your instance in the `~/quality/` 
folder. If not, you can get this data by running:

```
curl -L -o
tar xvf ...
```

## Introduction to Salmon (adapted from salon documentation)

[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) is a tool for fast
transcript quantification from RNA-seq data. It requires a set of target 
transcripts (either from a reference or de-novo assembly) to quantify and 
FASTA/FASTQ file(s) containing your reads. 

Salmon runs in two phases, indexing and quantification. The indexing step is 
independent of the reads, and only need to be run one for a particular set of 
reference transcripts. The quantification step is specific to the set of RNA-seq
reads and is thus run more frequently. 

## Install software

Salmon is installed through conda:

```
conda install -y salmon
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
| --seqBias | learn and correct for sequence-specific biases in the input data |
| --gcBias | learn and correct for fragment-level GC biases in the input data |
| --validateMappings | Enables selective alignment, which improves salmon's sensitivity |

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
less aux_info/meta_info.json
```

We see information about our run parameters and performance. To exit out of the 
file, press `q`.

We can also look at the count file:
```
less -S quant.sf
```

We see our transcript names, as well as the number of reads that aligned to 
each transcript. In our next lesson, we will be reading these quant files into 
R and performing differential expression with them.
"Salmon provides accurate, fast, and bias-aware transcript expression estimates using dual-phase inference" [Patro et al., 2016](http://biorxiv.org/content/early/2016/08/30/021592).

Also see [seqanswers](http://seqanswers.com/) and [biostars](https://www.biostars.org/).

