# Short read quality and trimming

Learning objectives:
- Install software (fastqc, multipc) via conda
- download data
- visualize read quality
- quality filter and trim reads

TODO CTB:

* add instructions on downloading/viewing output HTML with RStudio
* put fastqc/multiqc output files in this repo

Start up a Jetstream m1.small or larger
[as per Jetstream startup instructions](jetstream/boot.html).

---

You should now be logged into your Jetstream computer!  You should see
something like this

```
titus@js-17-71:~$ 
```

## Getting started

Change to your home directory:

```
cd ~/
```

and install FastQC, MultiQC, and trimmomatic:

```
conda install -y fastqc multiqc trimmomatic
```


## Data source

Make a "data" directory:

```
cd ~/
mkdir -p data
cd data
```

and download download some data from the
[Schurch et al, 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/)
yeast RNAseq study:

```
curl -L https://osf.io/5daup/download -o ERR458493.fastq.gz
curl -L https://osf.io/8rvh5/download -o ERR458494.fastq.gz
curl -L https://osf.io/2wvn3/download -o ERR458495.fastq.gz
curl -L https://osf.io/xju4a/download -o ERR458500.fastq.gz
curl -L https://osf.io/nmqe6/download -o ERR458501.fastq.gz
curl -L https://osf.io/qfsze/download -o ERR458502.fastq.gz
```

Let's make sure we downloaded all of our data using md5sum.

```
md5sum *.fastq.gz
```

You should see this:

```
2b8c708cce1fd88e7ddecd51e5ae2154  ERR458493.fastq.gz
36072a519edad4fdc0aeaa67e9afc73b  ERR458494.fastq.gz
7a06e938a99d527f95bafee77c498549  ERR458495.fastq.gz
107aad97e33ef1370cb03e2b4bab9a52  ERR458500.fastq.gz
fe39ff194822b023c488884dbf99a236  ERR458501.fastq.gz
db614de9ed03a035d3d82e5fe2c9c5dc  ERR458502.fastq.gz
```

Now if you type:

```
ls -l
```

you should see something like:

```
-rw-rw-r-- 1 titus titus  59532325 Jun 29 09:22 ERR458493.fastq.gz
-rw-rw-r-- 1 titus titus  58566854 Jun 29 09:22 ERR458494.fastq.gz
-rw-rw-r-- 1 titus titus  58114810 Jun 29 09:22 ERR458495.fastq.gz
-rw-rw-r-- 1 titus titus 102201086 Jun 29 09:22 ERR458500.fastq.gz
-rw-rw-r-- 1 titus titus 101222099 Jun 29 09:22 ERR458501.fastq.gz
-rw-rw-r-- 1 titus titus 100585843 Jun 29 09:22 ERR458502.fastq.gz
```

These are six data files from the yeast study.
of the file.

One problem with these files is that they are writeable - by default,
UNIX makes things writeable by the file owner.  This poses an issue
with creating typos or errors in raw data.  Let's fix that before we
go on any further:

```
chmod a-w *
```

Take a look at their permissions now -- 
```
ls -l
```

and you should see that the 'w' in the original permission string
(`-rw-rw-r--`) has been removed from each file and now it should look like `-r--r--r--`.

We'll talk about what these files are below.

### 1. Copying data into a working location

First, make a working directory; this will be a place where you can futz
around with a copy of the data without messing up your primary data:

```
mkdir -p ~/quality
cd ~/quality
```

Now, make a "virtual copy" of the data in your working directory by
linking it in -- :

```
ln -fs ~/data/* .
```

and you will see that they are in the current directory when you do an
`ls`.

These are FASTQ files -- let's take a look at them:

```
less ERR458493.fastq.gz
```

(use the spacebar to scroll down, and type 'q' to exit `less`)

Question:

* where does the filename come from?
* why are there 1 and 2 in the file names?

Links:

* [FASTQ Format](http://en.wikipedia.org/wiki/FASTQ_format)

### 2. FastQC


We're going to use 
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
summarize the data. We already installed 'fastqc' above, with the
conda command.

Now, run FastQC on two files:

```
fastqc ERR458493.fastq.gz
fastqc ERR458500.fastq.gz
```

Now type 'ls':

```
ls -d *fastqc.zip*
```

to list the files, and you should see:

```
ERR458493_fastqc.zip
ERR458500_fastqc.zip
```

Inside each of the fastqc directories you will find reports from the fastqc program. You can download these files using your RStudio Server console, if you like. (@CTB)


or you can look at these copies of them:

* [ERR458493_fastqc.html](_static/ERR458493_fastqc.html)
* [ERR458500_fastqc.html](_static/ERR458500_fastqc.html)

(UPDATE THESE - put on OSF @CTB)

Questions:

* What should you pay attention to in the FastQC report?
* Which is "better", file 1 or file 2? And why?

Links:

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FastQC tutorial video](http://www.youtube.com/watch?v=bz93ReOv87Y)

There are several caveats about FastQC - the main one is that it only
calculates certain statistics (like duplicated sequences) for subsets
of the data (e.g. duplicate sequences are only analyzed for the first
100,000 sequences in each file


### 3. Trimmomatic

Now we're going to do some trimming!  We'll be using
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), which
(as with fastqc) we've already installed via conda.

The first thing we'll need are the adapters to trim off:

```
cp /opt/miniconda/pkgs/trimmomatic-*/share/trimmomatic-*/adapters/TruSeq2-PE.fa .
```

(you can look at the contents of this file with `cat TruSeq2-PE.fa`)

Now, to run Trimmomatic on both of them:

```
trimmomatic SE ERR458493.fastq.gz \
        ERR458493.qc.fq.gz \
        ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 \
        LEADING:2 TRAILING:2 \
        SLIDINGWINDOW:4:2 \
        MINLEN:25
        
trimmomatic SE ERR458500.fastq.gz \
        ERR458500.qc.fq.gz \
        ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 \
        LEADING:2 TRAILING:2 \
        SLIDINGWINDOW:4:2 \
        MINLEN:25
        
```

You should see output that looks like this:

```
...
Input Reads: 1093957 Surviving: 1092715 (99.89%) Dropped: 1242 (0.11%)
TrimmomaticSE: Completed successfully
```

Questions:

* How do you figure out what the parameters mean?
* How do you figure out what parameters to use?
* What adapters do you use?
* What version of Trimmomatic are we using here? (And FastQC?)
* Do you think parameters are different for RNAseq and genomic data sets?
* What's with these annoyingly long and complicated filenames?
* why are we running R1 and R2 together?

For a discussion of optimal trimming strategies, see 
[MacManes, 2014](http://journal.frontiersin.org/Journal/10.3389/fgene.2014.00013/abstract) -- it's about RNAseq but similar arguments should apply to metagenome
assembly.

Links:

* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

### 4. FastQC again

Run FastQC again on the trimmed files:

```
fastqc ERR458493.qc.fq.gz
fastqc ERR458500.qc.fq.gz
```

And now view my copies of these files: 

* [ERR458493.qc_fastqc.html](_static/ERR458493.qc_fastqc.html)
* [ERR458500.qc_fastqc.html](_static/ERR458500.qc_fastqc.html)

Let's take a look at the output files:

```
less ERR458493.qc.fq.gz                                                       
```

(again, use spacebar to scroll, 'q' to exit less).

### 5. MultiQc
[MultiQC](http://multiqc.info/) aggregates results across many samples into a single report for easy comparison.

Run Mulitqc on both the untrimmed and trimmed files

```
multiqc .
```

And now you should see output that looks like this:

```
[INFO   ]         multiqc : This is MultiQC v1.0
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching '.'
Searching 15 files..  [####################################]  100%
[INFO   ]          fastqc : Found 4 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete
```

You can view the output html file
[multiqc_report.html](_static/multiqc_report.html) by going to RStudio, selecting the file, and saying "view in Web browser."

Questions:

* is the quality trimmed data "better" than before?
* Does it matter that you still have adapters!?
