# Short read quality and trimming

Start up a Jetstream m1.small or larger
[as per Jetstream startup instructions](https://angus.readthedocs.io/en/2017/jetstream/boot.html).

---

You should now be logged into your Jetstream computer!  You should see
something like this

```
titus@js-17-71:~$ 
```

## Installing some software

Run:

```
sudo apt-get -y update && \
sudo apt-get -y install trimmomatic fastqc python-pip \
   samtools zlib1g-dev ncurses-dev python-dev
pip install multiqc   
```
`apt-get install` doesn't work properly for `fastqc`. So we will update the default `fastqc` version using the following commands

```
cd ~/
wget https://launchpad.net/ubuntu/+archive/primary/+files/fastqc_0.11.5+dfsg-3_all.deb && \
sudo dpkg -i fastqc_0.11.5+dfsg-3_all.deb && \
sudo apt-get install -f
```

## Data source

We're going to be using a subset of metagenomic data from [Hu et al.,
2016](http://mbio.asm.org/content/7/1/e01669-15.full). This paper
from the Banfield lab samples some relatively low diversity environments
and finds a bunch of nearly complete genomes.

### 1. Copying in some data to work with.

We've loaded subsets of the data onto a public location for you, to
make everything faster for today's work.  We're going to put the
files on your computer locally under the directory ``~/data``:

```
mkdir ~/data
```

Next, let's grab part of the data set:

```
cd ~/data
curl -O -L https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/metagenomics-scripps-2016-10-12/SRR1976948_1.fastq.gz
curl -O -L https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/metagenomics-scripps-2016-10-12/SRR1976948_2.fastq.gz
```
 
Let's make sure we downloaded all of our data using md5sum.

```
md5sum SRR1976948_1.fastq.gz SRR1976948_2.fastq.gz
```

You should see this:

```
37bc70919a21fccb134ff2fefcda03ce  SRR1976948_1.fastq.gz
29919864e4650e633cc409688b9748e2  SRR1976948_2.fastq.gz
```

Now if you type:

```
ls -l
```

you should see something like:

```
total 346936
-rw-rw-r-- 1 ubuntu ubuntu 169620631 Oct 11 23:37 SRR1976948_1.fastq.gz
-rw-rw-r-- 1 ubuntu ubuntu 185636992 Oct 11 23:38 SRR1976948_2.fastq.gz
```

These are 1m read subsets of the original data, taken from the beginning
of the file.

One problem with these files is that they are writeable - by default, UNIX
makes things writeable by the file owner.  This poses an issue with creating typos or errors in raw data.  Let's fix that before we go
on any further:

```
chmod u-w *
```

We'll talk about what these files are below.

### 1. Copying data into a working location

First, make a working directory; this will be a place where you can futz
around with a copy of the data without messing up your primary data:

```
mkdir ~/work
cd ~/work
```

Now, make a "virtual copy" of the data in your working directory by
linking it in -- :

```
ln -fs ~/data/* .
```

These are FASTQ files -- let's take a look at them:

```
less SRR1976948_1.fastq.gz
```

(use the spacebar to scroll down, and type 'q' to exit 'less')

Question:

* where does the filename come from?
* why are there 1 and 2 in the file names?

Links:

* [FASTQ Format](http://en.wikipedia.org/wiki/FASTQ_format)

### 2. FastQC


We're going to use 
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
summarize the data. We already installed 'fastqc' on our computer for
you.

Now, run FastQC on two files:

```
fastqc SRR1976948_1.fastq.gz
fastqc SRR1976948_2.fastq.gz
```

Now type 'ls':

```
ls -d *fastqc.zip*
```

to list the files, and you should see:

```
SRR1976948_1_fastqc.zip
SRR1976948_2_fastqc.zip
```

Inside each of the fatqc directories you will find reports from the fastqc. You can download these files using your RStudio Server console, if you like. To install and run an RStudio Server, go [here](https://angus.readthedocs.io/en/2017/visualizing-blast-scores-with-RStudio.html#installing-and-running-rstudio-on-jetstream). 


or you can look at these copies of them:

* [SRR1976948_1_fastqc/fastqc_report.html](http://2017-ucsc-metagenomics.readthedocs.io/en/latest/_static/SRR1976948_1_fastqc/fastqc_report.html)
* [SRR1976948_2_fastqc/fastqc_report.html](http://2017-ucsc-metagenomics.readthedocs.io/en/latest/_static/SRR1976948_2_fastqc/fastqc_report.html)

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
(as with fastqc) we've already installed via apt-get.

The first thing we'll need are the adapters to trim off:

```
curl -O -L http://dib-training.ucdavis.edu.s3.amazonaws.com/mRNAseq-semi-2015-03-04/TruSeq2-PE.fa
```

Now, to run Trimmomatic:

```
TrimmomaticPE SRR1976948_1.fastq.gz \
                 SRR1976948_2.fastq.gz \
        SRR1976948_1.qc.fq.gz s1_se \
        SRR1976948_2.qc.fq.gz s2_se \
        ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 \
        LEADING:2 TRAILING:2 \
        SLIDINGWINDOW:4:2 \
        MINLEN:25
```

You should see output that looks like this:

```
...
Input Read Pairs: 1000000 Both Surviving: 885734 (88.57%) Forward Only Surviving: 114262 (11.43%) Reverse Only Surviving: 4 (0.00%) Dropped: 0 (0.00%)
TrimmomaticPE: Completed successfully
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
[MacManes, 2014](http://journal.frontiersin.org/Journal/10.3389/fgene.2014.00013/abstract)
-- it's about RNAseq but similar arguments should apply to metagenome
assembly.

Links:

* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

### 4. FastQC again

Run FastQC again on the trimmed files:

```
fastqc SRR1976948_1.qc.fq.gz
fastqc SRR1976948_2.qc.fq.gz
```

And now view my copies of these files: 

* [SRR1976948_1.qc_fastqc/fastqc_report.html](http://2016-metagenomics-sio.readthedocs.io/en/work/_static/SRR1976948_1.qc_fastqc/fastqc_report.html)
* [SRR1976948_2.qc_fastqc/fastqc_report.html](http://2016-metagenomics-sio.readthedocs.io/en/work/_static/SRR1976948_2.qc_fastqc/fastqc_report.html)

Let's take a look at the output files:

```
less SRR1976948_1.qc.fq.gz
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

You can view output html file
[multiqc_report.html](_static/multiqc_report.html)

Questions:

* is the quality trimmed data "better" than before?
* Does it matter that you still have adapters!?
