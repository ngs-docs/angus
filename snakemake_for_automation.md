# Snakemake for Automation

## Objectives

+ Identify cases where workflow managers are helpful for automation
+ Understand the components of a Snakefile: rules, inputs, outputs, and actions.
+ Write and run a Snakefile.

## Getting started

Start up a Jetstream m1.medium or larger
[as per Jetstream startup instructions](jetstream/boot.html).

---

You should now be logged into your Jetstream computer!  You should see
something like this

```
diblynn@js-17-71:~$
```

Make sure you're in your home directory:

```
cd ~
```

Install `snakemake` using conda. 

```
conda install -y -c bioconda snakemake-minimal
```

## Automation

In both our RNA-seq workflow and our mapping and variant calling workflow, we 
performed many steps. We performed steps like quality control and analysis using
`fastqc` and `trimmomatic`. We performed these steps on 6 files using for loops.

In our last lesson, we automated these steps using a bash script. We put all of 
our commands into one file so that we only had to run one command to orchestrate
our quality control workflow. Bash scripting for automation is really powerful!

Let's revisit our bash script for running and organizing our fastqc results:

```
set -e
cd ~/data/

echo "Running FastQC ..."
fastqc *.fastq.gz

mkdir -p ~/fastqc_untrimmed

echo "Saving FastQC results..."
mv *.zip ~/fastqc_untrimmed/
mv *.html ~/fastqc_untrimmed/

cd ~/fastqc_untrimmed_reads/

echo "Unzipping..."
for filename in *.zip
do
    unzip $filename
done

echo "Saving summary..."
cat */summary.txt > ~/fastqc_summaries.txt
```

We can run it using this command:

```
bash read_qc.sh
```

Oh crap! We realize just after we've finished `Running FastQC` that we wanted 
to move our summary file to a subdirectory! Quick, press `ctrl - C` and cancel
the run!

Even if we made this change though, we're in a bit of a pickle. We want to
re-start our script that automates the runs and the file movements for us, but
we already ran the first part of our file! Let's add comment characters to the
lines we know already ran and then re-run the script:

```
set -e
cd ~/data/

#echo "Running FastQC ..."
#fastqc *.fastq.gz

#mkdir -p ~/fastqc_untrimmed

echo "Saving FastQC results..."
mv *.zip ~/fastqc_untrimmed/
mv *.html ~/fastqc_untrimmed/

cd ~/fastqc_untrimmed_reads/

echo "Unzipping..."
for filename in *.zip
do
    unzip $filename
done

echo "Saving summary..."
cat */summary.txt > ~/fastqc_summaries.txt
```

Now we can re-run the script:

```
bash read_qc.sh 
```

This (maybe) worked out ok this time. However, it's hard to be sure we know 
where we stopped our command. For this reason and many others, we use workflow
managers to keep track of the things that have and haven't run yet! 

We'll be using snakemake for automation. 

## Introduction to Snakemake

What is Snakemake and why are we using it?

The Snakemake workflow management system is a tool to create reproducible and 
scalable data analyses. It orchestrates and keeps track of all the different
steps of workflows that have been run so you don't have to! It has a lot of 
wonderful features that can be invoked for different applications, making it
very flexible while maintaining human interpretability.  

There are many different tools that researchers use to automate computational
workflows. We selected snakemake for the following reasons:

+ It was written by a bioinformatician for bioinformatics workflows.
+ It’s free, open-source, and conda-installable
+ Snakemake works cross-platform (Windows, MacOS, Linux) and is compatible with 
all HPC schedulers. It works on laptops, the cloud, and clusters without
modification to the main workflow (as long as you have enough compute 
resources!).
+ Snakemake is written using Python, but supports bash and R code as well.
+ Anything that you can do in Python, you can do with Snakemake (since you can 
pretty much execute arbitrary Python code anywhere).

Our goal is to automate our example workflow using snakemake! 

## Starting with Snakemake



rule all:
    input:
        expand("outputs/fastqci/{samples}.html", sample = SAMPLES)
    
SAMPLES=["ERR..", "ERR.."]

rule fastqc_raw:
    input:
    
