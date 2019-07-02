# Snakemake for Automation

## Objectives
+
+
+

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
conda install -y -c bioconda snakemake
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
we already ran the first part of our file! Let's add a comment 
However, while we can use bash scripting for automation, there are also other
software tools that are written specifically for this task. We'll be using
snakemake for automation. 

 

For example, in our mapping and variant calling lesson, we
used quality trimming, 

