# Genome assembly - some basics

In [variant calling](https://angus.readthedocs.io/en/2017/variant-calling.html), we mapped reads to a reference and looked systematically for differences.

But what if you don't have a reference? How do you construct one?

The answer is *de novo* assembly, and the basic idea is you feed in
your reads and you get out a bunch of *contigs*, that represent
stretches of DNA present in the reads that don't have any long repeats
or much significant polymorphism.  Like everything else, the basic idea
is that you run a program, feed in the reads, and get out a pile of
assembled DNA.

MEGAHIT, used below, works well for assembly short-read data sets from
genomes and metagenomes.  For transcriptomes, you might use Trinity -
see
[the eel-pond protocol](https://eel-pond.readthedocs.io/en/latest/)
for our guide to doing RNA-seq assembly.

### Start up a Jetstream instance

Goal: provide a platform to run stuff on.

[Start up an m1.medium instance running Ubuntu 16.04 on Jetstream.](jetstream/boot.html)

### Install the MEGAHIT assembler

Check out and build [MEGAHIT](https://www.ncbi.nlm.nih.gov/pubmed/27012178):

    git clone https://github.com/voutcn/megahit.git
    cd megahit
    make -j 6

### Download an E. coli data set

Grab the following E. coli data set:

    mkdir ~/work
    cd ~/work
    
    curl -O -L https://s3.amazonaws.com/public.ged.msu.edu/ecoli_ref-5m.fastq.gz
    
### Run the assembler

Assemble the E. coli data set with MEGAHIT:

    ~/megahit/megahit --12 ecoli_ref-5m.fastq.gz -o ecoli

(This will take about 3 minutes.)  You should see something like:

    --- [STAT] 117 contigs, total 4577284 bp, min 220 bp, max 246618 bp, avg 39122 bp, N50 105708 bp
    --- [Fri Feb 10 14:33:59 2017] ALL DONE. Time elapsed: 342.060158 seconds ---

at the end.

Questions while we're waiting:

* how many reads are there?

* how long are they?

* are they paired end or single-ended?

* are they trimmed?

...and how would we find out?

Also, what expectation do we have for this genome in terms of size,
content, etc?

### Looking at the assembly

First, save the assembly:

    cp ecoli/final.contigs.fa ecoli-assembly.fa
    
Now, look at the beginning:

    head ecoli-assembly.fa
    
It's DNA! Yay!

So this is the curse and the benefit of assembly - you go through
some amount of work to get your data, QC it, clean it up, and assemble it,
but then you're faced with a pile of assembled but unannotated results!
(We'll talk about annotation next tutorial.)

But before you put effort into annotating the assembly, you should think
about whether it's any good...

### Measuring the assembly

Install [QUAST](http://quast.sourceforge.net/quast):

    cd ~/
    git clone https://github.com/ablab/quast.git -b release_4.5
    export PYTHONPATH=$(pwd)/quast/libs/

Run QUAST on your assembly:

    cd ~/work
    ~/quast/quast.py ecoli-assembly.fa -o ecoli_report
    
and now take a look at the report:

    cat ecoli_report/report.txt
    
You should see something like:

```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                    ecoli-assembly
# contigs (>= 0 bp)         117           
# contigs (>= 1000 bp)      91            
# contigs (>= 5000 bp)      69            
# contigs (>= 10000 bp)     64            
# contigs (>= 25000 bp)     52            
# contigs (>= 50000 bp)     32            
Total length (>= 0 bp)      4577548       
Total length (>= 1000 bp)   4565216       
Total length (>= 5000 bp)   4508381       
Total length (>= 10000 bp)  4471170       
Total length (>= 25000 bp)  4296203       
Total length (>= 50000 bp)  3578898       
# contigs                   101           
Largest contig              246618        
Total length                4572094       
GC (%)                      50.75         
N50                         105709        
N75                         53842         
L50                         15            
L75                         30            
# N's per 100 kbp           0.00  
```

This is a set of summary stats about your assembly. Are they good?
Bad? How would you know?

## What are other metrics you could use to evaluate your assembly?

This is a good opportunity for brainstorming and group thinking :)

## End of day

Question: why so many contigs?!

And what do you do with a bunch of assembled DNA sequence anyway?
