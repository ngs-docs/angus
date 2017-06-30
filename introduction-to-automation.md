# Introduction to automation

## Keeping a log of what you ran

Let's revisit the day's commands for assembly and annotation: assuming
you've run all of the installations already, the commands to download
the data, [run the assembly](genome-assembly.html), and
[run the annotation](prokka_genome_annotation.html), are:

```
# download the data
curl -O -L https://s3.amazonaws.com/public.ged.msu.edu/ecoli_ref-5m.fastq.gz

# run the assembler
~/megahit/megahit --12 ecoli_ref-5m.fastq.gz -o ecoli

# save the assembly
cp ecoli/final.contigs.fa ecoli-assembly.fa

# get some basic metrics on the assembly
~/quast/quast.py ecoli-assembly.fa -o ecoli_report

# run prokka
prokka ecoli-assembly.fa --outdir prokka_annotation --prefix myecoli
```

(here the '#' represent comments that the shell will ignore.)

Let's try running all of that in one go in a new directory:

```
cd ~/
mkdir work2
cd work2
```

and then copy and paste the assembly workflow above. (This will take a few
minutes.)

## While it's running...

A few points here -- 

First, suppose you wanted to edit this workflow, e.g.
change the name you gave to the prokka ``--prefix``?  You could put
this in a Word document, and edit it there, and *then* copy/paste everything
with that one change -- and it should still work.

Second, if you wanted to communicate to a collaborator or a computer
help person what you were doing, you could just send them all the
commands.

Third, this is a pretty good addition to your Methods section, don't
you think?

## After it's done running: put it in a shell script

(Use RStudio to put that into a file in your home directory, and then
create a new directory

## Long-running jobs more generally

(Introduce screen for long-running jobs)

## Scripts: 'bash' shell scripts, R scripts, and Python scripts

(These are all the same things but in different languages.)

## Calling R scripts from bash scripts to produce figures

(You can automate most things.)
