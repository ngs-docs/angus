# Automating workflows using bash

## Objectives

+ Write a shell script with multiple variables.
+ Incorporate a for loop into a shell script.
+ Understand how to make a workflow more efficient and less error-prone


We have now worked through two workflows, RNA-seq analysis and variant calling.
In both cases, we pasted all of our commands into the terminal, and used for 
loops to run each step of the workflow over our six samples. This is a great way
to understand how each tool works and to establish a workflow you want to use.

Let's imagine now that we decided to download all 672 samples from this dataset
instead of only working with our six samples. We already know how each tool in
our workflow works, and we know how to run them. We can now use bash scripting
to automate our workflow. 

We are going to automate the quality control portion of our workflow. To write 
a script to run our FastQC analysis, we’ll take each of the commands we entered 
to run FastQC and process the output files and put them into a single file with 
a `.sh` extension. The `.sh` is not essential, but serves as a reminder to 
ourselves and to the computer that this is a shell script.

## Starting with automation

We will be constructing our bash script in RStudio. We will construct our bash
script in a text file in the "script" pane of RStudio, and we will execute it 
in the "Terminal" tab of the "console" pane. We could use nano, but RStudio
makes it easier to navigate a text file. 

First, from the terminal, let's make a new directory called `scripts`.


```
cd ~
mkdir scripts
```

## Constructing a bash script

When we write a bash script, we need to add *all* commands that we ran in our
workflow. This includes making and changing directories, moving files, and 
running analysis programs like `fastqc`. 

Let's start by changing into the right directory and running `fastqc`. 

```
cd ~/data/

fastqc *fastq.gz
```

That's it! We now have a bash script that automates `fastqc`! 

Let's run the bash script from our terminal in RStudio

```
bash qc.sh
```

Let's add more to our script. Next we'll organize our output files. 

```
cd ~/data/

fastqc *.fastq.gz

mkdir -p ~/fastqc_untrimmed

mv *.zip ~/fastqc_untrimmed/
mv *.html ~/fastqc_untrimmed/
``` 

We can also add a few pointers that
will print to our screen to let us know what point in the workflow we are in.

```
cd ~/data/

echo "Running FastQC ..."
fastqc *.fastq.gz

mkdir -p ~/fastqc_untrimmed

echo "Saving FastQC results..."
mv *.zip ~/fastqc_untrimmed/
mv *.html ~/fastqc_untrimmed/
```

Let's run this again.

```
bash qc.sh
```

We now see that our echo messages print to the terminal and tell us where we
are in the workflow. If we also run `ls`, we see our `fastqc` output files
are in the new directory we created, `fastqc_untrimmed`.

```
ls ~/fastqc_untrimmed
```

**Challenge** In our script, change directories into `fastqc_untrimmed` and 
then add a command to run `multiqc`. In case you forget, the command to run
`multiqc` in the directory with `fastqc` files is `multiqc .`

With the multiqc step added in, our script should look something like this:

```
cd ~/data/

echo "Running FastQC ..."
fastqc *.fastq.gz 

mkdir -p ~/fastqc_untrimmed

echo "Saving FastQC results..."
mv *.zip ~/fastqc_untrimmed/
mv *.html ~/fastqc_untrimmed/

cd ~/fastqc_untrimmed/

echo "Running multiqc..."
multiqc .
```

Let's add one more line. To the top of the script, add `set -e`. This tells
bash to stop running the script if there is a failure, and it's a good
safe-gaurd to include.

```
set -e
cd ~/data/

echo "Running FastQC ..."
fastqc *.fastq.gz 

mkdir -p ~/fastqc_untrimmed

echo "Saving FastQC results..."
mv *.zip ~/fastqc_untrimmed/
mv *.html ~/fastqc_untrimmed/

cd ~/fastqc_untrimmed/

echo "Running multiqc..."
multiqc .
```

Let's run our script one last time

```
bash qc.sh
```

We should now see the `multiqc` output files in our directory.

```
ls ~/fastqc_untrimmed
```

**Optional Challenge** Transfer the `multiqc` html file to your local computer
using `scp` and open it in your internet browser.
