================================================
Assembly Quest
================================================

I am going to give you some data and a genome (well, 500Mb of the genome), I want you to quality trim and map using the Trimming lecture and the Mapping Lecture.

** Launch an AMI. For this exercise, try a **c4.8xlarge** with a 200Gb EBS volume. Wow yes this is a BIG machine!!


You will need to install `abyss`

**Download the reads** They have already been trimmed for you :)

::

    Paired-end data

    curl -LO https://s3.amazonaws.com/macmanes_share/ERR1328557_1.fastq.gz
    curl -LO https://s3.amazonaws.com/macmanes_share/ERR1328557_2.fastq.gz

    Mate-pair data

    curl -LO https://s3.amazonaws.com/macmanes_share/ERR1328558_1.fastq.gz
    curl -LO https://s3.amazonaws.com/macmanes_share/ERR1328558_2.fastq.gz


**Use ABySS to assemble reads**

1. For people whose 1st name starts with between A-F, skip this step to allow us to evaluate the performance of mapping untrimmed reads.
2. Everybody else, choose a random number between 1 and 30 and trim at that level..


**Evaluate your assembly data**

::

    samtools flagstat sorted.bam

**Add your data** about mappability

https://docs.google.com/spreadsheets/d/1svCJMefIy-BTBSWl5Vixf0UPgLIGOrp5u6YRoyMbSng/edit?usp=sharing


**Super bonus points to anybody that maps with a different assembler** e.g., SPAdes, AllPaths, Velvet, ...

**TERMINATE YOUR INSTANCE!!!**
