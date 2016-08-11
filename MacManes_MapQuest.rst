================================================
Map Quest
================================================

I am going to give you some data and a genome (well, 500Mb of the genome), I want you to quality trim and map using the Trimming lecture and the Mapping Lecture.

** Launch an AMI. For this exercise, try a **c4.4xlarge** with a 200Gb EBS volume.


You will need to install `skewer`, `bwa`, and `samtools`. Maybe try `brew`???

**Download the genome and the reads**

::

    https://s3.amazonaws.com/macmanes_share/Cavia.fa
    curl -LO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR636/SRR636898/SRR636898_1.fastq.gz
    curl -LO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR636/SRR636898/SRR636898_2.fastq.gz


**Use Skewer to trim reads**

1. For people whose 1st name starts with between A-F, skip this step to allow us to evaluate the performance of mapping untrimmed reads.
2. Everybody else, choose a random number between 1 and 30 and trim at that level..

**Map the reads with BWA** You will need to figure out how to map paired end reads. May luck be ever in your favor.. Also, there is google and the help info.

Type `bwa mem` to get a list of the options. You will need to index the genome using `bwa index` before you map.


::

    bwa mem [you need to figure out the command] \
    | samtools view -@10 -Sb - \
    | samtools sort -T sort -O bam -@10 -o sorted.bam -


**Evaluate your mapping data**

::

    samtools flagstat sorted.bam

**Add your data** about mappability

https://docs.google.com/spreadsheets/d/1svCJMefIy-BTBSWl5Vixf0UPgLIGOrp5u6YRoyMbSng/edit?usp=sharing


**Super bonus points to anybody that maps with a differnt mapper** e.g., Bowtie, STAR, VelociMapper, HiSat2, etc..

**TERMINATE YOUR INSTANCE!!!**
