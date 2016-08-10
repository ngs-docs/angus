================================================
Map Quest
================================================

I am going to give you some data and a genome, I want you to quality trim and map using the Trimming lecture and the Mapping Lecture.

** Launch and AMI. For this exercise, we will use a **c4.4xlarge** with a 200Gb EBS volume.


You will need to install `skewer`, `bwa`, and `samtools`.

**Download the genome and the reads**

::

    curl -LO ftp://ftp.ensembl.org/pub/release-85/fasta/cavia_porcellus/dna/Cavia_porcellus.cavPor3.dna.nonchromosomal.fa.gz
    curl -LO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR636/SRR636898/SRR636898_1.fastq.gz
    curl -LO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR636/SRR636898/SRR636898_2.fastq.gz


**Use Skewer to trim reads** - for people whose name ends with between A-F, skip this step to allow us to evaluate the performance of mapping untrimmed reads.

**Map the reads with BWA** You will need to figure out how to map paired end reads. May luck be ever in your favor.. Also, there is google and the help info.

**Evaluate your mapping data**
