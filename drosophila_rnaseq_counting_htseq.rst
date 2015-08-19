===========================================================================
RNA-seq: mapping to a reference genome with tophat and counting with HT-seq
===========================================================================

In this tutorial, we'll use some sample data from a project we did on flies (Drosophila melanogaster) to illustrate how you can use RNA-seq data to look for differentially expressed genes. Here's some brief background on the project: we're trying to understand how different wild-type genetic backgrounds can influence the phenotypic effects of mutations, using the developing fly wing as our model system. We have several mutations that disrupt wing development, and we've backcrossed them into the genetic backgrounds of two different wild-type fly strains (Samarkand/SAM and Oregon-R/ORE). For this tutorial, we've taken a subset of the data--we'll look for expression differences in developing wing tissues between non-mutant and mutant (sd[E3]) flies in each of the two genetic backgrounds, and in flies with a "hybrid" genetic background (i.e., crosses between SAM/ORE flies, again both with and without the mutation). To make things run a little bit faster, we've included only sequence reads that map to X-linked genes.

First, launch an EC2 instance and log in (ami-d05e75b8). Next, create an EBS volume from our snapshot (snap-9b14b2f3). Make sure to create your EC2 instance and your EBS volume in the same availability zone! The snapshot has the raw reads, as well as pre-computed results files so we don't need to wait for every step to finish running before we proceed.

Attach the EBS volume to your instance and mount it in /mnt/ebs/. If you don't know how to do this, ask for a demonstration.

First, we'll need to install a bunch of software. Some of these tools can be installed using apt-get. Note that apt-get does not necessarily always install the most up-to-date versions of this software! You should always double check versions when you do this. For instance, when I was writing this tutorial, apt-get gave me a warning that cufflinks might be out of date, so we're going to install by downloading it directly from the authors.

::

    #Housekeeping stuff for a new instance
    sudo bash
    apt-get update
    apt-get -y install screen git curl gcc make g++ python-dev unzip default-jre pkg-config libncurses5-dev r-base-core r-cran-gplots python-matplotlib sysstat
    mkdir /mnt/ebs
    mount /dev/xvdf /mnt/ebs

    cd /mnt/ebs/
    
    #Install samtools, bowtie, tophat, pysam
    apt-get -y install samtools
    apt-get -y install bowtie
    apt-get -y install tophat
    apt-get -y install python-pip
    pip install pysam
     
    #Download and install HTSeq
    sudo apt-get install build-essential python2.7-dev python-numpy python-matplotlib
    cd /mnt/ebs
    mkdir temp
    cd temp
    curl -O https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1.tar.gz
    tar -xzvf HTSeq-0.6.1.tar.gz
    cd HTSeq-0.6.1/
    python setup.py build
    python setup.py install
     

Next, we need to get our reference genome. This is another area where you want to be careful and pay attention to version numbers--the public data from genome projects are often updated, and gene ID's, coordinates, etc., can sometimes change. At the very least, you need to pay attention to exactly which version you're working with so you can be consistent throughout all your analyses.

In this case, we'll download the reference Drosophila genome and annotation file (which has the ID's and coordinates of known transcripts, etc.) from ensembl. We'll put it in its own directory so we keep our files organized. Then we'll prepare the genomes for use with our software tools by indexing them::

    cd /mnt/ebs/drosophila
    mkdir references
    cd references
    curl -O ftp://ftp.ensembl.org/pub/release-75/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.75.dna.toplevel.fa.gz
    gunzip Drosophila_melanogaster.BDGP5.75.dna.toplevel.fa.gz
    curl -O ftp://ftp.ensembl.org/pub/release-75/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.75.gtf.gz
    gunzip Drosophila_melanogaster.BDGP5.75.gtf.gz
    
    bowtie-build Drosophila_melanogaster.BDGP5.75.dna.toplevel.fa Drosophila_melanogaster.BDGP5.75.dna.toplevel
    samtools faidx Drosophila_melanogaster.BDGP5.75.dna.toplevel.fa


Running this block of code may take a bit of time (especially the bowtie-build command). If you don't want to wait for it, we've included some pre-computed reference files for you. If you want to use our pre-computed references, in the appropriately named precomputed_references directory. If you want to use that instead, simply run the following lines to rename the precomputed_references directory to references::

    cd /mnt/ebs/drosophila
    mv precomputed_references references


Now we're ready for the first step: mapping our RNA-seq reads to the genome. We will use tophat+bowtie1, which together are a splicing-aware read aligner. The raw sequencing reads are in /mnt/ebs/drosophila_reads/. Feel free to take a look at how we've named and organized the files::

    cd /mnt/ebs/drosophila/cleaned_reads
    ls -lh
    
Don't forget that with your reads, you'll probably want to take care of the usual QC steps before you actually begin your mapping. The cleaned_reads directory contains reads that have already been filtered (adapter removal, etc.) and lightly trimmed.

Since we have a lot of files to map, it would take a long time to re-write the mapping commands for each one. And with so many parameters, we might make a mistake or typo. It's usually safer to use a simple shell script with shell variables to be sure that we do the exact same thing to each file. Using well-named shell variables also makes our code a little bit more readable::

    #Create shell variables to store the location of our reference genome and annotation file
    #Note that we are leaving off the .fa from the reference genome file name, because some of the later commands will require just the base of the file name
    reference=/mnt/ebs/drosophila/references/Drosophila_melanogaster.BDGP5.75.dna.toplevel
    annotation=/mnt/ebs/drosophila/references/Drosophila_melanogaster.BDGP5.75.gtf
    
    #Make sure we are in the right directory
    #Let's store all of our mapping results in /mnt/ebs/rnaseq_mapping/ to make sure we stay organized
    cd /mnt/ebs/drosophila
    mkdir rnaseq_mapping
    cd rnaseq_mapping

    #Create an array to hold the names of all our samples
    #Later, we can then cycle through each sample using a simple foor loop
    samples[1]=ORE_wt_rep1
    samples[2]=ORE_wt_rep2
    samples[3]=ORE_sdE3_rep1
    samples[4]=ORE_sdE3_rep2
    samples[5]=SAM_wt_rep1
    samples[6]=SAM_wt_rep2
    samples[7]=SAM_sdE3_rep1
    samples[8]=SAM_sdE3_rep2
    samples[9]=HYB_wt_rep1
    samples[10]=HYB_wt_rep2
    samples[11]=HYB_sdE3_rep1
    samples[12]=HYB_sdE3_rep2
    
    
    #Now we can actually do the mapping
    for i in 1 2 3 4 5 6 7 8 9 10 11 12
    do
        sample=${samples[${i}]}
        #Map the reads
        tophat -p 4 -G ${annotation} -o ${sample} ${reference} /mnt/ebs/drosophila/cleaned_reads/${sample}_1.fastq /mnt/ebs/drosophila/cleaned_reads/${sample}_2.fastq
        #Count the number of reads mapping to each feature using HTSeq
        htseq-count --format=bam --stranded=no --order=pos ${sample}/accepted_hits.bam ${annotation} > ${sample}_htseq_counts.txt
    done


We now have count files for each sample. Take a look at one of the count files using less. You'll notice there are a lot of zeros, but that's partially because we've already filtered the dataset for you to include only reads that map to the X chromosome.

Now we'll need to import them into R to use additional analysis packages to look for differentially expressed genes--in this case, DESeq2. At this point I usually prefer to download these data files and run the analyses locally, because I like using R interactively. I would suggest copying the files using CyberDuck or WinSCP, or through a synchronized Dropbox folder. Once you've got them downloaded, we're now ready to start crunching some numbers.
