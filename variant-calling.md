# BWA and samtools and variant calling

Here we have two previous lessons:

* [Labs from CTB's 2017 IGG course at UC Davis](https://github.com/ctb/2017-ucdavis-igg201b/blob/master/lab2/README.md) (+ [using bedtools and more samtools](https://github.com/ctb/2017-ucdavis-igg201b/blob/master/lab3/README.md)).  There's also [a small homework exercise](https://github.com/ctb/2017-ucdavis-igg201b/blob/master/hw1/README.md).

* [the 2016 ANGUS tutorial](https://angus.readthedocs.io/en/2016/variant.html)

Proposal:
* merge some of the nice discussion from the 2016 tutorial with CTB's
materials from 2017;
* eliminate the Python code from CTB's materials until maybe later, and just use bedtools;
* use samtools too;
* maybe demo the [Artemis viewer](http://www.sanger.ac.uk/science/tools/artemis) at end. Probably too much to expect people to install it, and the VNC stuff with Jetstream is ...challenging.

-----

## Getting started

[Start up an m1.small instance running Ubuntu 16.04 on Jetstream.](jetstream/boot.html)

## Download data

Goal: get the sequence data!

1. Run:

        mkdir ~/data
        cd ~/data
        curl -O http://dib-training.ucdavis.edu.s3.amazonaws.com/2017-ucdavis-igg201b/SRR2584857.fq.gz

## Map data

Goal: execute a basic mapping

1. Run the following commands to install bwa:

        cd
        curl -L https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2/download > bwa-0.7.15.tar.bz2

        tar xjvf bwa-0.7.15.tar.bz2
        cd bwa-0.7.15
        make

        sudo cp bwa /usr/local/bin
        
        echo 'export PATH=$PATH:/usr/local/bin' >> ~/.bashrc
        source ~/.bashrc

2. Make & change into a working directory:

        mkdir ~/work
        cd ~/work

3. Copy and gunzip the reference:

        wget https://github.com/ctb/2017-ucdavis-igg201b/raw/master/lab1/ecoli-rel606.fa.gz
        gunzip ecoli-rel606.fa.gz
        
4. Prepare it for mapping:

        bwa index ecoli-rel606.fa
        
5. Map! (This will take about 2 minutes.)

        bwa mem -t 6 ecoli-rel606.fa ~/data/SRR2584857.fq.gz > SRR2584857.sam
        
6. Observe!

        head SRR2584857.sam
        
## Visualize mapping

Goal: make it possible to go look at a specific bit of the genome.

1. Install samtools:

        sudo apt-get -y install samtools
        
2. Index the reference genome:

        samtools faidx ecoli-rel606.fa
        
3. Convert the SAM file into a BAM file:

        samtools view -hSbo SRR2584857.bam SRR2584857.sam
        
4. Sort the BAM file by position in genome:

        samtools sort SRR2584857.bam SRR2584857.sorted
        
5. Index the BAM file so that we can randomly access it quickly:

        samtools index SRR2584857.sorted.bam
        
6. Visualize with `tview`:

        samtools tview SRR2584857.sorted.bam ecoli-rel606.fa
        
   `tview` commands of relevance:
   
   * left and right arrows scroll
   * `q` to quit
   * CTRL-h and CTRL-l do "big" scrolls
   * `g ecoli:3931002` will take you to a specific location.
   
## Call variants!

Goal: find places where the reads are systematically different from the
genome.
   
Now we can call variants using
[samtools mpileup](http://samtools.sourceforge.net/mpileup.shtml):

```
samtools mpileup -uD -f ecoli-rel606.fa SRR2584857.sorted.bam | \
    bcftools view -bvcg - > variants.raw.bcf
    
bcftools view variants.raw.bcf > variants.vcf
```

## Look at the VCF file

The output VCF file contains a list of all the variants that samtools
thinks are there. What's in it?

The [official VCF specification](https://samtools.github.io/hts-specs/VCFv4.1.pdf) is a great read...if you're suffering from insomnia.
Let's skip this and just take a quick look at the file.

1. Look at the non-commented lines:

        grep -v ^# variants.vcf

   The first five columns: `CHROM  POS     ID      REF     ALT`.
   It's a little easier to see if you run

        grep -v ^# variants.vcf | less -S

    Use your left and right arrows to scroll, and 'q' to quit.

2. Examine one of the variants with tview:

        samtools tview SRR2584857.sorted.bam ecoli-rel606.fa -p ecoli:920514

   'q' to quit, left arrow to scroll a bit left.
   
Well, at least that variant looks real...

## Look at the VCF file with bedtools.

[bedtools docs](https://bedtools.readthedocs.io/en/latest/)

1. Download and build bedtools:

        cd ~/
        curl -O -L https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
        tar -xzf bedtools-2.26.0.tar.gz

        cd bedtools2
        make
        sudo make install

2. Go back to work:

        cd ~/work
        
3. Download a GFF3 file with annotations for E. coli:

        wget https://github.com/ctb/2017-ucdavis-igg201b/raw/master/lab3/ecoli-rel606.gff.gz .

4. Run bedtools intersect:

        bedtools intersect -a ecoli-rel606.gff.gz -b variants.vcf -wa -u

   [Documentation for bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)

## Extract reads with samtools.

1. Execute:

        samtools view SRR2584857.sorted.bam 'ecoli:920514-920514' > out.bam
        wc -l out.bam

and this will give you the coverage of the relevant position.

## Discussion points / extra things to cover

* What are the drawbacks to mapping-based variant calling? What are
  the positives?

* Where do reference genomes come from?
