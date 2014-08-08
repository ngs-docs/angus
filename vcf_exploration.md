Variant calling and exploration of polymorphisms
================================================
Now that we have some experience in R, we will check out a vcf file with polymorphisms from 

## Getting the data and installing extra packages
----------------------------------------------
Installing a bunch of stuff:

get bwa
```  
  cd /root
  wget -O bwa-0.7.10.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download
```

untar and compile (via make) bwa  
 
```  
  tar xvfj bwa-0.7.10.tar.bz2
  cd bwa-0.7.10
  make

  cp bwa /usr/local/bin
```
install some tools
```
  apt-get update
  apt-get -y install samtools screen git curl gcc make g++ python-dev unzip \
        default-jre pkg-config libncurses5-dev r-base-core \
        r-cran-gplots python-matplotlib sysstat libcurl4-openssl-dev libxml2-dev

  git clone https://github.com/schimar/ngs2014_popGen.git

  cd ngs2014_popGen/var_call2/
```

Let's do another round of variant calling
--------------------------------

index the reference genome
```
  bwa index read_file.fq
```
map our reads to the indexed reference genome
```
  bwa aln read_file.fq ref_genome.fna  > mapped_reads.sai
```
Create the SAM file 
```
  bwa samse ref_genome.fna mapped_reads.sai ref_genome.fna > mapped_reads.sam
```
Index the reference genome
```
  samtools faidx ref_genome.fna
```
Convert from SAM to BAM
```
  samtools view -b -S -o mapped_reads.bam mapped_reads.sam
```
Sort the BAM
```
  samtools sort mapped_reads.bam mapped_reads.sorted
```
And index again, but now the sorted BAM file
```
  samtools index mapped_reads.sorted.bam
```
Visualize the alignment
```
  samtools tview SRR098038.sorted.bam REL606.fa
```


Variant exploration with Bioconductor
-------------------------------------

Now simply type R in the shell and:
```r  
  source("http://bioconductor.org/biocLite.R")
  biocLite()
  biocLite("VariantAnnotation")
  biocLite("SNPlocs.Hsapiens.dbSNP.20101109")
  biocLite("BSgenome.Hsapiens.UCSC.hg19_1.3.1000")
```
Quality control
---------------

Now we load the VariantAnnotation package as well as the data. The objective of this exercise is to compare the quality of called SNPs that are located in dbSNP, versus those that are novel.

```r
  library(VariantAnnotation)
  fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
```
Locate the sample data in the file system. Explore the metadata (information about the content of
the file) using scanVcfHeader . Discover the ‘info’ fields VT (variant type), and RSQ (genotype imputation
quality).

```r
  (hdr <- scanVcfHeader(fl))
  info(hdr)[c("VT", "RSQ"),]
```
Input the data and peak at their locations:
```r
  (vcf <- readVcf(fl, "hg19"))

  head(rowData(vcf), 3)
```
SNPs were called with MaCH/thunder (part of GotCloud) , for more info, see :doc: http://genome.sph.umich.edu/wiki/Thunder and http://genome.sph.umich.edu/wiki/MaCH_FAQ. 
Notice that the seqnames (chromosome levels) are set to '22', we want to rename those
```r
  rowData(vcf) <- renameSeqlevels(rowData(vcf), c("22"="ch22"))
```

We now load the SNP database and discover whether our SNPs are in dbSNP
```
  library(SNPlocs.Hsapiens.dbSNP.20101109)

  destination <- tempfile()
  pre <- FilterRules(list(isLowCoverageExomeSnp = function(x) {
  grepl("LOWCOV,EXOME", x, fixed=TRUE)
  }))
  filt <- FilterRules(list(isSNP = function(x) info(x)$VT == "SNP"))
  snpFilt <- filterVcf(fl, "hg19", destination, prefilters=pre, filters= filt)
  vcf_filt <- readVcf(snpFilt, "hg19")

  rowData(vcf) 
  rowData(vcf_filt)
```
If we compare vcf and vcf_filt, we see that of the 10376 SNPs in our initial vcf file, 794 are in the database. 

```
  inDbSNP <- rownames(vcf) %in% rownames(vcf_filt)
  table(inDbSNP)
  metrics <- data.frame(inDbSNP = inDbSNP, RSQ = info(vcf)$RSQ)
```
Let's finally visualize it:
```
  library(ggplot2)
  ggplot(metrics, aes(RSQ, fill=inDbSNP)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(name="MaCH / Thunder Imputation Quality") +
  scale_y_continuous(name="Density") +
  theme(legend.position="top")
```
(This won't work in R on EC2, simply because we can't run X11 through an ssh connection)




