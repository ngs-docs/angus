---
title: "Association Mapping: GWAS and Sequencing Data"
author: "Tiffany Timbers"
date: "August 20, 2015"
output: html_document
---

## Introduction to GWAS
A genome-wide association study (GWAS) is a method which tests whether genetic variants are associated with a trait (*e.g.,* phenotype or
disease outcome). The method compares the DNA sequences of two groups: individuals with the phenotype or disease outcome (cases)
and individuals without the phenotype or disease outcome (controls). Individuals DNA is read via SNP arrays or sequencing and
if one variant appears significantly more frequently in the pool of individuals with the disease then the variant is said to be 
"associated" with the disease. 

## New considerations when performing GWAS on sequencing data
GWAS was developed for data from SNP microarrays and has begun to be subsequently adapted for next-generation sequencing data. Thus there
are several considerations when doing GWAS with sequencing data, including:
* how to deal with rare variants (these didn't exist in SNP array experiments)?
* focus only on the exome? or look at whole-genome?
* how to deal with synonomous mutations?

## Advantages when performing GWAS on sequencing data
* potentially increased power due to rare variants potentially being more sever than common variants 
* sequencing gives a complete picture of genetic variation (*e.g.,* SNPs, Copy Number Variants (CNVs) and insertions-deletions)
* can predict severity of variants from sequencing data and use this to weight analysis

## Outline of project workflow for GWAS of sequencing data using the Sequence Kernel Association Test (SKAT)
1. Create a single `.vcf` file containing all the variants (as rows) and all the individuals (as columns)
2. Choose to perform association analysis via either gene sets or windows (*e.g.* SNP/variant set is defined to be those 
variants which fall within a chosen region size).
3. Create a SNP set ID (`.SSID`) file which indicates which variants belong to which set (*e.g.,* gene)
4. (optional) Create a custom weights file to assign weights to each variant
5. Create binary `plink` files from the .vcf file
6. Append the phenotype/disease outcome to the `plink` `.fam` file
7. (optional) Create co-variate file
8. Apply association test (*e.g.,* `SKAT`)
9. Apply multiple testing correction to determine which genes/regions are significantly associated with the
phenotype/disease outcome

## Let's do it!

###Dependencies
Our program dependencies are: the `Bash Shell`, `Perl`, `bgzip`, `tabix`, `R` and `R` packages `SKAT`, `dplyr`, `stringr` and `fdrtool`, as well as the program `plink`. Our data 
dependencies are: one `.vcf` file for each individual and a single phenotype/disease outcome file
with 2 columns (the first being the individual ID and the second being the phenotype/disease outcome).

### Create a merged `.vcf` file containing all the variants and all the individuals

In the `data` directory we have a collection of 480 `.vcf` files, each of which represent the whole-genome sequences of individual 
isogenic nematode, *Caenorhabditis elegans*, strains from the [Million Mutation Project](http://genome.sfu.ca/mmp/about.html) (Thompson *et. al., Genome Research*, 2013) called 
against the VC2010 reference strain (a derivative of N2). Given that these are isogenic strains, we will treat each strain as an 
individual for this analysis. 
```
# Bash Shell
ls data/*.vcf | head -5
```

Next, we use the VCFtools program `vcf-merge` to create a single merged `.vcf` file containing all the variants and all 
the individuals.
```
# Bash Shell
## CAUTION: LONG-RUN TIME PROCESS

# pre-process .vcf files so they are compatible with vcf-merge
for file in data/*.vcf; do bgzip -c $file > $file.gz; tabix -p vcf $file.gz; done
ls data/*.vcf.gz | head -5

# merge them together (note - if you run this twice you will get an error message because it will try to merge MMP_vcf-merge.vcf.gz)
vcf-merge data/*.vcf.gz | bgzip -c > merged
mv merged data/MMP_vcf-merge.vcf.gz
ls data/MMP_vcf-merge.vcf.gz
```

Then we uncompress `MMP.vcf.gz` so we can use it to create a snp set ID (`SSID`) file, as well as binary `plink` files.
```
# Bash Shell

gunzip data/MMP_vcf-merge.vcf.gz
grep -v -E '##' data/MMP_vcf-merge.vcf | cut -f 10-15 | head -5
```

When we merged the files, we see that for each individual we either have "1/1" or "." for a genotype. vcf-merge assumed
that if there was no variant info for an individual when it merged the file that the data was missing. In
the case of this data this is not true. For this dataset, each individual .vcf files simply did not contain all variants for this whole dataset
because each file would be ~ 860,000 lines long and filled with mostly "0/0". But we need those "0/0" for our analysis and so we need to replace 
the "." with "0/0". We can use the gsub function in awk to do this. 
```
# Bash Shell

# replace  "." with "0/0""
awk '{gsub("\t[.]","\t0/0",$0); print;}' data/MMP_vcf-merge.vcf > data/MMP.vcf
grep -v -E '##' data/MMP.vcf | cut -f 10-15 | head -5
```

### Create a SNP set ID (`.SSID`) file which indicates which variants belong to which set (*e.g.,* gene)

For this analysis, we will define our SNP/variant sets by genes, thus we will omit all other DNA outside of exons and 
introns from our analysis. We also want to omit silent mutations in coding regions To do this we need to create a snp 
set ID (`.SSID`) file which indicates which variants belong to which set (*e.g.,* gene). `SSID` files contain 2 columns, 
the first being the gene name and the second being the variant name. Note - you can define these SNP sets however you 
want for your own analysis. For example, you could define them by biochemical pathways as opposed to genes to potentially 
increase power.

First we extract all the lines from the `.vcf` file with gene names using a series of `grep` commands. Note - these 
commands depend on how you defined your gene/sequence names in the`INFO` column. In the case of the `.vcf` files
we are using here, gene/sequence names are found either after `SN=` or `CODING=`.
```
# Bash Shell

grep -E -v '#' data/MMP.vcf | grep -E 'SN=[A-Z0-9.]{4,}|CODING=' | grep -E -v 'SN=I;' | grep -E -v 'GRANT=0|GRANT=NA' | cut -f 1-8 > data/gene_variants.txt
head -5 data/gene_variants.txt
```

Next we need to extract the gene and the sequence names from `data/gene_variants.txt` and save them in a file called `MMP.SSID`.
We will do this in R for simplicity, but it could be done faster with another language (*e.g.,* Perl)
```
# R

## load data/gene_variants.txt into R
gene_variants <- read.table(file = "data/gene_variants.txt", stringsAsFactors = FALSE)
head(gene_variants)

# load stringr library to perform regular expressions
library(stringr)

# use str_extract() to capture gene/sequence name from INFO column (V8 in the gene_variants data frame)
pattern  <- "SN=[a-zA-Z0-9.]{1,}|CODING=[a-zA-Z0-9.]{1,}"
gene  <- str_extract(gene_variants$V8, pattern)
gene  <- sub("SN=|CODING=", "", gene)
head(gene)

# make a new data frame containing 2 columns, gene and variants (gene_variants$V3)
SSID <- data.frame(gene, gene_variants$V3)
head(SSID)

# remove NA's (sometimes CODING=NA and we want to remove those)
SSID <- subset(SSID, gene!="NA")

# save the SSID file
write.table(x = SSID, file = "data/MMP.SSID", row.names = FALSE, col.names = FALSE, append = FALSE, quote = FALSE)
```

### Create binary `plink` files from the .vcf file

`plink` (Purcell *et al., American Journal of Human Genetics*, 2007) is one of the most widely used programs to perform traditional GWAS (*e.g.,* using SNP array data) 
and thus most newly developed tools accept `plink` files as input. In addition to association analysis tools, `plink` now
also contains tools to convert `.vcf` files to `plink` files.
```
# Bash Shell

./plink --vcf data/MMP.vcf --allow-extra-chr --no-fid --no-parents --no-sex --no-pheno --out data/MMP
ls data/MMP*
```

### Append the phenotype/disease outcome to the `plink` `.fam` file

The `.fam` file does not yet contain the phenotype/disease outcome data. We need to manually add this to the `.fam` by removing
column 6 from the `.fam` file and joining it by strain/individual with the file containing our phenotype/disease outcome data.
Because these files are smaller (2-6 column and hundreds of to tens of thousands of rows) we can load these files into memory and
will take advantage of the powerful `dplyr` package in `R` to do this.
```
# R

# load dplyr library (we will use this to join the files)
library(dplyr)

# read in phenotypes file
pheno_file <- tbl_df(read.table(file = "data/phenotype_dichotomous.csv", header = TRUE, stringsAsFactors = FALSE))
head(pheno_file)

# read in fam file
fam_file <- tbl_df(read.table(file = "data/MMP.fam", stringsAsFactors = FALSE))
head(fam_file)

# name columns and drop column V6 (we will replace this with the phenotype column)
colnames(fam_file)[1] <- "strain"
fam_file$V6 <- NULL
head(fam_file)

#  use dplyr to join columns 1-5 of fam_file with columns 1-2 of pheno_file, join by column 1
fam_w_phenos <- left_join(x = fam_file, y = pheno_file)
head(fam_w_phenos)
write.table(x = fam_w_phenos, file = "data/MMP.fam", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE)
```

```
# Bash Shell

head -5 data/MMP.fam
```

### Apply association test (*e.g.,* `SKAT`)

Now that we have our `SSID` file and our binary `plink` files we can do the association analysis. We will work in R
to do this and use the SKAT package (Wu *et al., American Journal of Human Genetics*, 2011).
```
# R
## CAUTION: LONG-RUN TIME PROCESS

# load SKAT library (to do association analysis)
library(SKAT)

# Generate a SNP set data file (`.SSD`) and an `.info` file from binary plink formated data files using user specified SNP sets.
Generate_SSD_SetID('data/MMP.bed', 'data/MMP.bim', 'data/MMP.fam', 'data/MMP.SSID', 'data/MMP.SSD', 'MMP.info')

# read in fam file
fam_file <- read.table(file = "data/MMP.fam")

# open .SSD and .info files we just created
SSD.info <- Open_SSD('data/MMP.SSD', 'MMP.info')

# create null model based on phenotypes (and co-variates if you have any)
set.seed(100)
Null_Model <- SKAT_Null_Model(formula = fam_file$V6 ~ 1, out_type="D")

# perform SKAT on all sets of variants (we wont run this because it takes too long)
#All_SKAT_Data  <- SKAT.SSD.All(SSD.INFO = SSD.info, obj = Null_Model) 
```

SKAT.SSD.All returns a dataframe that contains `SetID`, p-values (`P.value`), the number of markers in the SNP sets 
(`N.Marker.All`), and the number of markers to test for an association after excluding non-polymorphic or high
missing rates markers (`N.Marker.Test`).
```{r}
# R

# Checkout SKAT results
# All_SKAT_Data$results[1:20,]

# open saved file
All_SKAT_Data <- read.table(file = "data/SKAT_all-pvals.results", header =TRUE)
All_SKAT_Data[1:20,]
```

The results are sorted via SetID, and not p-value. Thus, to find which genes are most highly associated with our
phenotype/disease outcome we must sort the results by p-value.
```
# R

# sort SKAT results by p-value
head(All_SKAT_Data[order(All_SKAT_Data$P.value),])

# save SKAT results
# write.table(x = All_SKAT_Data$results, file = "data/SKAT_all-pvals.results", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
```

### Apply multiple testing correction

We have just performed association tests on > 1000 genes... We need to correct for multiple testing. There
are many ways to do this, including Bonferroni correction, False discovery rate and resampling. We will choose
a "medium" stringency multiple correction adjustment, the false discovery rate.
```
# R

# load fdrtool library
library(fdrtool)

# calculate q-values
qvals <- fdrtool(All_SKAT_Data$P.value, statistic = "pvalue", cutoff.method="fndr", plot = FALSE)
All_SKAT_Data$Q.value <- qvals$qval

# sort SKAT results by q-value
All_SKAT_Data <-All_SKAT_Data[order(All_SKAT_Data$P.value),]
All_SKAT_Data[1:20,]

#save SKAT results
# write.table(x = All_SKAT_Data$results, file = "data/SKAT_all-qvals.results", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
```

We find that if we choose a 30% false-discovery rate cut-off that we find that a lot of genes are significantly associated with
the phenotype. This doesn't make much sense? This is due to a skewed p-value distrubution from these samples which arises from 
the distribution of minor allele counts in the Million Mutation project strains which exponentially decreases from 1 to N.
```
# R

# load plyr library so we can count how many variants are in each snp set
library(plyr)

# load snp set ID (.SSID) file
SSID <- read.table(file = "data/MMP.SSID", stringsAsFactors = FALSE)

# count how many variants are in each snp set
MACs <- count(SSID[,1])

# plot this as a histogram
pdf("MAC_hist.pdf")
MAC_hist <- hist(x <- (MACs[MACs$freq < 50,2]), breaks = 50, xlab = "Minor Allele Count", ylab = "Number of Genes", main = "Allele distribution")
dev.off()
```

View the file on your local machine by opening a terminal and typing:
```
# Bash Shell

cd Desktop
# scp -i <Path_to_your_key_pair_file> ubuntu@<your_ec2_number>:/home/ubuntu/MAC_hist.pdf .
```
Next open the .pdf on your desktop from Finder/Explorer.


Thus, to use false-discovery rate we need to decide on a minimum minor allele count for SNP sets we want 
to include our analysis. This is important beyond correctiong for multiple comparisons, for example, how
meaningful is a gene that only has a single variant in the population being looked at? How much confidence
do we have that it would be associated (or not associated) with the individuals phenotype/disease outcome?  
What is the minimum number of variants in SNP set to include it in the analysis? Unfortuneatly, this is a 
question that as of yet has no solid answer and researchers must make their own decision on the minimum
minor allele count for a SNP set must be. 


From the histogram of the minor allele count per gene for this dataset, a minimum minor allele count of 7 
seems to be reasonable. Thus, we will reduce our SNP set ID (`.SSID`) file to only those genes which
meet that criteria.
```
# R

# load .SSID file and give meaningful column names
SSID <- read.table(file = "data/MMP.SSID", header=FALSE)
colnames(SSID) <- c("gene", "variant")
head(SSID)

# make a list of the variants to keep (i.e. plyr::count the number of occurrences of each gene)
SSID_count <- count(SSID[1])
head(SSID_count)
  
# make SSID file with only those variants (i.e. reduce the SSID list to only those with >= 6 variants)
genes_to_include <- data.frame(SSID_count[SSID_count$freq >= 5, 1])

# give the data frame meaningful column names  
colnames(genes_to_include) <- ("gene")
head(genes_to_include)

# reduce SSID to contain only those genes from the list genes_to_include via dplyr::semi_join
reduced_SSID <- semi_join(SSID, genes_to_include)
head(reduced_SSID)

# save SSID file
write.table(x = reduced_SSID, file = "data/MMP-reduced.SSID", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

Now can run SKAT with this reduced SNP set ID (`.SSID`) file:
```
# R
## CAUTION: LONG-RUN TIME PROCESS

# load SKAT library (to do association analysis)
library(SKAT)

# Generate a SNP set data file (`.SSD`) and an `.info` file from binary plink formated data files using user specified SNP sets.
Generate_SSD_SetID('data/MMP.bed', 'data/MMP.bim', 'data/MMP.fam', 'data/MMP-reduced.SSID', 'data/MMP-reduced.SSD', 'MMP-reduced.info')

# open .SSD and .info files we just created
SSD.info_reduced <- Open_SSD('data/MMP-reduced.SSD', 'MMP-reduced.info')

# null model doesn't take in info from .SSID, so we don't need to run that again
# perform SKAT on all sets of variants (we wont run this because it takes too long)
# All_SKAT_Data_reduced  <- SKAT.SSD.All(SSD.INFO = SSD.info_reduced, obj = Null_Model) 

# open saved file
All_SKAT_Data_reduced <- read.table(file = "data/SKAT_all_reduced-pvals.results", header =TRUE)
All_SKAT_Data_reduced[1:20,]

# sort SKAT results by p-value
head(All_SKAT_Data_reduced[order(All_SKAT_Data_reduced$P.value),])

# calculate q-values
qvals <- fdrtool(All_SKAT_Data_reduced$P.value, statistic = "pvalue", cutoff.method="fndr", plot = FALSE)
All_SKAT_Data_reduced$Q.value <- qvals$qval

# sort SKAT results by q-value
All_SKAT_Data_reduced <- All_SKAT_Data_reduced[order(All_SKAT_Data_reduced$P.value),]
All_SKAT_Data_reduced[1:20,]

#save SKAT results
# write.table(x = All_SKAT_Data_reduced$results, file = "data/SKAT_all_reduced.results", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
```

Now we can see that only a handful of genes are associated with the phenotype - a much more reasonable number. So
what's next? The gold standard in human studies would to be to perform an independent data collection and 
analysis on another sample of the population and if any genes show up as signficant in both sets of analyses
then we would have a fair anount confidence that their disruption may cause the phenotpye/disease outcome.


In genetically manipuable organisms you could instead experimentally confirm that mutations in these 
genes are causative for the phenotpye/disease outcome via Crispr/Casp, genetic rescue, RNAi 
and/or mapping experiments.


**References:**
<br/>**1.** Kumar P, Henikoff S, Ng PC. 2009. Predicting the effects of coding non-synonymous variants on protein function using the SIFT algorithm. Nat Protoc 4: 1073-1081.
<br/>**2.** Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics, 81.
<br/>**3.** Ramensky V, Bork P, Sunyaev S. 2002. Human non-synonymous SNPs: server and survey. Nucleic Acids Res 30: 3894-3900.
<br/>**4.** Thompson O, Edgley M, Strasbourger P, Flibotte S, Ewing B, Adair R, Au V, Chaudry I, Fernando L, Hutter H et al. 2013. The Million Mutation Project: A new approach to genetics in Caenorhabditis elegans. Genome Res doi:gr.157651.113 [pii] 10.1101/gr.157651.113.
<br/>**5.** Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. American journal of human genetics 89: 82-93.
