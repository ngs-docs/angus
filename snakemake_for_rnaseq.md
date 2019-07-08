# Using snakemake to automate an RNA-seq workflow

In our previous lesson, we use snakemake to automate the first two rules of our
quality control/quality analysis. Given the things you learned during that 
lesson, we will now automate an RNA-seq workflow through read quantification
 using snakemake. For this lesson, we would like you to work with others in 
your room to make the snakefile. We have provided the basic structure for each 
of the rules, and would like you to fill in the remaining necessary details for 
the full workflow to run on all samples. 

We've added in the rules we created from our previous lesson as well.

```
SAMPLES=['ERR458493', 'ERR458494', 'ERR458495', 'ERR458500', 'ERR458501', 
'ERR458502']

rule all:
    input:
        "fastqc_raw/multiqc_report.html"


rule fastqc_raw:
    input: "data/{sample}.fastq.gz"
    output: 
        "fastqc_raw/{sample}_fastqc.html",
        "fastqc_raw/{sample}_fastqc.zip"
    shell:'''
    fastqc -o fastqc_raw {input}
    '''

rule multiqc_raw:
    input: expand("fastqc_raw/{sample}_fastqc.html", sample = SAMPLES)
    output: "fastqc_raw/multiqc_report.html"
    shell:'''
    multiqc -o fastqc_raw fastqc_raw
    '''

rule trim:

rule fastqc_trimmed:

rule multiqc_trimmed:

rule download_transcriptome:

rule index_transcriptome:

rule quantify_reads:
  
```

When complete, this snakemake workflow automates all of the bash steps of our
RNA-seq workflow. 

We could also integrate R scripts into our workflow. The syntax is a little 
funky, but if you would like to see what an R script looks like inside of
Snakemake, you can download this [script](). 
