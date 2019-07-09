# Using snakemake to automate an RNA-seq workflow

In our previous lesson, we use snakemake to automate the first two rules of our
quality control/quality analysis. Given the things you learned during that 
lesson, we will now automate the rest of the quality control  workflow through
using snakemake. For this lesson, we would like you to work with others in 
your room to make the snakefile. We have provided the basic structure for each 
of the rules, and would like you to fill in the remaining necessary details for 
the full workflow to run on all samples. 

We've added in the rules we created from our previous lesson as well.

Let's make sure we're in a good working directory and we have the necessary
software installed. If you've already ran the installation command on your 
current instance, you don't need to run it again. 

```
cd ~
conda install -y fastqc multiqc trimmomatic
```

```
SAMPLES=['ERR458493', 'ERR458494', 'ERR458495', 'ERR458500', 'ERR458501', 
'ERR458502']

rule all:
    input:
        "" # add the path to the final file 

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

# Add the trimmomatic commands from our trimming lesson. 
rule trim:
    input:
        reads = "",
        adapters = ""
    output: "" # the output should be the trimmed file
    shell: '''
    trimmomatic SE 
    '''

# Use the commands above as a template to fill in these rules, this time
running the analyses on the trimmed reads.
rule fastqc_trimmed:

rule multiqc_trimmed:

```

When complete, this snakemake workflow automates all of the bash steps of our
quality control workflow!

You can generate a dag of your workflow using the following command:

```
snakemake --dag | dot -Tpng > dag.png
```

Remember, "dag" stands for Directed Acyclic Graph. It shows the steps of your
workflow that are executed by Snakemake. This dag can be helpful to 
communicate your workflow to others, or to visualize the decisions snakemake is
making as you are constructing your workflow. 

You can open it to view using the "Files" tab in the RStudio pane and opening
the `dag.png` file you created.  
