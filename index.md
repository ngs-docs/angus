# Next-Gen Sequence Analysis Workshop (2017)

These are the schedule and classroom materials for the
[ANGUS workshop at UC Davis](http://ivory.idyll.org/dibsi/ANGUS.html),
which will run from June 26th to July 8th, 2017.

This workshop runs under a [Code of Conduct](code-of-conduct.html). Please
respect it and be excellent to each other!

Twitter hash tag: [#ngs2017](https://twitter.com/search?f=tweets&q=%23ngs2017&src=typd)

## The main workshop materials

### Monday, Day 1: Introduction

7:30pm - 9pm, Valley Hall.

Introductions!

[Booting a cloud computer from Jetstream!](jetstream/boot.html)

[Running BLAST at the command line](running-command-line-blast.html)

Assessment.

Student 30-second introductions.

### Tuesday, Day 2: Command line.

9am lecture: [C. Titus Brown, UC Davis - "On Biology and Data Analysis."](https://osf.io/nsab3/)

Morning: 10:15am.
* [Running large scale BLASTs and output.](running-blast-large-scale.html)

Afternoon: 1:15pm - 3pm.
* [Visualizing BLAST score distributions in RStudio](visualizing-blast-scores-with-RStudio.html)
* [Review and explore: Command line UNIX, and R/RStudio](command-line-and-rstudio.html)

Evening: 7pm-8:30: student presentations and questions! Ice cream social.

### Wednesday, Day 3: Mapping, aligning, and sequence files.

9am lecture:
[Melissa Wilson Sayres](http://www.wilsonsayreslab.org/), Arizona State University - ["Sex-biased genome evolution"](https://osf.io/czj42/)

Morning and afternoon: 10:15am-4pm.
* [Quality trimming your reads.](quality-trimming.html)
* [Mapping and samtools and variant calling.](variant-calling.html)
* [High Throughput sequencing. What could go wrong?](https://github.com/wltrimbl/whatcouldgowrong)

(Evening: free time - [Wed Farmers' market!](http://www.davisfarmersmarket.org/))

### Thursday, Day 4: Introduction to counting and differential expression.

9am lecture: [Chris Hamm](https://butterflyology.github.io/about-me.html), Monsanto - [Why are there so many butterflies?](https://osf.io/3j5yf/)

Morning: 10:15am-noon
* [Introduction to R and Data Analysis](introduction-to-R-and-dataframes.html)

Afternoon and evening: 1:15pm-3pm, 7-8:30.
* [Short reads and counting for variant calling and differential expression](counting.html)

### Friday, Day 5: Reference independent analyses: assembly

9am lecture: [Adrienne Roeder](http://roeder.wicmb.cornell.edu/), Cornell - [Reaching biological conclusions from RNA-seq: the good, the bad, and the ugly](https://osf.io/qz3m6/)

Morning: 10:15am-noon.
* [Microbial genome assembly](genome-assembly.html)

Afternoon: 1-3pm.
* [Prokka genome annotation](prokka_genome_annotation.html)

Evening: (also running on Saturday afternoon; you only need to come to one.)
* [Introduction to automation.](introduction-to-automation.html)

### Saturday, Day 6: Automation and repeatability.

(Morning: free time / Saturday market)

Afternoon (1-4pm): [Introduction to automation.](introduction-to-automation.html)

Note! this will be a repeat of Friday evening! You only need to come to one.

### Sunday, Day 7: Day of rest!

### Monday, Day 8: Genomes, and GWAS, and k-mers.

9am lecture: [Erich Schwarz](https://mbg.cornell.edu/people/erich-schwarz), Cornell - ["Assembling and biologically interpreting nematode genomes"](https://osf.io/9yu76/)

Morning: 10:15am-lunch, in lions & tigers & bears.
* [Reference independent analyses with k-mers; sourmash.](kmers-and-sourmash.html) - Luiz in Lions, Phil in Tigers, and Titus in Bears.

Lunch will be local, with a food truck

Afternoon: 1pm onwards, in 1030 Valley.
* [Variant calling and big genomes: GATK.](GATK_pipeline.html) (Tamer Mansour)
* [Genome Wide Association Studies.](GWAS.html)
* [Meta-analysis of GWAS studies](meta_GWAS.html) (Shannon Joslin)
   
Evening: Ice cream and demos! 1030 Valley.
* 6:30pm onwards: ice cream!
* 7pm-9pm: demonstrations and eye candy (1030 Valley)
  - CyVerse and what an allocation request looks like.
  - [Jupyter Notebook, R and Python for data science.](Jupyter-Notebook-Notes.html)
  - where should I put my data? One option: the [Open Science Framework](the_osf.html)
  - [GitHub](github.html)
  - [where do I find the data? NCBI, ENSEMBL, ENA; how to get FASTQ out of NCBI.](database_resources.html)
  - & other things

### Tuesday, Day 9: Nanopore sequencing; and a panel discussion!

9am tutorials:
* [Assessing & assembling nanopore data.](analyzing_nanopore_data.html) (Lisa Cohen and Jon Badalamenti)
* 11am: panel discussion: The Future of Biology, Bioinformatics, and Humanity.
   * [Titus' talk](https://osf.io/zbqtv/)
   * [Adrienne's talk](https://osf.io/8mb2y/)
   * [Erich's talk](https://osf.io/5sr93/)

4pm onwards: [July 4th celebration at Community Park](http://cityofdavis.org/city-hall/city-manager-s-office/community-events/fourth-of-july). We will have a sun cover and ice chests!

### Wednesday, Day 10: RNAseq.

9am lecture: [Megan Dennis](http://www.dennislab.org/), UC Davis - ["Complex genomic variation and its role in human evolution and disease"](https://osf.io/9d5ge/)

Morning: 10am-noon
* [RMarkdown](rmarkdown_rnaseq.html) - in lions/tigers/bears rooms.

Wed 1:15-4pm: [breakout sessions](https://hackmd.io/CwZg7AjAhgnDBGBaAZgDmRRwAmAGKiqUIAxovKhAKzI4xQ4hA===?view)
* meet in 1030 Valley at 1:15pm to choose breakout sessions - see [hackmd document](https://hackmd.io/CwZg7AjAhgnDBGBaAZgDmRRwAmAGKiqUIAxovKhAKzI4xQ4hA===?view) for breakout sessions!
* move to breakout rooms on special topics.

5pm onwards: Wed market.

### Thursday, Day 11: RNAseq.

9am lecture: [Michael I Love](https://mikelove.github.io/), UNC Chapel Hill - ["Statistics and bias correction in RNAseq differential expression analysis"](https://osf.io/gbjhn/)

Morning (10am-noon):
* Lions: [De novo RNAseq assembly](assembly-trinity.html) - Lisa
* Tigers: [DESeq2](deseq2-asthma.html) - Mike and Rob
* Bears: [ChIP-seq](chip-seq.html) - Fotis and Titus

Afternoon (1:15pm-3:15pm):
* Lions: [Transcriptome annotation](dammit_annotation.html) - Camille
* Tigers: [De novo RNAseq assembly](assembly-trinity.html) - Lisa
* Bears: [DESeq2](deseq2-asthma.html) - Mike and Rob

Evening: (7pm, 1030 Valley)
* Q&A with Mike and Rob; optional.

### Friday, Day 12: RNAseq.

9am lecture: [Robert Patro](http://www.robpatro.com/redesign/), Stony Brook University - "Don't count on it: Pragmatic and theoretical concerns and best practices for mapping and quantifying RNA-seq data"

Morning (10am-noon):
* Lions: [DESeq2](deseq2-asthma.html) - Mike and Rob
* Tigers: [ChIP-seq](chip-seq.html) - Titus
* Bears: [De novo RNASeq assembly](assembly-trinity.html) - Tessa

Afternoon (1:15pm-3:15pm):
* Lions: [ChIP-seq](chip-seq.html) - Fotis
* Tigers:  [Transcriptome annotation](dammit_annotation.html) - Camille
* Bears:  [Transcriptome annotation](dammit_annotation.html) - Lisa

Assessment.

### Saturday, last day: Summary

9am lecture: C. Titus Brown, UC Davis - "Effectively infinite: next steps in Data Intensive Biology."

Further resources; continuing your learning.

Questions and answers.

Course review, input, and feedback.
