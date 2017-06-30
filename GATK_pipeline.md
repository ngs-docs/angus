# Variant calling pipeline for a mammalian genome 

We will run a variant calling pipeline using Genome Analysis Toolkit (GATK) using a subset sample of dog WGS as a representative 
to large mammalian genomes.

## Getting started

[Start up an m1.medium instance running Ubuntu 16.04 on Jetstream.](jetstream/boot.html)

log in, and then make & change into a working directory:

      mkdir ~/GATK_tutorial && cd ~/GATK_tutorial
      
## Download trimmed Fastq files 

      wget https://de.cyverse.org/dl/d/3CE425D7-ECDE-46B8-AB7F-FAF07048AD42/samples.tar.gz
      tar xvzf samples.tar.gz
      rm samples.tar.gz

> Quick notes about read trimming for variant calling 
