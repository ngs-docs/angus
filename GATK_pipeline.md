# Variant calling pipeline for a mammalian genome 

We will run a variant calling pipeline using [Genome Analysis Toolkit (GATK)](https://software.broadinstitute.org/gatk/) using a subset sample of dog WGS as a representative 
to large mammalian genomes.

## Getting started

[Start up an m1.medium instance running Ubuntu 16.04 on Jetstream.](jetstream/boot.html)

log in, and then make & change into a working directory:

      mkdir ~/GATK_tutorial && cd ~/GATK_tutorial
      
## Download trimmed Fastq files 

      wget https://de.cyverse.org/dl/d/3CE425D7-ECDE-46B8-AB7F-FAF07048AD42/samples.tar.gz
      tar xvzf samples.tar.gz
      rm samples.tar.gz

> Quick notes about read trimming for variant calling:
> 1. Trimming is data loss so be careful.
> 2. Sequence trimming is complementary to variant filtration
> 3. Sources of errors: 
>    a) The call is suspcious ==> low quality score (variant filtration is better than quality trimming)
>    b) Technical problems (e.g. sequencing chemistry or physics) ==> systematic erros (can be removed by careful kmer based trimming but GATK recalibration is an alternative) 
> 4. Very mild quality trimming: SLIDINGWINDOW:4:2 ==> this means that the [Base call accuracy is ~ 40%](https://en.wikipedia.org/wiki/Phred_quality_score)

## Mapping

1.  Install [bwa](http://bio-bwa.sourceforge.net/bwa.shtml):

        cd
        curl -L https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2/download > bwa-0.7.15.tar.bz2
        tar xjvf bwa-0.7.15.tar.bz2
        cd bwa-0.7.15
        make
      
        sudo cp bwa /usr/local/bin
      
        echo 'export PATH=$PATH:/usr/local/bin' >> ~/.bashrc
        source ~/.bashrc

2.  change into a working directory:

        mkdir ~/GATK_tutorial && cd ~/GATK_tutorial

3.  download and prepare the reference for mapping

        wget https://de.cyverse.org/dl/d/A9330898-FC54-42A5-B205-B1B2DC0E91AE/dog_chr5.fa.gz
        gunzip dog_chr5.fa.gz
        bwa index -a bwtsw dog_chr5.fa

4.  Add [Read group information](Read_group_info.md) and do mapping

> Read group information is typically added during this step, but can also be added or modified after mapping using Picard AddOrReplaceReadGroups.

      for R1 in *_R1_001.pe.fq.gz;do
        SM=$(echo $R1 | cut -d"_" -f1)                                          ##smaple ID
        LB=$(echo $R1 | cut -d"_" -f1,2)                                        ##library ID
        PL="Illumina"                                                           ##platform (e.g. illumina, solid)
        RGID=$(zcat $R1 | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)       ##read group identifier 
        PU=$RGID.$LB                                                            ##Platform Unit
        echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"

        R2=$(echo $R1 | sed 's/_R1_/_R2_/')
        echo $R1 $R2
        bwa mem -t 4 -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" dog_chr5.fa $R1 $R2 > ${R1%_R1_001.pe.fq.gz}.sam
      done


