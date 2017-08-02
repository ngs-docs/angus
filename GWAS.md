## Genome Wide Association analysis (GWAS)
To test the association of of a genome-wide set of genetic variants with a given trait.

## install [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)

      cd /usr/local/bin/
      sudo wget https://www.cog-genomics.org/static/bin/plink170725/plink_linux_x86_64.zip
      sudo unzip -o plink_linux_x86_64.zip
      sudo rm -f plink_linux_x86_64.zip
      cd plink-1.07-x86_64/
      echo export PATH=$PATH:$(pwd) >> ~/.bashrc
      source ~/.bashrc

## install [vcftools](https://vcftools.github.io/)

      cd
      git clone https://github.com/vcftools/vcftools.git
      cd vcftools
      ./autogen.sh
      ./configure
      make
      sudo make install
     
## Install R and RStudio
      sudo apt-get update && sudo apt-get install -y gdebi-core r-base r-base-dev
      sudo gdebi -n rstudio-server-1.0.143-amd64.deb

## Make a working directory for the GWAS analysis

      mkdir ~/GWAS && cd ~/GWAS

## Download the sample VCF file and phenotype data
Genotyping of 476840 SNPs in 53 dogs (24 yellow coat and 29 dark coat)

      wget https://de.cyverse.org/dl/d/E0A502CC-F806-4857-9C3A-BAEAA0CCC694/pruned_coatColor_maf_geno.vcf.gz
      gunzip pruned_coatColor_maf_geno.vcf.gz
      wget https://de.cyverse.org/dl/d/3B5C1853-C092-488C-8C2F-CE6E8526E96B/coatColor.pheno

## convert VCF into Plink readable format (map,ped) then Plink binary format (fam,bed,bim)

      vcftools --vcf pruned_coatColor_maf_geno.vcf --plink --out coatColor
      plink --file coatColor --allow-no-sex --dog --make-bed --out coatColor.binary

## create list of alternative alleles

      cat pruned_coatColor_maf_geno.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{{if($3==".")$3=$1":"$2;}print $3,$5;}'  > alt_alleles

## Run a simple association analysis
> --assoc performs a standard case/control association analysis which is a chi-square test of allele frequency.

> By default, the minor allele is coded A1 and tested for being the risk allele. The minor allele is expected to be the alternative allele
but could happen to be the reference one. --reference-allele allow you to use your list of A1 alleles

> --adjust enables correction for multiple analysis and automatically calculates the genomic inflation factor  

      plink --bfile coatColor.binary --make-pheno coatColor.pheno "yellow" --assoc --reference-allele alt_alleles --allow-no-sex --adjust --dog --out coatColor

## Create Manhattan plot

Install qqman package

    sudo Rscript -e "install.packages('qqman',  contriburl=contrib.url('http://cran.r-project.org/'))"

Identify statistical cutoffs

    unad_cutoff_sug=$(tail -n+2 coatColor.assoc.adjusted | awk '$10>=0.05' | head -n1 | awk '{print $3}')
    unad_cutoff_conf=$(tail -n+2 coatColor.assoc.adjusted | awk '$10>=0.01' | head -n1 | awk '{print $3}')

Run the plotting function
```
Rscript -e 'args=(commandArgs(TRUE));library(qqman);'\
'data=read.table("coatColor.assoc", header=TRUE); data=data[!is.na(data$P),];'\
'bitmap("coatColor_man.bmp", width=20, height=10);'\
'manhattan(data, p = "P", col = c("blue4", "orange3"),'\
'suggestiveline = 12,'\
'genomewideline = 15,'\
'chrlabs = c(1:38, "X"), annotateTop=TRUE, cex = 1.2);'\
'graphics.off();' $unad_cutoff_sug $unad_cutoff_conf
```

The top associated mutation is a nonsense SNP in MC1R (c.916C>T) known to control pigment production. The file `coatColor_man.bmp` can be visualized with RStudio or downloaded to your local computer (from local commandline) via: `scp username@ipaddress:/home/username/GWAS/coatColor_man.bmp .` 
