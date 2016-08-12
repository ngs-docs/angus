============================
Population Genetics Tutorial
============================

========= 
Exercise
========= 
Before you guys got here
~~~~~~~~~~~~~~~~~~~~~~~~

Started with data from: "Genomic islands of speciation separate cichlid ecomorphs in an East African crater lake", Malinsky et al 2015. 

Downloaded VCF from http://datadryad.org/resource/doi:10.5061/dryad.770mc
	- http://datadryad.org/bitstream/handle/10255/dryad.101389/Massoko_Dryad_VCF_final.vcf.gz
	- These data had been filtered for quality
	- And only variable sites had been retained
	- And phased using the program `BEAGLE`, which relies on linkage disequilibrium to phase haplotypes

Made the VCF smaller so we could analyze it in this lifetime: 36 individuals and no indels.::

	vcftools --gzvcf Massoko_Dryad_VCF_final.vcf.gz --keep inds_to_keep.txt --stdout --recode --recode-INFO-all --remove-indels | gzip -c > Massoko_Dryad_VCF_final_subset_noIndels.vcf.gz

Made the VCF smaller still to remove low frequency sites and then local linkage disequilibrium. We will use these files for many of our analyses.::

	vcftools --gzvcf Massoko_Dryad_VCF_final_subset_noIndels.vcf.gz --maf 0.05 --max-maf 0.95 --stdout --recode --recode-INFO-all | gzip -c > Massoko_Dryad_VCF_final_subset_noIndels_maf05.vcf.gz
	vcftools --gzvcf Massoko_Dryad_VCF_final_subset_noIndels_maf05.vcf.gz --thin 1000 --stdout --recode --recode-INFO-all | gzip -c > Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.vcf.gz

Used the thinned VCF to make input files for phylogenetic inference and population structure analyses.::

	python vcf_to_phy.py --infile Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.vcf.gz --thin 5
	python vcf_to_geno.py
	EIG-6.1.3/bin/convertf -p convert.par
	cat Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.map | awk -F'\\\s+' '{print $1,$2, $3,$4}' > map
	mv map Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.map
	plink --file Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K --recode
	mv plink.ped Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.ped
	mv plink.map Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.map


Now you start
~~~~~~~~~~~~~

Start your instance with at least 50 Gb.::

	sudo apt-get update
	sudo apt-get install build-essential
	sudo apt-get install ruby git
	ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"
	export PATH="/home/ubuntu/.linuxbrew/bin:$PATH"
	export MANPATH="/home/ubuntu/.linuxbrew/share/man:$MANPATH"
	export INFOPATH="/home/ubuntu/.linuxbrew/share/info:$INFOPATH"
	export LD_LIBRARY_PATH="/home/ubuntu/.linuxbrew/lib/"
	brew tap homebrew/science

Install many useful packages.::

	brew install zlib
	brew install gsl
	brew install homebrew/science/raxml
	brew install homebrew/science/vcftools
	brew install homebrew/science/openblas

Install some other software.::

	wget https://data.broadinstitute.org/alkesgroup/EIGENSOFT/EIG-6.1.3.tar.gz
	tar -xvzf EIG-6.1.3.tar.gz
	wget https://www.genetics.ucla.edu/software/admixture/binaries/admixture_linux-1.3.0.tar.gz
	tar -xzvf admixture_linux-1.3.0.tar.gz

You can download the smaller data set and ancillary files from here.::

	wget https://www.dropbox.com/s/ra4yqix0jfe1fgn/tutorial_files.tar.gz
	tar -xzvf tutorial_files.tar.gz
	cd tutorial_files

Calculate nucleotide diversity (:math:`\pi`). Use `VCFtools` to figure out how to calculate it. We want to calculate it for 'benthic' and 'littoral' morphs separately.::

	vcftools --gzvcf Massoko_Dryad_VCF_final_subset_noIndels.vcf.gz --keep littoral.txt --window-pi 100000 --out littoral_pi
	vcftools --gzvcf Massoko_Dryad_VCF_final_subset_noIndels.vcf.gz --keep benthic.txt --window-pi 100000 --out benthic_pi

Calculate linkage disequilibrium.::

	vcftools --gzvcf Massoko_Dryad_VCF_final_subset_noIndels.vcf.gz --keep littoral.txt --ld-window-bp 500000 --chr scaffold_0 --hap-r2 --out littoral_scaffold_0_ld --min-r2 0.001
	vcftools --gzvcf Massoko_Dryad_VCF_final_subset_noIndels.vcf.gz --keep benthic.txt --ld-window-bp 500000 --chr scaffold_0 --hap-r2 --out benthic_scaffold_0_ld --min-r2 0.001

Summarize linkage disequilibrium data files so that they are smaller and easier to plot.::

	python summarize_ld.py --infile littoral_scaffold_0_ld.hap.ld --win 10
	python summarize_ld.py --infile benthic_scaffold_0_ld.hap.ld --win 10

Calculate :math:`F_{ST}` between benthic and limnetic forms.::

	vcftools --gzvcf Massoko_Dryad_VCF_final_subset_noIndels.vcf.gz --weir-fst-pop littoral.txt --weir-fst-pop benthic.txt --fst-window-size 100000 --out benthic_limnetic_fst

Make a phylogeny.::

	raxmlHPC-PTHREADS -T 12 -m GTRGAMMA -n Massoko -s Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K_thin5.phy -p 123 -o A_calliptera_Chitimba,A_calliptera_Bua,A_calliptera_Chizumulu

Run `ADMIXTURE` for up to 6 populations.::

	~/admixture_linux-1.3.0/admixture Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.ped 1
	~/admixture_linux-1.3.0/admixture Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.ped 2
	...

Run `smartpca`.::

	~/EIG-6.1.3/bin/smartpca -p Massoko_smartpca.par > Massoko_smartpca.out

Now that we have all the different pieces, let's start to plot the data and see what we find. Put all the results into one folder and download them locally so that we can plot and visualize them using `R`.

Just to be sure, here are all the files you should have. Should things be taking too long, you can borrow my results that I generated earlier: https://www.dropbox.com/s/czrru76ku2kqwt2/results.tar.gz?dl=0 

- benthic_limnetic_fst.windowed.weir.fst
- benthic_pi.windowed.pi
- benthic_scaffold_0_ld.hap.ld_summary_window10
- littoral_pi.windowed.pi
- littoral_scaffold_0_ld.hap.ld_summary_window10
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.1.P
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.1.Q
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.2.P
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.2.Q
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.3.P
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.3.Q
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.4.P
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.4.Q
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.5.P
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.5.Q
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.6.P
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.6.Q
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.eval
- Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.evec
- RAxML_bestTree.Massoko

We have the following data types.

#. Genetic diversity.
#. Genetic differentiation. (:math:`F_{ST}`)
#. Decay of linkage disequilibrium.
#. A tree.
#. PCA results.
#. `ADMIXTURE` population clustering results.

We will be using `R` to plot all these results. I will get you started on how to start thinking about some of these. I would recommend setting your working directory to be the directory that has all your results. For example,::

	setwd("/Users/sonal/Desktop/results/")

Note that this is generally considered bad programming practice for scripts that will be publicly shared, but it is convenient when doing exploratory data analysis.

Genetic diversity
~~~~~~~~~~~~~~~~~
To load the genetic diversity results,::

	b = read.table("benthic_pi.windowed.pi", header=T)
	l = read.table("littoral_pi.windowed.pi", header=T)

Look at how the data is structured and summarize it quickly,::

	head(b)
	summarize(b)

To answer some of the questions below, it might be useful to combine across both data-frames::

	x = merge(b, l, by=c("CHROM", "BIN_START", "BIN_END"))

To answer some of the questions below, it might be useful to combine across both data-frames in another way::

	all = data.frame(c(b$PI, l$PI), c(rep("benthic", nrow(b)), c(rep("littoral", nrow(l)))))
	names(all) = c("PI", "MORPH")

You might want to also explore the following functions to answer the questions::

	cor.test()
	boxplot()
	aov() 
        # if you store the results of aov() in a variable and then run summary() on the variable, you get more info


Some questions:

#. What is min, max, and mean levels of genetic diversity in each morph?
#. Is genetic diversity between the two morphs significantly different?
#. Why might genetic diversity be higher in one morph than another? How could you test this?
#. How correlated is genetic diversity between the two morphs?
#. Why would genetic diversity be correlated between the two morphs?

Genetic differentiation
~~~~~~~~~~~~~~~~~~~~~~~
To load the genetic differentiation results,::

	fst = read.table("benthic_limnetic_fst.windowed.weir.fst", header=T)

To select rows that have certain values,::

	x = fst[fst$CHROM == 'scaffold_0', ]
	x = fst[fst$WEIGHTED_FST >= 0.1, ]

You might want to explore the functions::

	dim()
	nrow()

Which allow you to quickly figure out how big these dataframes are.

Some questions:

#. What is the mean :math:`F_{ST}` between these two morphs?
#. Is there a correlation between the number of variants in a window and :math:`F_{ST}`? If so, it would suggest we need to be cautious of these results.
#. In this paper, the authors emphasize the importance of genomic regions that are highly differentiated. How many 100 kb windows are differentiated above :math:`F_{ST}`>0.1? :math:`F_{ST}`>0.2? :math:`F_{ST}`>0.3?
#. Plot :math:`F_{ST}` along BIN_START on scaffold_15. The authors originally identified 10 peaks (see Fig. 3D) that are highly differentiated. How many do you identify? Why might our results be different?
#. How might you determine if windows with high :math:`F_{ST}` are significant?

Decay of Linkage Disequilibrium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can read in the tables for linkage disequilibrium just like you did for nucleotide diversity.

Having done that, we can now plot the data. Plot distance on the x-axis and :math:`r^2` on the y-axis (a measure of linkage disequilibrium that looks at the correlation coefficient between pairs of loci -- higher values means that two loci "travel" together more than you would expect under random assortment).

Try plotting both morphs at once. You will want to use the `points()` function.

Some questions:

#. Do the two morphs have different decay patterns? 
#. A key aspect of linkage disequilibrium is how quickly it decays. At what physical distance is the level of linkage disequilbrium halved? You can estimate this visually or using R.
#. These points are very very noisy. How might you do this exercise again to reduce some of this noise? If you have time, try it!

Plot the phylogeny
~~~~~~~~~~~~~~~~~~
To plot the phylogeny, you will need to install the library ape.::

	install.packages("ape")
	library(ape)

Then, you can read in and plot tree.::

	t = read.tree("RAxML_bestTree.Massoko")
	# makes the tree easier to visualize by ladderizing it
	t = ladderize(t)
	plot(t)

Some questions:

#. What do you think is going on with the "small" morph?
#. Looking at this tree, would you say that the "littoral" and "benthic" morphs are differentiated? Why or why not?
#. Before we use this tree for any formal analysis, what else might you want to check about the tree?

Plot the PCA
~~~~~~~~~~~~
To read in the PCA data::

	d = read.table("Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.evec")

Note that the eval file has the data we would need to calculate the eigenvalues for each PCA axis.

Look at the data file using `head()` -- how is it structured? What does each column mean? 

You can plot it by::

	plot(d$V2, d$V3, col=as.factor(d$V12), pch=16)

This isn't such an informative plot. Why? How would you subset the data to make it more informative? Hint: look at column V12.::

	s = d[d$V12 %in% c("Massoko_benthic", "Massoko_littoral", "Massoko_small"),]

This still isn't as informative as it could be. It likely would have been much more informative if we removed the outgroups before doing the PCA. That said, are these morphs differentiated? How do these results compare to what we saw with the phylogeny? Why might these results be different?

ADMIXTURE results
~~~~~~~~~~~~~~~~~
To read in the `ADMIXTURE` results::

	d1 = read.table("Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.1.Q")
	d2 = read.table("Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.2.Q")
	d3 = read.table("Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.3.Q")
	d4 = read.table("Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.4.Q")
	d5 = read.table("Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.5.Q")
	d6 = read.table("Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.6.Q")

To plot the results::

	par(mfrow=c(6,1), mar=c(1,4,1,1))
	barplot(t(as.matrix(d1)), col=rainbow(1), border=NA)
	barplot(t(as.matrix(d2)), col=rainbow(2), border=NA)
	barplot(t(as.matrix(d3)), col=rainbow(3), border=NA)
	barplot(t(as.matrix(d4)), col=rainbow(4), border=NA)
	barplot(t(as.matrix(d5)), col=rainbow(5), border=NA)
	par(mar=c(3,4,1,1))
	x = barplot(t(as.matrix(d6)), col=rainbow(6), border=NA)
	inds = c(rep('A_cal', 3), rep('Ita', 3), rep('B', 10), rep('L', 10), rep('S', 10))
	mtext(inds, 1, at=x, las=2)

What's going on here? Based on all the results you have seen from the phylogeny, the PCA, and this, how would you characterize the differentiation between these morphs?

========= 
Resources
========= 
Population Genetics Books
~~~~~~~~~~~~~~~~~~~~~~~~~
- Coop's Class Notes: http://cooplab.github.io/popgen-notes/
- Felsenstein's Book: http://evolution.genetics.washington.edu/pgbook/pgbook.html
- Gillespie's *Population Genetics: A Concise Guide*
- Hartl and Clark's *Principles of Population Genetics*
- Nielsen and Slatkin's *An Introduction to Population Genetics*
- Wakeley's *Coalescent Theory*
- Yang's *Computational Molecular Evolution*

Great set of tutorials
~~~~~~~~~~~~~~~~~~~~~~
- http://evomics.org/learning/population-and-speciation-genomics/
- http://grunwaldlab.github.io/Population_Genetics_in_R/Preface.html

Papers on population genomics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- *A framework for variation discovery and genotyping using next-generation DNA sequencing data*, DePristo et al 2010; 10.1038/ng.806
- *Genome sequencing and population genomics in non-model organisms*, Ellegren 2014; 10.1016/j.tree.2013.09.008
- *Genotype and SNP calling from next-generation sequencing data*, Nielsen et al 2011; 10.1038/nrg2986
- *Methods and models for unravelling human evolutionary history*, Schraiber and Akey 2015; 10.1038/nrg4005
- *Population Genomics: Whole-Genome Analysis of Polymorphism and Divergence in Drosophila simulans*, Begun et al 2007; 10.1371/journal.pbio.0050310
- *The power and promise of population genomics: from genotyping to genome typing*, Luikart et al 2003; 10.1038/nrg1226

Software & Programs for working with data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml; great for quality filtering and simple parsing of variants 
- https://github.com/thibautjombart/adegenet/wiki; R package that can parse variant data
- https://vcftools.github.io/index.html; can generate many useful statistics from VCF files
- https://cran.r-project.org/web/packages/PopGenome/index.html; R package that calculates statistics from VCFs, note not very transparent in how it handles missing data
- http://vcf.iobio.io/; allows quick visualization of VCFs
- http://popgen.dk/wiki/index.php/ANGSD; ideal for low coverage data

Learn Python
~~~~~~~~~~~~
- https://github.com/singhal/python_workshop/blob/master/Python.Md
- http://learnpythonthehardway.org/
- https://www.coursera.org/course/pythonlearn
- http://rosalind.info/problems/locations/

Learn R
~~~~~~~
- http://tryr.codeschool.com/
- https://www.coursera.org/learn/r-programming
- https://www.edx.org/course/introduction-r-data-science-microsoft-dat204x-1
- http://swirlstats.com/students.html
- http://r4ds.had.co.nz/

Learn Shell / Unix
~~~~~~~~~~~~~~~~~~
- https://www.codecademy.com/learn/learn-the-command-line
- http://korflab.ucdavis.edu/unix_and_Perl/
- http://www.learnshell.org/

Learn Perl
~~~~~~~~~~
- http://korflab.ucdavis.edu/unix_and_Perl/
