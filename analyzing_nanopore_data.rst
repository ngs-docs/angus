=======================
Analyzing nanopore data
=======================

We sequenced *Escherichia coli* K-12 MG1655 yesterday on the MinION. Conveniently this is the same strain that you assembled earlier in the week. Therefore it should be possible to improve the Illumina assembly with this dataset.

The goals of this tutorial are to:
   *  demonstrate real-time base-calling in the cloud
   *  assess a nanopore run
   *  perform a hybrid assembly using nanopore and Illumina data with SPAdes.

Acquiring E. coli nanopore data
===============================

The run we put on last night is available, we got about 1000 reads. You can download them and take a look: ::

   wget http://microbesng.uk/filedist/nanopore/Ecoli_KBS_2D.fasta
   wget http://microbesng.uk/filedist/nanopore/Ecoli_KBS_complement.fasta
   wget http://microbesng.uk/filedist/nanopore/Ecoli_KBS_template.fasta

In addition, here is a better run with around 8000 'two-direction' reads (http://microbesng.uk/filedist/nanopore/FC20.wf1.9.2D.pass.fasta). Retrieve it from::

   wget http://microbesng.uk/filedist/nanopore/FC20.wf1.9.2D.pass.fasta

You could analyse a) one of them b) both of them c) the two of them combined (hint use *cat* to concatenate two files). Why not do something different from your neighbor and then we can compare notes?

Installing SPAdes
=================

This time rather than building from source, we will use the SPAdes Linux binaries.

Installing SPAdes binary (v3.6.0)::

   wget http://spades.bioinf.spbau.ru/release3.6.0/SPAdes-3.6.0-Linux.tar.gz
   tar xvfz SPAdes-3.6.0-Linux.tar.gz
   export PATH=$PATH:`pwd`/SPAdes-3.6.0-Linux/bin

Check it is working::

   spades.py -h

Get the E. coli MDA Illumina files (these are the same you generated before)::

   wget http://public.ged.msu.edu.s3.amazonaws.com/ecoli_ref-5m-trim.se.fq.gz
   wget http://public.ged.msu.edu.s3.amazonaws.com/ecoli_ref-5m-trim.pe.fq.gz

SPAdes has support for nanopore reads using the --nanopore option which expects FASTA formatted files.

You can give 2D reads (preferred), but I have had good results giving both the 1D and 2D reads. More reads is usually better in this case because SPAdes uses the long nanopore reads for gap closure and repeat resolution, rather than for building the assembly graph.::

   spades.py --sc --12 ecoli_ref-5m-trim.pe.fq.gz -s ecoli_ref-5m-trim.se.fq.gz --nanopore FC20.wf1.9.2D.pass.fasta -o nanopore-ecoli-sc

Genome assessment
=================

First, run QUAST on your new assembly and compare against the reference and the Illumina-only assembly. If you don't have QUAST installed anymore, go back and look at :doc:`assembling-ecoli` for instructions.

How many contigs do you have? What is the N50? Where are the discontiguities (hint: find and look at the diagonal plot).

Installing BWA
==============

We want a version of BWA that supports Oxford Nanopore reads: ::

   git clone https://github.com/lh3/bwa.git
   cd bwa
   make
   export PATH=`pwd`:$PATH
   cd ..

Let's now look in more detail. Map reads back to the assembly and assess coverage evenness. The steps are:

   * indexing the reference genome - in this case the reference genome is our de novo assembly
   * aligning, converting SAM to BAM, then sorting the BAM file
   * indexing the BAM file

The commands to perform these tasks are: ::

   bwa index nanopore-ecoli-sc/scaffolds.fasta
   bwa mem -t4 -x ont2d nanopore-ecoli-sc/scaffolds.fasta FC20.wf1.9.2D.pass.fasta | samtools view -bS - | samtools sort - mapped_reads.sorted
   samtools index mapped_reads.sorted

Download the resulting mapped_reads.sorted.bam, mapped_reads.sorted.bam.bai and nanopore-ecoli-sc/scaffolds.fasta files and open in IGV.

What does it look like? What's the coverage like? Can you spot any problems? What is the Oxford Nanopore error profile? Does it do badly in any regions, which ones? Why?


