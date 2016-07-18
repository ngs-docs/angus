==========================================
Next-Gen Sequence Analysis Workshop (2016)
==========================================

.. Warning These documents are not maintained and their instructions may be
	out of date. However the GED Lab does maintain the `khmer protocols
	<http://khmer-protocols.readthedocs.org/>`__ which may cover similar
	topics. See also the `installation instructions for the current version
	of the khmer project <https://khmer.readthedocs.org/en/latest/install.html>`__.

This is the schedule for the `2016 MSU NGS course <http://bioinformatics.msu.edu/ngs-summer-course-2016>`__.


.. , which will run from August 8th to August 19st, 2015.  If you're interested iin seeing the 2015 description, please see `the 2015 announcement <http://bioinformatics.msu.edu/ngs-summer-course-2015>`__.

This workshop has a :doc:`code-of-conduct`.

`Download all of these materials <https://github.com/ngs-docs/angus/archive/2016.zip>`__ or `visit the GitHub repository <https://github.com/ngs-docs/angus/tree/2016>`__.

Meal Times: Breakfast 7-9, `Lunch <_static/McCrary_17_21.pdf>`__ 12-1, Dinner 6-7 (Unless noted below)

===============  =============================================================
Day              Schedule
===============  =============================================================
Monday 8/8       * 1:30pm lecture: `Welcome! <_static/2015Lecture1Welcome.pptx.pdf>`__ (Meg and Matt)
                 * 3pm: Getting Started with AWS :doc:`amazon/index` (Matt)
                 * 4pm: Intro to Linux (Matt) :doc:`MacManes_UNIX`
                 * 7pm: :doc:`running-command-line-blast` (Meg)
                 * 8pm: Instructor and TA Research presentations, socialize

Tuesday 8/9      * Tutorial :doc:`day2`
                 * 9:15am lecture: `Sequencing considerations <_static/2015-lecture2-sequencing.pptx.pdf>`__ (Matt)
                 * 11:00am Assessment (Julie Libarkin or somebody else)
                 * *Lunch at McCrary 12pm - 1pm*
                 * 1:15pm tutorial, `background PDF <_static/MacManes_Trimming.pdf>`__ and :doc:`MacManesTrimming` (Matt)
                 * Evening *firepit social*

Wed 8/10         * 9:15am lecture, mapping and variant calling lecture (Meg)
                 * 10:30am  tutorial, BASH for genomics, :doc:`GenomicsShell` (Amanda)
                 * *Lunch at McCrary 12pm - 1pm*
                 * 1:15pm tutorial, :doc:`variant` (Meg)
                 * 7:15pm lecture, Teach me scripting :download:`final script <files/Still_script.sh>`
                 * 8:30pm student presentations

Thursday 8/11    * 9:15am lecture, `Eukaryotic Genome assembly and analysis <_static/NGS_Schwarz_2015.08.13.pdf>`__ (Shaun)
                 * 10:15am lecture, Genome assembly exercise (Shaun)
                 * *Lunch at McCrary 12pm - 1pm*
                 * 1:15pm Prokaryotic Genome assembly and analysis, :doc:`assembling-ecoli` (Torst)
				 * 8PM Free Time or Extra help sessions or ...

Friday 8/12      * 9:15am, bacterial genome annotation lecture and practical (Torst) :doc:`prokka_genome_annotation`
                 * *Lunch at McCrary 12pm - 1pm*
                 * 1:15pm sourcing NCBI data including SRA :doc:`howe-ncbi` (Meg); intro to R (Amanda)
				 * 6PM: BBQ on the Island
                 * 7:15pm student research presentations (~3 minutes)

Saturday 8/13    * 9:15am, pop gen lecture and practical (Sonal)
                 * Lunch at McCrary 12pm - 1pm
                 * 1:15pm, lecture and practical, `long read sequencing <_static/Torsten_Seemann_LRS.pdf>`__ (Torsten Seeman)
                 * Takeout (Thai??) Dinner on the island

Sunday 8/14      * Free Day
                 * Brunch at McCrary 12pm - 1pm*
                 * *BBQ Dinner on the island*

Monday 8/15      * 9:15am lecture/tutorial, Transcriptome assembly and evaluation :doc:`MacManes_Trinity` (Matt)
                 * 1:15pm `Mapping / Transrate <https://github.com/ngs-docs/angus/blob/2015/transrate.rst>`__ (Rob)
				 * 7pm: twitter, blogging and bioinformatics by way of social media (All)


Tuesday 8/16     * 9:15am transcriptome read counting lecture and practical (Rob - Salmon)
                 * 10:30am: `R tutorial <https://github.com/ngs-docs/angus/blob/2015/R_Introductory_tutorial_2015.md>`__ (Meg and Ian)
                 * 1:15pm: lecture, `mRNAseq differential expression <_static/NGS2015_RNAseq_2.pptx>`__ (Ian) and lecture, `mRNA stats <_static/NGS2015_RNAseq_ID_1.pptx>`__
                 * 7pm: Journal Club - need article idea


Wed 8/17         * 9:15am lecture, :doc:`functional_annotation`
                 * 10:15am Assembly quality assessment (Transdecoder) (Meg & Matt)
                 * 1:00pm practical, :doc:`analyzing_nanopore_data` (Nick Loman)
                 * Evening, Ask the Expert  (All)

Thursday 8/18    * 9:15am activity, :doc:`MacManes_adventure` (Matt) and Prokka (Torsten)
                 * 1:15pm Assessment (Julie)
                 * 2pm: GitHub practical: :doc:`CTB-github`
                 * *BBQ Dinner on the island*
                 * social

Friday 8/19      * 9:15-9:45 closing lecture (Meg & Matt)
                 * 10am discussion about class; more stuff


                 * Links:
                   `Opinionated guides to NGS <http://davis-assembly-masterclass-2013.readthedocs.org/en/latest/outputs/index.html>`__ /
                   `Software Carpentry <http://software-carpentry.org>`__

===============  =============================================================

Dramatis personae
=================

Instructors:

* Meg Staton
* Matt MacManes
* Ian Dworkin
* Rob Patro
* Torst Seeman
* Sonal Singhal
* Shaun Jackman



TAs:

* Amanda Charbonneau
* Lisa Cohen
* Ming Tang


She Who Drives Many Places:

* Kate MacManes

Papers and References
=====================

Books
-----

* `Practical Computing for Biologists <http://practicalcomputing.org/>`__

  This is a highly recommended book for people looking for a systematic
  presentation on shell scripting, programming, UNIX, etc.

RNAseq
------

* `Differential gene and transcript expression analysis of RNA-seq
  experiments with TopHat and Cufflinks
  <http://www.ncbi.nlm.nih.gov/pubmed/22383036>`__, Trapnell et al.,
  Nat. Protocols.

  One paper that outlines a pipeline with the tophat, cufflinks, cuffdiffs and
  some associated R scripts.

* `Statistical design and analysis of RNA sequencing
  data. <http://www.ncbi.nlm.nih.gov/pubmed/20439781>`__, Auer and
  Doerge, Genetics, 2010.

* `A comprehensive comparison of RNA-Seq-based transcriptome analysis from reads to differential gene expression and cross-comparison with microarrays: a case study in Saccharomyces cerevisiae. <http://www.ncbi.nlm.nih.gov/pubmed/?term=22965124>`__ Nookaew et al., Nucleic Acids Res. 2012.

* `Challenges and strategies in transcriptome assembly and differential gene expression quantification. A comprehensive in silico assessment of RNA-seq experiments <http://www.ncbi.nlm.nih.gov/pubmed/?term=22998089>`__ Vijay et al., 2012.

* `Computational methods for transcriptome annotation and quantification using RNA-seq <http://www.ncbi.nlm.nih.gov/pubmed/21623353>`__, Garber et al., Nat. Methods, 2011.

* `Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. <http://www.ncbi.nlm.nih.gov/pubmed/?term=20167110>`__, Bullard et al., 2010.

* `A comparison of methods for differential expression analysis of RNA-seq data <http://www.biomedcentral.com/1471-2105/14/91>`__, Soneson and Delorenzi, BMC Bioinformatics, 2013.

* `Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. <http://www.ncbi.nlm.nih.gov/pubmed/?term=22872506>`__, Wagner et al., Theory Biosci, 2012.  Also see `this blog post <http://blog.nextgenetics.net/?e=51>`__ explaining the paper in detail.

Computing and Data
------------------

* `A Quick Guide to Organizing Computational Biology Projects <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000424>`__, Noble, PLoS Comp Biology, 2009.

* `Willingness to Share Research Data Is Related to the Strength of the Evidence and the Quality of Reporting of Statistical Results <http://software-carpentry.org/blog/2012/05/the-most-important-scientific-result-published-in-the-last-year.html>`__, Wicherts et al., PLoS One, 2011.

* `Got replicability? <http://econjwatch.org/articles/got-replicability-the-journal-of-money-credit-and-banking-archive>`__, McCullough, Economics in Practice, 2007.

Also see this great pair of blog posts on `organizing projects <http://nicercode.github.io/blog/2013-04-05-projects/>`__ and `research workflow <http://carlboettiger.info/2012/05/06/research-workflow.html>`__.

Links
=====

Humor
-----

* `Data Sharing and Management Snafu in 3 Short Acts <http://www.youtube.com/watch?v=N2zK3sAtr-4&feature=youtu.be>`__

Resources
---------

* `Biostar <http://biostars.org>`__

  A high quality question & answer Web site.

* `SEQanswers <http://seqanswers.com/>`__

  A discussion and information site for next-generation sequencing.

* `Software Carpentry lessons <http://software-carpentry.org/4_0/index.html>`__

  A large number of open and reusable tutorials on the shell, programming,
  version control, etc.

Blogs
-----

* http://www.genomesunzipped.org/

  Genomes Unzipped.

* http://ivory.idyll.org/blog/

  Titus's blog.

* http://bcbio.wordpress.com/

  Blue Collar Bioinformatics

* http://massgenomics.org/

  Mass Genomics

* http://blog.nextgenetics.net/

  Next Genetics

* http://gettinggeneticsdone.blogspot.com/

  Getting Genetics Done

* http://omicsomics.blogspot.com/

  Omics! Omics!

* http://lab.loman.net/

  Nick Loman's lab notebook

Complete table of contents
==========================

.. toctree::
   :maxdepth: 5

   day1
   day2
   variant

   assembling-ecoli-with-velvet
   interval-analysis-and-visualization
   sam-format-tutorial
   R_Introductory_tutorial_2014
   IntroductionControlFlowR
   vcf_exploration
   eel-pond

   SOAPdeNovoTrans_count_eXpress
   analyzing_drosophila_htseq
   drosophila_rnaseq1
   drosophila_rnaseq_bwa_htseq
   genome-comparison-and-phylogeny
   git-koans
   howe-mgrast
   howe-ncbi
   kmers-etc
   long-read
   mount_chris_snapshot
   short-read-quality-evaluation-ctb

   amazon/index
   instructors-guide
   code-of-conduct
