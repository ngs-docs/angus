==========================================
Next-Gen Sequence Analysis Workshop (2015)
==========================================

.. Warning These documents are not maintained and their instructions may be
	out of date. However the GED Lab does maintain the `khmer protocols
	<http://khmer-protocols.readthedocs.org/>`__ which may cover similar
	topics. See also the `installation instructions for the current version
	of the khmer project <https://khmer.readthedocs.org/en/latest/install.html>`__. 

This is the schedule for the `2015 MSU NGS course <http://bioinformatics.msu.edu/ngs-summer-course-2015>`__.


.. , which ran from August 10th to August 21st, 2015.  If you're interested in this course in 2015, please see `the 2016 announcement <http://bioinformatics.msu.edu/ngs-summer-course-2014>`__.

This workshop has a :doc:`code-of-conduct`.

`Download all of these materials <https://github.com/ngs-docs/angus/archive/2015.zip>`__ or `visit the GitHub repository <https://github.com/ngs-docs/angus/tree/2015>`__.

`Meal schedule <_static/NGS_meal_schedule.pdf>`__ 


===============  =============================================================
Day              Schedule
===============  =============================================================
Monday 8/10      * 1:30pm lecture: `Welcome! <_static/2014-lecture1-welcome.pptx.pdf>`__ (Titus)
                 * Tutorial: :doc:`day1`
                 * 7pm: research presentations
                 
Tuesday 8/11     * Tutorial:doc:`day2`
                 * 9:15am lecture: `Sequencing considerations <_static/2014-lecture2-sequencing.pptx.pdf>`__ (Titus)
                 * 10:30am: tutorial, UNIX command line (Amanda?)
                 * 1:15pm: tutorial, :doc:`short-read-quality-evaluation-ctb` (Titus?)
                 * Evening: *firepit social*
                 
Wed 8/12         * 9:15am lecture: `Mapping and Assembly <_static/2014-lecture3-mapping-assembly.pptx.pdf>`__ (Titus)
                 * 10:30am: tutorial, :doc:`variant` (Titus)
                 * 1:15pm: tutorial, :doc:`running-command-line-blast` (Adina?)
                 * 7:15pm: lecture: More Amazon (Titus)

Thursday 8/13    * 9:15am lecture: `Genome assembly and analysis <_static/NGS_Acey_2013.06.20.01.pdf>`__ (Erich)
                 * 10:15am lecture: exploring short read data sets with k-mer analyses (Titus)
                 * 1:15pm: tutorial, :doc:`assembling-ecoli-with-velvet` (Titus)
                 * 5:30pm: *leave for Kalamazoo* `Bells <http://bellsbeer.com/eccentric-cafe/menu>`__ (No dinner at McCrary)

Friday 8/14      * 9:15am-noon `lecture <_static/2014-lecture-mrnaseq-protocol.pptx.pdf>`__/tutorial, :doc:`eel-pond` (Titus)
                 * 1:15pm: tutorial, `Advanced UNIX <https://github.com/datacarpentry/shell-genomics/>`__ (Amanda?)

Saturday 8/15    * 9:15am-3pm `lecture <_static/howe_mgrast.pptx>`__ / tutorial :doc:`howe-ncbi`  (Adina)
                 * *Lunch at McCrary 12pm - 1pm*
                 * *Takeout Dinner on the island*

Sunday 8/16      * Free Day
                 * *"Brunch" at McCrary 12pm - 1pm*
                 * *BBQ Dinner on the island*
            
Monday 8/17      * 9:15am lecture, `mRNAseq and counting I <_static/NGS2014_RNAseq_ID_1.pdf>`__ (Ian)
                 * 10:30am tutorial, :doc:`drosophila_rnaseq1` (Chris)
                 * 11:30am tutorial, :doc:`mount_chris_snapshot` (Meg)
                 * 2:00pm  :doc:`git-koans` (Adina)                
                 * 8:15pm: *firepit and gin social*
                                  
Tuesday 8/18     * 9:15am lecture, `mRNAseq and counting II <_static/NGS2014_RNAseq_2.pdf>`__ (Ian)
                 * 11:00am tutorial, :doc:`analyzing_drosophila_htseq` (Chris?)
                 * 1:15pm tutorial, :doc:`drosophila_rnaseq_bwa_htseq` (Meg)
                 * 2:15pm tutorial, :doc:`SOAPdeNovoTrans_count_eXpress` (Matt)
                 * 7:15pm tutorial, Mapping reads to transcriptomes and counting: Trinity and SOAP (Matt) 

Wed 8/19         * 9:15am-noon: lecture/tutorial, :doc:`R_Introductory_tutorial_2014` R etc (Ian)
                 * 1:15pm: lecture/discussion, `mRNAseq assembly with Trinity <_static/MegStaton_NGS_KBS_Staton_RNASeq.pdf>`__ (Meg)
                 * 3pm lecture, `A tableside discussion on transcriptome assembly <_static/ngs2014-trimming.pdf>`__ (Matt)
                 
Thursday 8/20    * 9:15-10:15
                 * 10:30am-noon 
                 * 1:15pm `lecture <_static/2014-kmers.pptx.pdf>`__ / tutorial :doc:`kmers-etc` (Titus)
                 * 3pm long reads `lecture <_static/Long_reads.pdf>`__ / tutorial :doc:`long-read` (Matt)
                 * *BBQ Dinner on the island*

Friday 8/21      * 9:15-9:45 `closing lecture <_static/2014-final-lecture.pptx.pdf>`__ (Titus)
                 * 10am discussion about class; more stuff

                 * Links:
                   `Opinionated guides to NGS <http://davis-assembly-masterclass-2013.readthedocs.org/en/latest/outputs/index.html>`__ /
                   `Software Carpentry <http://software-carpentry.org>`__

===============  =============================================================

Dramatis personae
=================

Instructors:

* C Titus Brown
* Ian Dworkin

TAs:

* Amanda Charbonneau
* Lisa Cohen
* Ryan Williams

Lecturers:

* Chris Chandler
* Adina Chuang Howe
* Matt MacManes
* Meg Staton
* Erich Schawrz
* Nick Loman
* Torsten Seemann

She Who Drives Many Places:

* Jessica Mizzi

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

