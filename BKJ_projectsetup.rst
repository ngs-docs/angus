.. _BKJ_projectsetup:

Setting up your RNA-seq project
===============================

So you've decided to do an RNA-seq experiment. Some things you will need to think about are:

	* What is the end-goal? (differential gene expression, variant calling, *de novo* assembly)
	
	* How many replicates do I need per condition? (related to first question; 3 is good, 6 is better)
	
	* How am I going to analyze this data?

You will save yourself a *lot* of headache if you take the time to set up your project in
a coherent manner. 

Some tips:

#. Put your project in a single directory/folder (named something informative without spaces)

#. :ref:`notebooks` for your analysis, commands, parameters, etc. - plain text or markdown

	- I like markdown because you can readily export to HTML to share with collaborators

#. Keep your raw data *raw* (e.g. make a working copy)

#. Name files something concise (NO SPACES!) but informative (tab auto-complete is your friend)

#. Make several *separate* folders within the main project folder to house different steps of the analysis

#. In each subfolder create a README to describe the contents of that folder, how they got there, etc. (metadata)

#. If possible, automate (script) analysis to allow you to re-run analyses on (many) files quickly and reproducibly

To put it another way, treat it as you would setting up a bench experiment and documenting
progress in a lab notebook. You will write down dates, times, reagents (lot numbers, purity, etc.),
buffer/media recipes, purpose of the experiment, experimental procedure, results (and if you believe them),
and what you plan to do next to continue to investigate the question or if there are new questions
to be answered.

You want to approach bioinformatics analyses in a similar manner. Document everything as best
as you can. It will be useful for when you share data with others (publish, collaborate, etc.)
and when you have to go back in 6 months to reassess what it was you did to generate graphs
x, y, and z.

Also, be aware that tools and technology evolve quickly in this field. It's worth getting a Twitter account to interact
with people who develop/perform these analyses and probing the literature for new/improved
tools (they come out *frequently*). Further, Google is your friend and there are all kinds of resources (forums, etc.)
to ask questions and get answers (usually pretty quickly).

.. _notebooks:

Start a virtual notebook
------------------------

Let's have a look at the Markdown syntax and how you can use it to document your work.

Markdown is a plain text format that can be rendered to HTML and is quite nice if working collaboratively,
like on GitHub.

I will use an example from Vince Buffalo's book "Bioinformatics Data Skills" (*highly recommended*)::

	# *Zea Mays* SNP Calling

	We sequenced three lines of *zea mays*, using paired-end
	sequencing. This sequencing was done by our sequencing core and we
	received the data on 2013-05-10. Each variety should have **two**
	sequences files, with suffixes `_R1.fastq` and `_R2.fastq`, indicating
	which member of the pair it is.

	## Sequencing Files

	All raw FASTQ sequences are in `data/seqs/`:

    	$ find data/seqs -name "*.fastq"
    	data/seqs/zmaysA_R1.fastq
    	data/seqs/zmaysA_R2.fastq
    	data/seqs/zmaysB_R1.fastq
    	data/seqs/zmaysB_R2.fastq
    	data/seqs/zmaysC_R1.fastq
    	data/seqs/zmaysC_R2.fastq

	## Quality Control Steps

	After the sequencing data was received, our first stage of analysis
	was to ensure the sequences were high quality. We ran each of the
	three lines' two paired-end FASTQ files through a quality diagnostic
	and control pipeline. Our planned pipeline is:

	1. Create base quality diagnostic graphs.
	2. Check reads for adapter sequences.
	3. Trim adapter sequences.
	4. Trim poor quality bases.

	Recommended trimming programs:

 	- Trimmomatic
 	- Scythe
 	
When this is rendered, it looks like `this <https://github.com/vsbuffalo/bds-files/blob/master/chapter-02-bioinformatics-projects/notebook.md>`__

Here are some Markdown editors for

	* Mac - `<http://25.io/mou/>`__
	
	* Windows - `<http://markdownpad.com/>`__
	
If you use GitHub to work collaboratively on things, you can copy and paste things right
into a `Gist <http://gist.github.com>`__ and it will render for you (if you add the .md file extension when you name your file) and you can share privately or
publicly!

`Markdown syntax guide <https://en.support.wordpress.com/markdown-quick-reference/>`__

