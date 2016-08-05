================================================ 
Self Guided Transcriptome Assembly Exercise
================================================

I want you to develop a pipeline to do some stuff. Think about automation - e.g. what if I gave you not just 1 dataset, but 100. 

1. Quality control the sequences
2. Assembly using Trinity
3. Run Busco

--------------

Here are some data:

::

  https://www.dropbox.com/s/gah3q0gt3z708yu/rnaseq.1.fq.gz
  https://www.dropbox.com/s/jt6bapq2kpjnoox/rnaseq.2.fq.gz

--------------

**Launch an AMI and update it**. For this exercise, we will use a **c4.4xlarge** We need to mount a hard  drive of size 300Gb. There are a few different ways to do this, but I like to add the drive during the launching of the instance using the "add storage" tab. See ``http://angus.readthedocs.org/en/2015/MacManes_Trinity.html``

--------------

**Install other software**: Follow directions at ``http://angus.readthedocs.org/en/2015/MacManes_Trinity.html``

----------------

**Quality Trim**: See ``http://angus.readthedocs.org/en/2015/MacManesTrimming.html``. What Quality threshold will you use?

----------------

**Assemble using Trinity**: See ``http://angus.readthedocs.org/en/2015/MacManes_Trinity.html`` 

----------------

**Run Busco**: You will need to download the Arthropod reference from the BUSCO website. See ``http://angus.readthedocs.org/en/2015/MacManes_Trinity.html`` 


**Optional Extra Credit Exercise** Do some annotation stuff like Meg showed you. 

================================================ 
Terminate your instance
================================================
