Transrate
###############

`Transrate <http://hibberdlab.com/transrate/index.html>` is a tool for assessing the quality of a de
novo transcriptome assembly.

We will use the assembly produced from the earlier Trinity tutorial.

Overview
~~~~~~~~~~~~~~

Transrate analyses a transcriptome assembly in three key ways:

- by inspecting the contig sequences
- by mapping reads to the contigs and inspecting the alignments
- by aligning the contigs against proteins or transcripts from a related species and inspecting the alignments

We'll try out the first two levels.
 
Using Trimmed Reads
-------------------

Lets map the reads and see how that looks. Transrate will
- map the provided reads to the assembly using `SNAP <http://snap.cs.berkeley.edu/>`_
- infer the most likely contig of origin for any multi-mapping reads with `Salmon <https://github.com/kingsfordgroup/sailfish/releases/tag/v0.3.0>`_
- inspect the resulting alignments with transrate-tools and use them to evaluate each contig in the assembly

Run transrate, pointing to the assembly and the input fastq files. 

 /root/transrate/transrate-1.0.1-linux-x86_64/transrate \
 --assembly /mnt/assembly/trinity_out_dir/Trinity.fasta \
 --threads 16 \
 --left /mnt/trimming/subsamp.Phred30_1P.fq \
 --right /mnt/trimming/subsamp.Phred30_2P.fq
 

Lots of output files. Lets focus on the main metrics. `Explanations are here <http://hibberdlab.com/transrate/metrics.html#read-mapping-metrics>`_

The main output file is transrate_assemblies.csv (CSV = comma separated value)

Transfer via dropbox or scp and look at it in Excel.

Mac Example with scp::

scp -i ~/Desktop/files/teaching/kbs-ngs-2014/2015/amazon.pem ubuntu@ec2-??-??-??-??.compute-1.amazonaws.com:/mnt/transrate/transrate_results/assemblies.csv .
scp -i ~/Desktop/files/teaching/kbs-ngs-2014/2015/amazon.pem ubuntu@ec2-??-??-??-??.compute-1.amazonaws.com:/mnt/transrate/transrate_results/Trinity/contigs.csv .

Lets see what each of these mean by looking at the `documentation page <http://hibberdlab.com/transrate/metrics.html#contig-metrics>`_.




