Understanding the SAM format
============================

Log into your instance, create a new directory, navigate to that directory::

    cd /mnt
	mkdir sam
	cd sam

	# Get the makefile.
	wget https://raw.githubusercontent.com/ngs-docs/angus/2014/files/Makefile-samtools -O Makefile

A series of exercises will show what the SAM format is and how it changes when
the query sequence is altered and how that reflects in the output.

Also, for the speed of result generation here is a one liner to generate a bamfile::

	# One line bamfile generation.
	bwa mem index/sc.fa query.fa | samtools view -bS - | samtools sort - results

This will produce the ``results.bam`` output.



