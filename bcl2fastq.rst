bcl2fastq demultiplexing
================================

http://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html

For each sequencing run, Illumina sequencers create base call files (.bcl) from each image on each tile for each cycle. These bcl files need to be stiched together with phred quality probability scores to create .fastq files. Most sequencing facilities provide .fastq files as the output of their services. If your facility does not do this, or if you're interested in how this process works, here are instructions for how to work with the software.

The conversion software used to be called **Cassava** until version 1.8.2. From version 1.8.4, it is now called **bcl2fastq**. Depending on the type of sequencer you have used, you will choose the version of bcl2fastq:

Note: Illumina instruments create a complex, albeit consistent directory structure of files that are all required for the bcl2fastq conversion process. Below is the directory structure for a run where YYMMDD is date, machine name is the ID for the instrument, e.g. D00595, XXXXX is the name of the project or PI or identifying information from the sequencing facility, and FC is the ID for the disposable flow cell used for that run, e.g. BC6UHTANXX. You need all of the directories and files output from the sequencing instrument. They will look like this (borrowed from `genomics-bcftbx <http://genomics-bcftbx.readthedocs.org/en/latest/protocols/prep_illumina.html>`_):

.. code::

        <YYMMDD>_<machinename>_<XXXXX>_FC/
             |
             +-- Data/
             |     |
             |     +------ Intensities/
             |                  |
             +                  +-- .pos files
             |                  |
             |                  +-- config.xml
             +-- RunInfo.xml    |
                                +-- L001(2,3...)/  (lanes)
                                |
                                +-- BaseCalls/
                                       |
                                       +-- config.xml
                                       |
                                       +-- SampleSheet.csv
                                       |
                                       +--L001(2,3...)/  (lanes)
                                             |
                                             +-- C1.1/   (lane and cycle)
                                                   |
                                                   +-- .bcl(.gz) files
                                                  |
                                                  +-- .stats files



According to the `software instructions <http://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html>`_ : 

    For Illumina sequencing systems running RTA version 1.18.54 and later, use bcl2fastq2 Conversion Software v2.17.
    For Illumina sequencing systems runnings RTA versions earlier than 1.18.54, use bcl2fastq Conversion Software v1.8.4.

`Real Time Analysis (RTA) analysis software <https://support.illumina.com/sequencing/sequencing_software/real-time_analysis_rta.html>`_ versions differ between Illumina instrument models. HiSeq and MiSeq instruments can usually run on v1.8.4. NextSeq is usually on v2.17.

We will install and use v1.8.4 here. The user manual for bcl2fastq v1.8.4 is here:
http://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf

Install bcl2fastq 
=================

Go to AWS and start a Red Hat Linux m4.large instance, with 500 GB volume attached and login:

.. code::

    ssh -i blah-blah.pem ec2-user@amazon.blah

.. code::

Mount storage, change permissions, and software:

.. code::

    sudo mkfs -t ext4 /dev/xvdb
    sudo chown -R ec2-user:ec2-user /mnt
    df -h
    sudo yum install wget screen dos2unix mac2unix

Make a directory for your bcl2fastq program files:

.. code::

    cd
    mkdir bcl2fastq
    cd bcl2fastq

Download the bcl2fastq program rpm file:

.. code::

    curl -O ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/Software/bcl2fastq/bcl2fastq-1.8.4-Linux-x86_64.rpm

Install the rpm file:

.. code::

    sudo yum install bcl2fastq-1.8.4-Linux-x86_64.rpm

Test to see if this works, should output instructions for how to run:

.. code::

    configureBclToFastq.pl -h | less

Now we need to get the sequencing files. If your sequencing files are in your Dropbox, connect it to your AWS instance to access the files:

    * http://ged.msu.edu/angus/tutorials-2011/installing-dropbox.html
    * https://www.dropbox.com/en/install?os=lnx

.. code::

    cd ~ && wget -O - "https://www.dropbox.com/download?plat=lnx.x86_64" | tar xzf -
    mkdir /mnt/Dropbox/
    cd /mnt/Dropbox/
    curl -O https://linux.dropbox.com/packages/dropbox.py
    python dropbox.py start -i
    python dropbox.py start
    python dropbox.py status

    ~/.dropbox-dist/dropboxd

..and you should see a message like this:

    >    This client is not linked to any account... Please visit https://www.dropbox.com/cli_link?host_id=XXXXX to link this > machine.

Copy/paste that URL into your Web browser; log into dropbox; and voila! The directory ~/Dropbox will be linked into your home directory!

    >    This computer is now linked to Dropbox. Welcome __!!

(NOTE: This might take a while if your Dropbox has a lot of files in it. It is easier to create a new Dropbox account with only these files.)

Copy the files from the Dropbox directory to your /mnt directory:

.. code::
                cp /Dropbox/<sequencing run files>  /mnt/run_files

If your files are coming from a high performance computing server, you can use scp to securely transfer the files:

.. code::

             cd /mnt/run_files
             screen
             scp -rp ligh@143.107.29.100:/home/terra/ngs/dados/corridasHiSeq/150804_D00549_0041_Bh3ntvbcxx/Data/Intensities/BaseCalls/*.*  .

Type Ctrl-A-D to detach from screen. Take a break. This will take several hours to copy ~5 GB of files. 

Configure SampleSheet.csv
=========================

A file SampleSheet.csv is required for the conversion program. It will contain your barcodes and sample ID information. It must be in a specific format with 10 column headers: "FCID", "Lane", "SampleID", "SampleRef", "Index", "Description", "Control", "Recipe", "Operator", "SampleProject". Not all of these fields are required. The Flow Cell ID (FCID), e.g. BC6UHTANXX must match the same FCID as the run. The "Index" column contains the barcode sequences. Dual index barcodes are separated by a "-" character. No spaces or special characters should be used in the sample sheet. Do not use these characters: $%^&*()!@~"';:?/}{

An example SampleSheet.csv is here:

https://dl.dropboxusercontent.com/u/9205689/SampleSheet.csv

Run bcl2fastq
=============

Run these commands:

.. code::

    dos2unix SampleSheet.csv
    mac2unix SampleSheet.csv
    OUT_DIR="/mnt/run_files/Unaligned/"
    IN_DIR="/mnt/run_files/"
    configureBclToFastq.pl \
    --input-dir $IN_DIR \
    --output-dir $OUT_DIR \
    --fastq-cluster-count 0 \
    --mismatches 1
    

If you don't have all the appropriate files, you will see an error message similar to this:

.. code::

        [2015-08-17 21:18:28]	[configureBclToFastq.pl]	INFO: Creating directory '/mnt/demultiplexing/Unaligned'
            ERROR: /mnt/demultiplexing/config.xml: file does not exist
            at /usr/local/lib/bcl2fastq-1.8.4/perl/Casava/Demultiplex.pm line 116.


Other configurations
====================

If you have different length barcodes or need to modify your SampleSheet.csv, here are some additional configurations for bcl2fastq. 

If things go bad (indecipherable errors), try adding one or all of these flags to the configuration above if: 

.. code::

        --ignore-missing-control --ignore-missing-stats --ignore-missing-bcl \

Instead of demultiplexing with barcodes, if you want to generate an index read containing all barcodes (if you have dual index barcodes, nextera):

.. code::

        mv -v ${BASE_CALLS_DIR}/SampleSheet.csv ${BASE_CALLS_DIR}/SampleSheet.0.csv
        /local/apps/bcl2fastq/1.8.4/bin/configureBclToFastq.pl \
        --input-dir ${BASE_CALLS_DIR} \
        --output-dir ${BASE_CALLS_DIR}/Unaligned \
        --fastq-cluster-count 0 \
        --use-bases-mask y*,y*,y*,y*

If you have single index, replace last line of above with this:

.. code::

        --use-bases-mask y*,y*,y*

Dual 8bp index read (nextera)

.. code::

        /local/apps/bcl2fastq/1.8.4/bin/configureBclToFastq.pl \
        --input-dir ${BASE_CALLS_DIR} \
        --output-dir ${BASE_CALLS_DIR}/Unaligned \
        --fastq-cluster-count 0 \
        --mismatches 0 \
        --use-bases-mask y*,i8,i8,y*


In-read barcodes

.. code::

        /local/apps/bcl2fastq/1.8.4/bin/configureBclToFastq.pl \
        --input-dir ${BASE_CALLS_DIR} \
        --output-dir ${BASE_CALLS_DIR}/Unaligned \
        --fastq-cluster-count 0 \
        --use-bases-mask i6y*,n*


More than one length barcode in same run

.. code::

        /local/apps/bcl2fastq/1.8.4/bin/configureBclToFastq.pl \
        --input-dir ${BASE_CALLS_DIR} \
        --output-dir ${BASE_CALLS_DIR}/Unaligned \
        --fastq-cluster-count 0 \
        --mismatches 0 \
        --use-bases-mask y*,i6n*,y*
        
Other references
================

* Many of these configurations are from Igor Dolgalev, the demultiplexing guru at GTC, NYUMC: igor.dolgalev@nyumc.org
* http://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html
* http://genomics-bcftbx.readthedocs.org/en/latest/protocols/prep_illumina.html  
* https://www.biostars.org/p/44927/


