bcl2fastq demultiplexing
================================

http://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html

Illumina sequencers create base call files (.bcl) from each image on each tile for each cycle of a sequencing run. These bcl files need to be stiched together with phred quality probability scores to create .fastq files. Most sequencing facilities provide .fastq files as the output of their services. If your facility does not do this, or if you're interested in how this process works, here are instructions for how to work with the software.

The software used to be called **Cassava** until version 1.8.2. From version 1.8.4, it is now called **bcl2fastq**. Depending on the type of sequencer you have used, you will choose the version of bcl2fastq:

According to the `software instructions <http://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html>`_ : 

    For Illumina sequencing systems running RTA version 1.18.54 and later, use bcl2fastq2 Conversion Software v2.17.
    For Illumina sequencing systems runnings RTA versions earlier than 1.18.54, use bcl2fastq Conversion Software v1.8.4.

The user manual is here:
http://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf


Install bcl2fastq 
=================

Go to AWS and start a Red Hat Linux m4.large instance and login:

.. code::

    ssh -i ????.pem ec2-user@amazon.blah

.. code::

Since running the bcl2fastq program will take some time, you will want to start a screen so that if you become disconnected, your work will not be lost: 

.. code::

    sudo yum install wget screen

Make a directory for your bcl2fastq program files:

.. code::

    mkdir bcl2fastq
    cd bcl2fastq

Download the bcl2fastq program rpm file:

.. code::

    curl -O ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/Software/bcl2fastq/bcl2fastq-1.8.4-Linux-x86_64.rpm

Install the rpm file:

    sudo yum install bcl2fastq-1.8.4-Linux-x86_64.rpm

Test to see if this works, should output instructions for how to run:

    configureBclToFastq.pl -h | less


If your bcl files are in your Dropbox, you will need to connect it to your AWS instance so you can access the files:

http://ged.msu.edu/angus/tutorials-2011/installing-dropbox.html
https://www.dropbox.com/en/install?os=lnx

    cd ~ && wget -O - "https://www.dropbox.com/download?plat=lnx.x86_64" | tar xzf -

..and you should see a message like this:

    This client is not linked to any account... Please visit https://www.dropbox.com/cli_link?host_id=XXXXX to link this machine.

Copy/paste that URL into your Web browser; log into dropbox; and voila! The directory ~/Dropbox will be linked into your home directory!

    This computer is now linked to Dropbox. Welcome __!!

This might take a while if your Dropbox has a lot of files in it.


Run bcl2fastq
=============

Configure SampleSheet.csv

Run these commands:

.. code::

    dos2unix
    mac2unix
    OUT_DIR="/"
    








