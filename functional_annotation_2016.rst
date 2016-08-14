=====================
Functional Annotation
=====================

Attach a Volume
===============
Lets add a volume - as a reminder, this is sort of like plugging in an external harddrive to your laptop. Its adding a set of data, in this case, the databases that we are going to use for comparison for functional annotation.

This requires a few new steps in our tried and true Amazon EC2 instance protocol.

1. Choose AMI - Still going to choose Ubuntu server 14.04. Click Select.
2. Choose c4.2xlarge - Click Next:Configure Instance.
3. **New Step** For Subnet, select the one that ends in "us-east-1c". Click Next: Add Storage.
4. Change Size 8Gb to 100Gb.
5. **New Step** Click add Volume. 

   - Enter snapshot 'snap-6df6088a'
   - Change to 100Gb 
   - And change to /dev/sdf

6. Click review and launch.

Now you can SSH into your instance as normal.

Make your Volume available. 
===========================

We made the mount point /dev/xvdf.  ::

  sudo mkdir /mnt/dammit_dbs
  sudo mount /dev/xvdf /mnt/dammit_dbs/

Check out all the files you now have ::

  ls /mnt/dammit_dbs/

Installs
========

Now start the install process.  First install linux brew, as described `here <http://angus.readthedocs.io/en/2016/linuxbrew_install.html>`__.

The rest of the installs will take ~10 minutes. Start by getting all the dependencies available in brew ::

	brew install last blast emboss
	brew install busco
	brew install hmmer infernal

Next, there is a reciprocal best hit blast software that can be installed with the ruby gems package management system::

	sudo gem install crb-blast

Transdecoder from source::

	curl -LO https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz
	tar -xvzf 2.0.1.tar.gz
	cd TransDecoder-2.0.1; make
	export PATH=$HOME/TransDecoder-2.0.1:$PATH
	cd

Now we'll use pip, the python package management system to get the last few dependencies and dammit::

	pip install -U setuptools
	pip install numpy
	pip install dammit
	ln -s /home/linuxbrew/.linuxbrew/bin/busco /home/linuxbrew/.linuxbrew/bin/BUSCO_v1.1b1.py


Whew. Lets make sure it all worked.::

	dammit dependencies

If all is well, it will say "All dependencies satisfied!"

And InterProScan (the second part of the lesson) needs java::

	sudo apt-get install default-jre
	sudo apt-get install default-jdk
	export JAVA_HOME=/usr/lib/jvm/java-8-oracle
	export PATH=$PATH:/usr/lib/jvm/java-8-oracle

Dammit
======
Dammit is now installed but it doesn't have any databases to compare anything to. These are also stored on the volume. Next step, tell dammit how to use the databases on the volume::

	dammit databases --database-dir /mnt/dammit_dbs/ --full --busco-group eukaryota

I think I missed one and I'm too lazy to make a new snapshot, lets grab it::

    dammit databases --database-dir /mnt/dammit_dbs/ --full --busco-group eukaryota --install

If you didn't have access to the volume, this would download all the dbs and install them (30-45 minutes of time)

Get example data and unzip::

	wget ftp://ftp.ebi.ac.uk/pub/databases/pombase/FASTA/cdna_nointrons_utrs.fa.gz
	gunzip cdna_nointrons_utrs.fa
	grep -c '^>' cdna_nointrons_utrs.fa

Dammit takes quite a while on the whole set, so lets extract a smaller set to practice with::
	head -500 cdna_nointrons_utrs.fa > practice.fa

Run annotation::
	dammit annotate practice.fa --database-dir /mnt/dammit_dbs/ --busco-group eukaryota --n_threads 15

You will get results from 6 different analyses:
#. Pfam-A
#. Rfam
#. OrthoDB
#. BUSCO
#. Uniref90
#. Transdecoder

InterProScan
============

Dammit runs a lot of good software, but you may want to also assign GO terms. InterProScan is one way to do this. It takes a while to install, so I put a copy on the volume. Lets add it to our path::

	export PATH=$PATH:/mnt/dammit_dbs/interproscan-5.19-58.0

And see the documentation::

	interproscan.sh | less



Notes on installing Interproscan
================================
If the volume is no longer available and you need to install IPS, then here's how to do it. Java needs to be installed (see above), then download and unpack the most recent version.::
	wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.19-58.0/interproscan-5.19-58.0-64-bit.tar.gz
	tar xvzf interproscan-5.19-58.0-64-bit.tar.gz



