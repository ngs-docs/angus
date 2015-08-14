===============
BUSCO
===============

http://busco.ezlab.org/

**INSTALL Augustus**

::

  cd $HOME
  curl -O http://augustus.gobics.de/binaries/augustus.2.5.5.tar.gz
  tar -zxf augustus.2.5.5.tar.gz
  cd augustus.2.5.5/
  make
  export PATH=$PATH:$(pwd)/bin
  export AUGUSTUS_CONFIG_PATH=$(pwd)/config

**INSTALL BUSCO**:

Some new system software is required::

   sudo -y apt-get install ncbi-blast+ hmmer

Now, install BUSCO::

  cd $HOME
  curl -O http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
  tar -zxf BUSCO_v1.1b1.tar.gz
  cd BUSCO_v1.1b1/
  export PATH=$PATH:$(pwd)

**Go to assembly directory and download BUSCO reference**: BUSCO reference needs to be in same DIR as where you are going to run BUSCO.

::

  cd /mnt/assembly/ #or wherever your assemblies are
  curl -LO https://www.dropbox.com/s/o8roa5ayt5dbnl4/bacteria_busco.tar.gz
  tar -zxf bacteria_busco.tar.gz

**Run BUSCO for assemblies**: There are Eukaryote, Metazoa, Arthropod, Vertebrate, Plant referneces for use with other genomes. 

::

  python3 /home/ubuntu/BUSCO_v1.1b1/BUSCO_v1.1b1.py -f \
     -m genome -in megahit-assembly.fa \
     --cpu 8 -l bacteria -o megahit
  python3 /home/ubuntu/BUSCO_v1.1b1/BUSCO_v1.1b1.py -f \
     -m genome -in spades-assembly.fa \
     --cpu 8 -l bacteria -o spades -f

**Results for SPAdes vs. Velvet vs. Reference**

Run::

  cat run*/short*
