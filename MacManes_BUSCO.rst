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
  PATH=$PATH:$(pwd)/bin
  export AUGUSTUS_CONFIG_PATH=$(pwd)/config

**INSTALL BUSCO**: Needs ``apt-get install ncbi-blast+ hmmer``

::

  cd $HOME
  curl -O http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
  tar -zxf BUSCO_v1.1b1.tar.gz
  cd BUSCO_v1.1b1/
  PATH=$PATH:$(pwd)

**Go to assembly directory and download BUSCO reference**: BUSCO reference needs to be in same DIR as where you are going to run BUSCO.

::

  cd /mnt/assembly/ #or wherever you're assemblies are
  curl -LO https://www.dropbox.com/s/o8roa5ayt5dbnl4/bacteria_busco.tar.gz
  tar -zxf bacteria_busco.tar.gz

**Run BUSCO for assemblies**: There are Eukaryote, Metazoa, Arthropod, Vertebrate, Plant referneces for use with other genomes. 

::

  cd /mnt/reads/
  python3 /home/ubuntu/BUSCO_v1.1b1/BUSCO_v1.1b1.py \
  -m genome -in [GENOME.ASSEMBLY.FASTA] \
  --cpu 8 -l bacteria -o [NAME]


**Results for SPAdes vs. Velvet vs. Reference**

::

  Ref:    C:97%[D:0.0%],F:2.5%,M:0.0%,n:40
  SPAdes: C:97%[D:2.5%],F:2.5%,M:0.0%,n:40
  Velvet: C:75%[D:2.5%],F:0.0%,M:25%,n:40
