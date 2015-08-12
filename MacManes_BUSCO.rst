===============
BUSCO
===============

**INSTALL Augustus**

::

  cd $HOME
  curl -O http://augustus.gobics.de/binaries/augustus.2.5.5.tar.gz
  tar -zxf augustus.2.5.5.tar.gz
  cd augustus.2.5.5/
  make
  PATH=$PATH:$(pwd)/bin
  export AUGUSTUS_CONFIG_PATH=$(pwd)/config

**INSTALL BUSCO**: Needs ncbi-blast+ and hmmer

::

  cd $HOME
  curl -O http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
  tar -zxf BUSCO_v1.1b1.tar.gz
  cd BUSCO_v1.1b1/
  PATH=$PATH:$(pwd)

**Go to assembly directory and download BUSCO reference**

::

  cd /mnt/reads/
  curl -O http://busco.ezlab.org/files/bacteria_buscos.tar.gz
  tar -zxf bacteria_buscos.tar.gz

**Run BUSCO for assemblies**

::

  cd /mnt/reads/
  python3 /home/ubuntu/BUSCO_v1.1b1/BUSCO_v1.1b1.py \
  -m genome -in Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
  --cpu 8 -l eukaryota -o cel2

  
