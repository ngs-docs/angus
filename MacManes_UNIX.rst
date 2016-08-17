===============
UNIX Toolkit
===============

The explicit goal of this tutorial is to provide you with some tricks and tools and to make your unix life happier.

|

See https://github.com/stephenturner/oneliners for a bunch of examples.

|

**Profile customization and aliases**: These are 'shortcuts' for make doing stuff more efficient.

::

    nano ~/.profile

    alias c='clear'
    alias gh='history | grep'
    alias ll='logout'
    alias l='ls -lth'
    alias mv='mv -i'
    alias cp='cp -i'


    alias tn="tmux new -s"
    alias ta="tmux attach -t"
    alias tl="tmux ls"
    alias tk="tmux kill-session -t"



**DOWNLOAD SOME DATA**

::

  mkdir /home/ubuntu/data && cd /home/ubuntu/data
  curl -LO https://s3.amazonaws.com/reference_assemblies/Oophaga/transcriptome/Oophaga_pumilio.1.0.0.fasta
  curl -L https://github.com/ngs-docs/angus/blob/2016/_static/quantification.tgz?raw=true > quantification.tgz
  tar -zxf *gz


**Find a file**

::

  cd $HOME
  find / -name quant.sf
  find / -type f -size +10M 2> /dev/null

**sed**

::


  less Oophaga_pumilio.1.0.0.fasta
  sed 's/_/+++bananas+++/g' Oophaga_pumilio.1.0.0.fasta | grep ^'>' | head
  sed -i 's/_/+++bananas+++/g' Oophaga_pumilio.1.0.0.fasta

`A helpful guide to sed <http://www.grymoire.com/Unix/Sed.html>`_

**Number fasta def line**

::

  awk '/^>/{print ">" ++i; next}{print}' < Oophaga_pumilio.1.0.0.fasta > renamed.fasta


**awk**

::

  file=$(find . -name quant.sf | head -1)
  awk '{print $1 "\t" $3}' $file | head
  awk '$1 == "FBtr0070078"' $file

**Random Stuff**

::

  cd -



========================
TERMINATE YOUR INSTANCE
========================
