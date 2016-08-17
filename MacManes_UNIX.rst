===============
UNIX Toolkit
===============

The explicit goal of this tutorial is to provide you with some tricks and tools and to make your unix life happier. Also, to possibly introduce you to samtools.

|

See https://github.com/stephenturner/oneliners for a bunch of examples.

|

**ALIASES**: These are 'shortcuts' for make doing stuff more efficient.

::

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

**.profile customization**

::

  nano ~/.profile

**DOWNLOAD SOME DATA**

::

  mkdir /home/ubuntu/data
  cd /home/ubuntu/data
  curl -LO https://www.dropbox.com/s/5fymuyb1f2l8kfj/ngsfile.tar.gz
  tar -zxf ngsfile.tar.gz


**Find a file**

::

  cd $HOME
  find / -name Trinity.fasta
  find / -type f -size +10M 2> /dev/null

**sed**

::

  cd /home/ubuntu/data/ngs2015
  less Trinity.fasta
  sed 's_|_-_g' Trinity.fasta | grep ^'>' | head
  sed -i 's_|_-_g' Trinity.fasta

`A helpful guide to sed <http://www.grymoire.com/Unix/Sed.html>`_

**Number fasta def line**

::

  awk '/^>/{print ">" ++i; next}{print}' < Trinity.fasta > Trinity.numbered.fasta


**awk**

::

  less Trinity.counts.RNAseq.txt
  awk '{print $1 "\t" $3}' Trinity.counts.RNAseq.txt | head
  awk '$1 == "c996_g1_i1"' Trinity.counts.RNAseq.txt

**Random Stuff**

::

  cd -
  tmux/screen



========================
TERMINATE YOUR INSTANCE
========================
