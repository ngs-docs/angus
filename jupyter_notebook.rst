======================
Jupyter notebook demo
=====================

Jupyter notebook examples: https://github.com/ipython/ipython/wiki/A-gallery-of-interesting-IPython-Notebooks
Jupyter install instructions: http://jupyter.readthedocs.io/en/latest/install.html

Anaconda package manager for Python libraries is strongly recommended.

Install on AWS, or other UNIX machiens:
::
  curl -OL https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-2.4.0-Linux-x86_64.sh
  bash Anaconda2-2.4.0-Linux-x86_64.sh
  
Install on Mac OS X:
::
  curl -OL http://repo.continuum.io/archive/Anaconda2-4.1.1-MacOSX-x86_64.sh
  bash Anaconda2-4.1.1-MacOSX-x86_64.sh

Type ``Enter`` when prompted, then press Enter through the instructions. Be careful not to keep pressing Enter without reading otherwise you will end up saying ```No``. You will need to type 'Yes' to continue with the installation.

After Anaconda install has finished, type:
::
  source ~/.bashrc
  conda install jupyter
  conda install -c r r r-essentials
  
This will install packages allowing you to open either a new Python .ipynb or an R .ipynb. 

Navigate to the directory on your computer with files you want to explore. Then type:
::
  jupyter notebook

This will open your browser with a list of files. Click on the "New" and bring down the pull-down menu. Under 'Notebooks', click on either R or Python language to start your new notebook!

You can also open existing an .ipynb. Today we'll explore an `existing .ipynb by Ben Langmead <https://github.com/BenLangmead/ads1-notebooks/blob/master/1.01_StringBasics.ipynb>`__.

See here for additional tutorial materials:

http://dib-training.readthedocs.io/en/pub/2016-03-09-jupyter-notebook.html

mybinder: https://github.com/ngs-docs/2016-mar-jupyter
  
Remote instructions: https://github.com/WhiteheadLab/Computational_Protocols/blob/master/install_jupyter_notebook_farm.md

Python References
===========
* http://rosalind.info/problems/locations/ 
* http://learnpythonthehardway.org/book/ 
* http://www.learnpython.org/
* http://www.pythontutor.com/visualize.html#mode=edit
