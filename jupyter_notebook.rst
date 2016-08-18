======================
Jupyter notebook demo
======================

`Jupyter <http://jupyter.org/>`__ notebook is an interactive web application that allows you to type and edit lines of code and see the output. The software requires Python installation, but currently supports interaction with over 40 languages. 

`Jupyter notebook examples <https://github.com/ipython/ipython/wiki/A-gallery-of-interesting-IPython-Notebooks>`__ 

Installation:
=============

`Jupyter install instructions <http://jupyter.readthedocs.io/en/latest/install.html>`__ 

First we will install Anaconda, a package manager for Python libraries.

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
  
If you are running windows or have problems on Mac, try these `remote instructions <https://github.com/WhiteheadLab/Computational_Protocols/blob/master/install_jupyter_notebook_farm.md>`__. Start an AWS (does not have to be big, try Ubuntu 14.04 t2.micro). Login and install packages:
::
  sudo apt-get update && \
  sudo apt-get -y install build-essential git curl gcc make g++ python-dev unzip \
    default-jre pkg-config libncurses5-dev r-base-core
  curl -OL https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-2.4.0-Linux-x86_64.sh
  bash Anaconda2-2.4.0-Linux-x86_64.sh
  source ~/.bashrc
  conda install jupyter
  conda install -c r r r-essentials
  source ~/.bashrc
  
Clone the repo:
::
  git clone https://github.com/ljcohen/ngs2016_important_files.git

Start a notebook by specifying a port:
::
  cd ngs2016_important_files
  jupyter notebook --port 8756

The notebook software will start running. Then open a new terminal window and type:
::
  ssh -i ~/xxx.pem -L 8888:localhost:8756 ubuntu@xxx.amazon.com

Just keep this open, you don't need to type anything in this terminal window. Then, open up your browser and type this into the url:
::
  http://localhost:8888/tree/

You should see the files in the repository. Click on ``1.01_StringBasics.ipynb`` and you are running a jupyter notebook!

Using Jupyter notebooks:
========================

The main keyboard command to remember is how to execute the code from a cell. Type code into a cell and then hit ``Shift-Enter``.

If you're in Python 2, type:
::
  print "Hello World!"

or for Python 3:
::
  print("Hello World!")

Then press ``Shift-Enter``

For more instructions, the Help menu has a good tour and detailed information. Notebooks can be downloaded locally by going to the File menu, then selecting Download and choosing a file type to download.

References for learning Python
=============================
* http://rosalind.info/problems/locations/ 
* http://learnpythonthehardway.org/book/ 
* http://www.learnpython.org/
* http://www.pythontutor.com/visualize.html#mode=edit
