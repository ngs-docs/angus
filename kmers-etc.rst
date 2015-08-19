Looking at k-mer abundance distributions
========================================

Start up an Ubuntu 14.04 instance and run ::

   sudo bash

   apt-get update
   apt-get -y install screen git curl gcc make g++ python-dev unzip \
           default-jre pkg-config libncurses5-dev r-base-core \
           r-cran-gplots python-matplotlib sysstat ipython-notebook

Install `khmer <http://khmer.readthedocs.org/en/v1.1/>`__::

   cd /usr/local/share
   git clone https://github.com/ged-lab/khmer.git
   cd khmer
   git checkout v1.1
   make install

Download the notebook for today's session::

   cd /mnt
   curl -O https://raw.githubusercontent.com/ngs-docs/angus/2014/files/kmer-abundance-distributions.ipynb

Run `IPython Notebook <http://ipython.org/notebook.html>`__::

   ipython notebook --no-stdout --no-browser --ip='*' --port=80 --notebook-dir=/mnt

(Note: you can use:

  --directory=/root/Dropbox

instead if you want to store notebooks in your Dropbox.)

To connect to your notebook server, you need to enable inbound traffic
on HTTP to your computer.  Briefly, go to your instance and look at
what security group you're using (should be 'launch-wizard-'
something).  On the left panel, under Network and Security, go into
Security Groups. Select your security group, and select Inbound, and
Edit.  Click "Add rule", and change "Custom TCP rule" to "http".  Then
click "save".  Done!

Now, open your EC2 machine's address in your Web browser.

---

