Looking at k-mer abundance distributions
========================================

Start up an Ubuntu 16.10 instance and run ::

   sudo bash

   apt-get update
   apt-get -y install screen git curl python3.5-dev python3.5-venv make \
       libc6-dev g++ zlib1g-dev
   python3.5 -m venv ~/py3 && \
   . ~/py3/bin/activate && pip3 install -U pip &&  \
   pip3 install -U jupyter jupyter_client ipython pandas matplotlib scipy scikit-learn

Install `khmer <http://khmer.readthedocs.org/en/v2.1.1/>`__::

   cd /usr/local/share
   git clone https://github.com/ged-lab/khmer.git
   cd khmer
   git checkout v2.1.1
   make install
   exit                 # and then leave the "sudo bash" environment

Create our workspace ::

   sudo mkdir -p /mnt
   sudo chown $USER /mnt
   echo '. ~/py3/bin/activate' >> ~/.bash_profile
   source ~/.bash_profile

Download the notebook for today's session::

   cd /mnt
   curl -O https://raw.githubusercontent.com/ngs-docs/angus/2017/files/kmer-abundance-distributions.ipynb

Setup the `Jupyter Notebook <http://ipython.org/notebook.html>`__::

    jupyter notebook --generate-config
    jupyter notebook passwd

Here, you will choose a password to authenticate to your notebook server.
Can anyone think about why you might not want your instance
configured to run arbitrary (shell) commands from anywhere?

And finally start the notebook server listening:
    echo Try to connect to http://$(hostname -I)
    jupyter notebook  --no-browser --ip='*' --port=8080 --notebook-dir=/mnt

If you are doing this on Amazon EC2, you need to enable inbound traffic
on HTTP to your computer before you can connect to the notebook server.
Briefly, go to your instance and look at
what security group you're using (should be 'launch-wizard-'
something).  On the left panel, under Network and Security, go into
Security Groups. Select your security group, and select Inbound, and
Edit.  Click "Add rule", and change "Custom TCP rule" to "http".  Then
click "save".  Done!

Now, open your machine's address in your Web browser.

---

