===============================================
Using GitHub repositories to store your scripts
===============================================

First, create a new repository
==============================

Go to github.com, sign in, and select "New repository".

.. image:: ./_static/github-1.png
   :width: 80%

Second, name your repository.
=============================

Something like 'rnaseq-scripts' is fine!  Be sure to check the
'initialize' box.

.. image:: ./_static/github-2.png
   :width: 80%

Third, clone your repository
============================

Get your GitHub repository url (https://github.com/ctb/rnaseq-scripts.git is
mine, in this example);

.. image:: ./_static/github-3.png
   :width: 80%

Then, on your remote UNIX machine, do::

   git clone https://github.com/ctb/rnaseq-scripts.git

This will create a new directory named ``rnaseq-scripts`` with a single
`'README.md`` file in it.

Fourth, add, commit, and push scripts
=====================================

The following commands are useful:

* ``git add script.txt`` will add the file script.txt into your local git
  repository.

* ``git commit -am "some message"`` will save the latest version of the script
  into your local git repository

* ``git push origin`` will send your committed versions to github, where
  they will be safe.

Other commands:

* ``git log`` will show you your history.

* ``git pull`` will update your command-line repository from your
  GitHub account.
