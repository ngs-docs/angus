===============================================
Using GitHub repositories to store your scripts
===============================================

By the end of this lesson, you will be able to:

1. Make a new github repo
2. Edit README.md (`markdown cheat sheet <https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet>`__ (`Why README.md is necessary <https://changelog.com/top-ten-reasons-why-i-wont-use-your-open-source-project/>`__
3. Fork and clone a repo
4. Make changes, commit and push changes
5. Submit a pull request
6. Submit an issue
7. Use Jupyter notebook to explore data


Create a new repository
==============================

Go to `github.com <https://github.com/>`__, sign in (sign up if you haven't already), go to your profile page (e.g. `https://github.com/ljcohen <https://github.com/ljcohen>`__) then select "New repository".

.. image:: ./_static/git_create_repo.png
   :width: 80%

Second, name your repository.
=============================

Something like 'super_awesome_killifish_data' is fine!  Be sure to check the
'initialize with Readme' box.

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
  
References:

* http://khmer.readthedocs.io/en/latest/dev/getting-started.html
* http://khmer.readthedocs.io/en/latest/dev/getting-started.html
* http://angus.readthedocs.io/en/2016/CTB-github.html
* https://monsterbashseq.wordpress.com/2016/03/08/intro-git-lab-meeting/
* https://education.github.com/guide/private_repos
* https://swcarpentry.github.io/git-novice/
* http://dib-training.readthedocs.io/en/pub/2016-02-05-intro-git.html
* https://classroom.github.com/
* http://stackoverflow.com/questions/19573031/cant-push-to-github-because-of-large-file-which-i-already-deleted
