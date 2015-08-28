======================================
GitHub, Pull Requests, and ReadTheDocs
======================================

This lesson will walk through the powerful combination of two Web
sites -- `github.com <http://github.com/>`__ for collaboratively
editing documents, and `readthedocs.org <http://readthedocs.org>`__
for hosting documentation and workshop Web sites.  You're probably all
familiar with either GitHub or a similar code-hosting Web site,
BitBucket; ReadTheDocs is what we use to host this site, and it's
based around the `Sphinx documentation generator
<http://sphinx-doc.org/index.html>`__.  Sphinx is used to build the
Python documentation, among other things.

Alternatives: You could do very similar things with Markdown and
`github pages <https://pages.github.com/>`__, or a variety of other
systems.  The overall workflow is largely the same, however!

**Learning goals:**

At the end of this lesson, students will understand:

* forking, online editing, and pull requests on github;

At the end of this lesson, students will be able to:

* connect a sphinx site to readthedocs;

* collaboratively update sphinx sites through pull requests;

The lesson
----------

:author: Titus Brown
:date: Aug 28, 2015

1. Go to http://github.com/ngs-docs/angus/ and click "Fork" in the upper right.

   Save this URL - this is YOUR ANGUS REPOSITORY.

2. Go to http://www.readthedocs.org/ and set up a ReadTheDocs site for
   YOUR ANGUS REPOSITORY.

   1. Log in
   2. Select "import a project";
   3. Click "connect to github";
   4. Authorize the application;
   5. Click "Import a github repository";
   6. Click "Sync projects;"
   7. Find YOUR ANGUS REPOSITORY;
   8. Click "create";
   9. Change the name to 'angus-USERNAME';
   10. click "next";

   Now, view the docs (upper right) and change the URL from 'latest'
   to 'stable'.  This is YOUR ANGUS READTHEDOCS site.

3. (Everybody wait while Titus adds something to ngs-docs/angus/.)

4. Merge in Titus's update into their own repository.

   1. Go to YOUR ANGUS REPOSITORY.
   2. Click "Compare" (above file listing, to the right);
      You should see "nothing to compare".
   3. Click on "switching the base";
   4. Click "Create pull request";
   5. Click "Create pull request" again;
   6. Click "Merge pull request";
   7. Click "Confirm merge"

5. Everyone should now be able to see MY change at YOUR ANGUS REPOSITORY,
   and, after a few minutes, at YOUR ANGUS READTHEDOCS SITE.

6. Now, everybody fix something in YOUR ANGUS REPOSITORY.

   1. Go find a file you want to edit (suggest ``week3.rst``, add affiliation);
   2. Edit, commit.
   3. Wait for YOUR ANGUS READTHEDOCS SITE to update.

7. Submit a pull request.

   1. Go to a directory view in YOUR ANGUS REPOSITORY.
   2. Click "compare";
   3. Click "Create pull request";
   4. Click "Create pull request" again;

8. (Wait for Titus to merge in your pull request.)

9. Changes will now appear on the `main angus RTD site
   <http://angus.readthedocs.org/en/2015/>`__

10. WINNAGE.

-----

A short guide to creating your own Sphinx repository
----------------------------------------------------

Create a virtualenv::

   python -m virtualenv env
   . env/bin/activate

Install sphinx::

   pip install sphinx

Create a new directory for your docs::

   mkdir my-project
   cd my-project

Run sphinx-quickstart::

   sphinx-quickstart

(accept default for most questions; answer those you can't.)

Build::

   make html

Marvel at the contents of ``_build/html/index.html``.

Now, put into git::

   git init
   git add Makefile index.rst conf.py
   git commit -am "initial commit"

Go to github, create a new github repo, and follow the instructions
to "push an existing repository from the command line".

Then connect that repository to ReadTheDocs as above.
