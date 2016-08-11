..  _01_Download_and_test:


Download and Test
=========================

Download the data
-----------------

Working with genomics files almost always requires using the shell, so
although this isn't the official shell lesson, we're going to use many
BASH functions. If you've already covered the shell lesson, consider
this a refresher.

First we're going to download a zip file that is saved in a public
Dropbox folder. We're going to download using either ``wget`` or
``curl`` and which one you need to use depends on your operating system,
as most computers will only have one or the other installed by default.
To see which one you have:

.. code:: bash

    which curl
    which wget

``which`` is a BASH program that looks through everything you have
installed, and tells you what folder it is installed to. If it can't
find the program you asked for, it returns nothing, i.e. gives you no
results.

On OSX, you'll likely get:

::

    $ which curl
    /usr/bin/curl
    $ which wget
    $ 

Once you know whether you have ``curl`` or ``wget`` use one of the
following commands to download the zip file:

.. code:: bash

    cd
    wget https://www.dropbox.com/s/d7zitckb5fz8494/GenomicsLesson.zip

or

.. code:: bash

    cd
    curl -LO https://www.dropbox.com/s/d7zitckb5fz8494/GenomicsLesson.zip

Checking the data
-----------------

Most genomics files that you download will be very large, which can make
them more prone to download errors. Flaky internet can make you think
that your file is finished downloading, even if it really just stopped
partway through. Furthermore, even the files that are filled with
familiar characters are nearly impossible to fact check by eye. Imagine
trying to make sure your favorite genome had downloaded correctly by
manually comparing each base on your computer to the one at NCBI.
Luckily theres a better way.

Whenever you download a large or important file, you should check to
make sure that it is an exact match to the copy online. The most common
way to do this is to run the file through a cryptographic hash function
which processes all of the information in the file through a complex
algorithm to get a hash value. A hash value looks a bit like an ideal
password: a random looking mix of letters and numbers. Because the hash
function is a very complex equation, in theory, for any given hash
function, every unique input will have a unique hash value. So if you
get the same hash value from two different files, those files are
identical.

There are many available hash functions, and as computers get more
sophisticated, older ones become easier to hack, so there will likely
alwasy be new ones. But right now, we're going to use MD5 because it is
pre-installed on most computers, and is most likely the one your
sequencing facility or collaborator will send you when they give you
large files.

There are two common versions of MD5, and we're going to use ``which``
again to see which one you have installed:

.. code:: bash

    which md5
    which md5sum

Then get get the hash value for the zip file you downloaded by running
either:

.. code:: bash

    md5 GenomicsLesson.zip

or

.. code:: bash

    md5sum GenomicsLesson.zip

The file I uploaded gave this answer:

.. code:: bash

    md5 GenomicsLesson.zip
    MD5 (GenomicsLesson.zip) = e237e30985867e6bea741949e42a0c3b

::

    Exercise 1: Download the best practices file from Data Carpentry:
    https://raw.githubusercontent.com/ACharbonneau/2016-01-18-MSU/gh-pages/Files/GoodBetterBest.md

:ref:`02_working_with_archives`
