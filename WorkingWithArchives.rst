Working with Compressed files
=============================

As previously mentioned, genomics data files tend to be *large*. Since
larger files are slower and more costly to move around, you will often
encounter files that have been *compressed* to save time/space/money.
The two most commonly encountered types of compressed files are Zip
archives (e.g. ``filename.zip``), Gzip archives (e.g. ``filename.gz``)
and Tarballs (e.g. ``filename.tar`` or ``filename.tar.gz``).

Zip 
-----------------

Once you've convinced yourself that the file you have is the file that
you *ought* to have, the next thing that you'll want to do is unzip it
(a.k.a. uncompress or decompress or extract). You can unzip your .zip
archive using the ``unzip`` program:

.. code:: bash

    unzip <filename.zip>

If you don't want to extract everything, but rather check the contents,
you can view what a zip contains using the ``-l`` flag ('list'):

.. code:: bash

    unzip -l <filename.zip>

Right now neither of these work, because unzip (and zip) are not
installed by default. However, Ubuntu recognizes that unzip is an
available program, and has given us a prompt telling us how to get unzip
if we want it

.. code:: bash

    sudo apt-get install unzip

When you want to go in the other direction and make your own archive the
command is simply ``zip``. It works like this:

.. code:: bash

    zip <mynewarchive.zip> <myfirstfile.txt> <mysecondfile.sam>

Note that you can also use the ``-r`` flag (recursive) to zip up a
folder and all its contents, including subfolders like so:

.. code:: bash

    zip -r <myproject.zip> myproject/

GZip
-----------------

If you have been sent a big bundle of data as a gzip archive, then
happily the same procedure applies for viewing and extracting as with
zip archives, but with the ``gunzip`` program:

.. code:: bash

    gunzip -l <bundle.gz>
    gunzip <bundle.gz>

Tarball
-----------------

Things are slightly different (read 'complex') if you encounter a
tarball: ``thisfile.tar`` or ``thatfile.tar.gz`` or ``tacofile.tgz``.

.. figure:: http://imgs.xkcd.com/comics/tar.png
   :alt: I'm so sorry...

You can view the contents of tarballs using the ``tar`` program:

.. code:: bash

    tar -tf <thisfile.tar>
    tar -ztvf <thatfile.tar.gz>
    tar -ztvf <tacofile.tgz>

...and extract them like this:

.. code:: bash

    tar -xf <thisfile.tar>
    tar -zxvf <thatfile.tar.gz>
    tar -zxvf <tacofile.tgz>

Other types of compressed files and archives do exist, but these are the
most common.
