Getting to know your data
=========================

Looking at data files
---------------------

Each step of processing and analysis of a genomics pipeline spawns many
new files, of many types. Some filetypes, like GFF are only found in a
single step of the pipeline, and so are relatively easy to keep track
of. However, most are more like FASTQ files, where any given file could
be from many different steps of the pipeline. These are the ones that
cause the most trouble, and need the most careful management.

.. figure:: ../Files/ngs_map_read_file_formats.png
   :alt: FileTypes

   FileTypes

First, lets see what data files we have available:

.. code:: bash

    ls

``ls`` stands for list, and if you call it all by itself, it just
returns a list of whatever is inside the folder you're currently looking
at. You should give you a fairly big list of files in alphabetical
order. However, they're hard to understand like this, so lets ask ``ls``
to make the list a little easier to read:

.. code:: bash

    ls -lah

Here we've added modifiers to ls. Computer people usually call these
modifiers 'flags' or 'arguments', and here we've added 3 flags: ``-l``
directs ``ls`` to give us the results in 'long format' so we get more
information ``-a`` tells ``ls`` to show us 'all' of the things in the
folder, even if they're usually hidden ``-h`` makes the output 'human
readable', so you see file sizes in kb or gb instead of bytes

Genomics Text files
~~~~~~~~~~~~~~~~~~~

FASTA
^^^^^

You're likely already familiar with FASTA files, as this is the most
common way to distribute sequence information. Let's look at one:

.. code:: bash

    head Raphanus.fa

``head`` is another program, and it shows you just the top few lines of
a file. By default, it shows ten, (so five sequences) but we can also
change that behavior with flags:

.. code:: bash

    head -4 Raphanus.fa

Now, you should see the first four lines of the Raphanus.fa file.

    Exercise Try looking at EV813540.fa

FASTA files always have at least one comment line, which almost always
begins with ">", but can start with ";". A given sequence in the file is
allowed to have multiple comment lines, but they usually don't. Extra
comment lines for sequences can break some downstream processes.

After the comment line is the sequence. Usually this is all on one line,
but you can see that this one is formatted so that each sequence line is
only 80 characters wide. This makes it easy to read, but makes it
slightly more difficult to search within the file. For searching, its
nice to have files where all of each sequence is on a single line. For
instance, lets see whether there are any EcoRI sites are in the
Raphanus.fa file:

.. code:: bash

    grep "GAATTC" Raphanus.fa

grep is a program that searches for any string, and by default returns
the entire line that your string is found in. For a file this big, this
isn't very helpful. So lets modify how grep reports it's findings:

.. code:: bash

    grep -B 1 "GAATTC" Raphanus.fa

``-B number`` grep will return the line with your string plus 'number'
lines of 'before context', so here we'll get one previous line...the
comment that tells us the sequence name

Now we know which of the sequences have the restriction site we're
looking for, but there's so many they've overfilled the screen. So lets
redirect the output from the screen into a file:

.. code:: bash

    grep -B 1 "GAATTC" Raphanus.fa > Raphanus_EcoRI.fa

The greater than sign takes everything that happens on this side of it
``>`` and dumps it into the place designated here. So, all of the output
from that ``grep`` command above got saved into a new file called
Raphanus\_EcoRI.fa Since we didn't specify a place to save it, the new
file is just saved in the same folder we're in, and we can see it by
using ``ls`` again:

.. code:: bash

    ls -latr

``-r`` makes the list print to our screen in reverse chronological
order, so the newest files are on the bottom. This makes it easier to
find what we're looking for.

``grep``, ``ls`` and ``head`` all have lots of useful flags, and we can find out what they
are by looking at the manual page:

.. code:: bash
	man grep

This opens the manual in the text viewer ``less``, which we'll talk about more in a few
minutes. For now, the important things to know are that you can scroll line by line
using the arrow keys, or go down one page at a time using the space bar. You can search 
for a keyword by typing ``/`` and text to search for. Let's look at the explanation for 
a flag we already used:

.. code:: bash
	/-B

I actually prefer to look at man pages online, because searching them is easier. 
Try `Googling 'man grep' <http://www.google.com/search?q=man%20grep>`_

.. code:: output
Exercise: How would you change ``grep -B 1 "GAATTC" Raphanus.fa > Raphanus_EcoRI.fa`` 
to add line numbers to the output? Hint1_.

So, now we can make a file that only has sequences with our cut site. Depending on what 
and why you're searching, this might be useful for making markers or primers. But maybe 
we just want to know how many sequences had our cut site:

.. code:: bash

    grep -c "GAATTC" Raphanus.fa

``-c`` grep 'counted' 88 instances of EcoRI

Grep happens to have a built in flag for counting matches, but many other programs don't. 
So there is a separate program just for counting that we could use by invoking a 'pipe':

.. code:: bash

	grep "GAATTC" Raphanus.fa | wc

``wc`` stands for word count, and actually gives us three numbers: number of lines, number
of words and number of characters, in that order. The first two are both 88 because there
are no spaces between the letters of the sequences, so each sequence is interpreted as one
big word.

If we only want one of those numbers, we can use the flags ``-l``, ``-w``, and ``-c`` 
respectively. 

A 'pipe' is a little like holding up a real-world pipe, everything you dump in the top 
comes out the bottom. Here, the answer from ``grep "GAATTC" Raphanus.fa`` goes in and 
becomes input for ``wc``. Notice that we only told the computer which file to use for 
``grep``, each pipe after that (there can be an as many as you want) gets its input from
the previous programs output. Also notice that we got rid of all of the grep flags. Why?

.. code:: output
Exercise: How would you get *just* the *names* of the sequences that match our restriction
site? Hint2_. And save that list to a file?

What if we want to do a 'fuzzy' search. Say we want to 





.. _Hint1:
Use a ``-n``
.. _Hint2: 
You'll need two greps



