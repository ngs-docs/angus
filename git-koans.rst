====================================
Automation, scripts, git, and GitHub
====================================

Start up an Ubuntu 14.04 instance and run ::

   sudo bash

   apt-get update
   apt-get -y install screen git curl gcc make g++ python-dev unzip \
           default-jre pkg-config libncurses5-dev r-base-core \
           r-cran-gplots python-matplotlib sysstat

Automation and scripts
======================

Suppose you want to run a few different commands in succession -- this is
called "automating tasks".  See http://xkcd.com/1205/

Key caveat: none of the commands require any user input ('y', 'n', etc)

#. Figure out what commands you want to run.

#. Put them in a file using a text editor, like nano or pico,
   or TextEdit, or TextWrangler, or vi, or emacs.

#. Run the file.

For example, try putting the following shell commands in a script called
'test.sh'::

   cd /root
   curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
   unzip Trimmomatic-0.32.zip
   cp Trimmomatic-0.32/trimmomatic-0.32.jar /usr/local/bin
   cd /mnt
   curl -O http://athyra.idyll.org/~t/subreads-A-R1.fastq.gz
   curl -O http://athyra.idyll.org/~t/subreads-A-R2.fastq.gz

   java -jar /usr/local/bin/trimmomatic-0.32.jar PE subreads-A-R1.fastq.gz subreads-A-R2.fastq.gz s1_pe s1_se s2_pe s2_se LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25

If you want to do this ghetto style, type::

   cd /root
   cat > test.sh

then paste in the above, put a newline at the end, and then type CTRL-D.

Now, run::

   bash test.sh

...and this will run all of the stuff that you put in test.sh.

This is called a shell script (because it is a *script* for running things,
written in *shell language*).  Adina has already told you about Python
scripts a bit.

You can do this kind of "scripting" with any set of commands.  Just
keep track of exactly what you're teaching at the command line by
pasting it into a document somewhere, then save the document in text
format, and that's your script!

----

Note: what happens when you run 'bash test.sh' again?

Another note: notice the explicit 'cd' steps... why?

But now! You're stuck keeping track of all of these files.

See: http://www.phdcomics.com/comics/archive.php?comicid=1531

This is what version control is for.

Here are some things to try.

Some git koans
==============

Forking a repository on github
------------------------------

1. **In a browser**, log into github.com.

2. Go to https://github.com/ngs-docs/ngs-scripts and click "fork" (upper right).  This will make a copy of that git repository under *your* account.  It should leave you at http://github.com/<YOUR ACCOUNT>/ngs-scripts.

3. Select the URL about midway down the page ('https://...') and copy
   it to your clipboard.  Hint: There's a handy little button on the
   right to do this.

4. Go to your EC2 command line, and type::

      git clone https://github.com/<YOUR ACCOUNT>/ngs-scripts.git

   where the last bit is pasted from what you copied in step 3.

5. Change into the ngs-scripts directory::

      cd ngs-scripts/

   and poke around.

6. Marvel.  Note that what is in your directory is the same as what
   you can see via the github interface.

7. **In a browser**, go back to your copy of ngs-scripts.  Select 'README.md'
   in the top-level directory.

8. Select 'Edit'.

9. Change something in the text box (e.g. add "Kilroy was here.")

10. Click "Commit changes".

11. Note that in the browser, README.md has been updated.

12. **In the command line**, note that README.md hasn't changed.  The
    repositories are distinct and separate.

13. Type::

       git pull https://github.com/<YOUR ACCOUNT>/ngs-scripts.git master

    to *pull* the changes from github into your local copy.

14. Now README.md is the same in both places!!

What you have done here is *cloned* your repository, then *edited*
your file in the original repository, and then *pulled* the changes
from the original repository into your new repository.

Create a new file on github and edit it, then pull
--------------------------------------------------

Note the '+' after the directory name (next to 'branch: ') just above
the list of files.

This gives you the opportunity to create and edit a *new* file.

Do so, 'commit new file', and then do step #13 above.

Now you've created a new file on github!

Edit local file and push to github
----------------------------------

**At the command line**,

1. Edit the README file (either with a local editor like 'pico', or with
   Dropbox, or something; e.g. do::

      cp README.md ~/Dropbox
      (edit it)
      cp ~/Dropbox/README.md .

   to update it remotely and copy it back over).  Use 'more' to make sure
   your local copy is different.

2. Type::

	git diff

   to see your changes.  The lines with '+' at the beginning are your
   new changes, the lines with '-' at the beginning are what they
   replaced.

2. Type::

	git commit -am "made some changes"

   to commit the changes as things you want to do.

   (At this point, you could also type 'git checkout README.md' to replace
   the changed file with the original.)

3. Type::

	git push https://github.com/<YOUR ACCOUNT>/ngs-scripts.git master

4. Marvel that the local changes are now viewable on github.com directly!

What you have done here is to edit files in one repository, and then
*pushed* the changes to another (remote) repository.

Create a new repository; add some files to it.
----------------------------------------------

Let's create a new repository, just for you.

**In a Web browser**,

1. Go to http://github.com/ and click on "New repository."

2. Make up a repository name (it will suggest one; ignore it.)

3. Select the "initialize this repo with a README."

4. Select 'Create repository."

5. Now, clone it to your EC2 machine::

      git clone https://github.com/<YOUR ACCOUNT>/<YOUR REPO NAME>.git

6. Change into the new repo directory::

      cd <YOUR REPO NAME>

7. Create a new file::

      echo hello world > greetings.txt

8. Add it to your repository::

      git add greetings.txt

9. Commit it::

      git commit -am "added greetigs"

10. Push it to your github repository::

      git push https://github.com/<YOUR ACCOUNT>/<YOUR REPO NAME>.git master

11. Go check it out on the Web -- do you see greetings.txt?
