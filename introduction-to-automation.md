# Introduction to automation

## Keeping a log of what you ran

Let's revisit the day's commands for assembly and annotation: assuming
you've run all of the installations already, the commands to download
the data, [run the assembly](genome-assembly.html), and
[run the annotation](prokka_genome_annotation.html), are:

```
# download the data
curl -O -L https://s3.amazonaws.com/public.ged.msu.edu/ecoli_ref-5m.fastq.gz

# run the assembler
~/megahit/megahit --12 ecoli_ref-5m.fastq.gz -o ecoli

# save the assembly
cp ecoli/final.contigs.fa ecoli-assembly.fa

# get some basic metrics on the assembly
~/quast/quast.py ecoli-assembly.fa -o ecoli_report

# run prokka
prokka ecoli-assembly.fa --outdir prokka_annotation --prefix myecoli
```

(Here the '#' represent comments that the shell will ignore. If you
haven't installed the annotation pipeline yet, then put a '#' in front
of the `prokka` command so that is ignored too!)

Let's try running all of that in one go in a new directory:

```
cd ~/
mkdir work2
cd work2
```

and then copy and paste the assembly and annotation workflow
above. (This will take a few minutes.)

## While it's running...

A few points here --

First, suppose you wanted to edit this workflow, e.g.
change the name you gave to the prokka ``--prefix``?  You could put
this in a Word document, and edit it there, and *then* copy/paste everything
with that one change -- and it should still work.

Second, if you wanted to communicate to a collaborator or a computer
help person what you were doing, you could just send them all the
commands.

Third, this set of information (what you ran, what parameters you gave
it) is a pretty good addition to a Methods section, don't you
think?

## After it's done running: put it in a shell script

There are several ways of putting the commands to a script file (which is just a fancy way of describing a text file that can be run as an executable). The usual approaches are `vim` or `emacs`, but really any editor that you like works fine.

In our case, and in order to keep things more convenient, we will use RStudio to create and save the file.

First of all, we need to install RServer on the JetStream instance using the following commands (the full description is in the earlier "[Installing and running RStudio on Jetstream](http://angus.readthedocs.io/en/2017/visualizing-blast-scores-with-RStudio.html#installing-and-running-rstudio-on-jetstream)" lesson).

*Step #1*: Install prerequisites
```
sudo apt-get update && sudo apt-get -y install gdebi-core r-base
```

*Step #2*: Download and install RStudio
```
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi -n rstudio-server-1.0.143-amd64.deb
```

*Step #3*: Figure out where it's running.

Because we're using the cloud to run things, everyone will have a different
computer that they're running RStudio on.  To find out what Web site to
connect to to access YOUR server, run:

```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

*Step #4*: If you don't remember the password you set, re-set it.

After running this, copy/paste the URL into your Web browser; you
should see login page. Enter the XSEDE username and password you were
given (should be `tx160085` or `dibtiger` or `diblion` username, with
associated password).

If the login is unsuccessful, return to the terminal and run:

```
sudo passwd tx160085
```
to change your password for this instance.

You will be prompted to enter a new password:

```
Enter new UNIX password:
Retype new UNIX password:
```
but note that the text will not echo to the screen (because it's a password!)

Return to the browser login page and enter your new password. Note
this will not change the global XSEDE login info (i.e. it only affects
this instance).


Now we are ready to create the script file. Before moving to RStudio, let's first create a new folder on the Web Terminal which will contain our script as well as the results of the script.

```
cd ~/
mkdir work3
cd work3
```

Go back to the RStudio page and go to the `Files` tab on the right, and let's move to the `work3` directory (you can just click on the folder name). If for any reason the folder is not showing yet, you can always refresh the file listing by clicking on the `refresh` button (curved arrow at the top right corner of the `Files` panel).

Now that we are in the correct folder, let's create a file by going on the Menu `File` -> `New file` -> `R Script`.

Next copy all commands that we want to run as a script into our new `Untitled` (still) file:

```
# download the data
curl -O -L https://s3.amazonaws.com/public.ged.msu.edu/ecoli_ref-5m.fastq.gz

# run the assembler
~/megahit/megahit --12 ecoli_ref-5m.fastq.gz -o ecoli

# save the assembly
cp ecoli/final.contigs.fa ecoli-assembly.fa

# get some basic metrics on the assembly
~/quast/quast.py ecoli-assembly.fa -o ecoli_report

# run prokka
prokka ecoli-assembly.fa --outdir prokka_annotation --prefix myecoli
```

After this, save the file within the `work3` directory, and let's name it `awesomeScript.sh`. At this point a warning will pop up, saying that you are changing the file type (from R Script to something else). It's fine, so you can go ahead and click `yes`.

Now that we have our script, let's actually run it. Go back to the Web Terminal and into the `work3` directory, and do a full listing of the files:

```
cd ~/work3
ls -la
```

You should see the following information:
```
drwxrwxr-x  2 dibtiger dibtiger 4096 Jun 30 12:49 .
drwxr-x--- 35 dibtiger dibtiger 4096 Jun 30 12:42 ..
-rw-r--r--  1 dibtiger dibtiger  412 Jun 30 12:49 awesomeScript.sh
```

Notice the attributes of our script file `-rw-r--r--`. Seen in groups of three, from left to right, the give the information of `user`, `group` and `other` permissions. In this case, the `user` has read (`r`) and write (`w`) permissions, whereas everyone else has just read (`r`) permissions (meaning that they can see, but they can't touch :) ).

However, we actually need to be able to run the script as an executable, so we'll change the permissions to the file, for the `user` to be able to execute (`x`), using the following command:

```
chmod u+x awesomeScript.sh
ls -la
```

The listing this time around should be almost identical, with only the following change:

```
-rwxr--r--  1 fpsom fpsom  412 Jun 30 12:49 awesomeScript.sh
```

Now we are ready to run the script, as follows (it will take about 6-8'):

```
./awesomeScript.sh
```

Well done! You've just run your first script!

## Passing variables to a script

There are several cases where the script may be general enough that can be applicable to data files other than the ones it was originally created for. Or, you may just want to make slight adjustments to your script (such as names, parameters, etc) but you don't want to change the actual script every time. The way to deal with this issue, is to use arguments (AKA input parameters) to your script file.

In the case of a `'bash'` shell script, the input arguments can be defined within the script using the symbol `$` followed by a number that corresponds to the relative order in which the input parameters are provided (`$1` is the first argument, `$2` is the second and so on).

Let's make a copy of our `awesomeScript.sh` and name the copy `awesomeScript_vars.sh`, as follows:

```
cp awesomeScript.sh awesomeScript_var.sh
```

We can now see and edit the file through RStudio (_remember to refresh the listing if you don't see it directly_).

Let's try to change our script, so that the user can provide (as an argument) the name of the prefix to be used. We only need to replace the current prefix (`myecoli`) to `$1`, as follows:

```
# download the data
curl -O -L https://s3.amazonaws.com/public.ged.msu.edu/ecoli_ref-5m.fastq.gz

# run the assembler
~/megahit/megahit --12 ecoli_ref-5m.fastq.gz -o ecoli

# save the assembly
cp ecoli/final.contigs.fa ecoli-assembly.fa

# get some basic metrics on the assembly
~/quast/quast.py ecoli-assembly.fa -o ecoli_report

# run prokka
prokka ecoli-assembly.fa --outdir prokka_annotation --prefix $1
```

You can save the file in RStudio, and make it an executable on the command line as follows:

```
chmod u+x awesomeScript_var.sh
```

Finally, let us run the same script, but provide a different prefix as input, as follows:

```
./awesomeScript_var.sh angus_ecoli
```

As you can see, the argument is provided directly after the script, separated by a white space. We can let the script run (it will take ~6-7 minutes), and we can check the change in the produced output by listing the contents of the `prokka_annotation` directory:

```
ls prokka_annotation/

angus_ecoli.err  angus_ecoli.fna  angus_ecoli.gff  angus_ecoli.tbl
angus_ecoli.faa  angus_ecoli.fsa  angus_ecoli.log  angus_ecoli.tsv
angus_ecoli.ffn  angus_ecoli.gbk  angus_ecoli.sqn  angus_ecoli.txt
```


## Long-running jobs more generally

Although scripts are really handy in general, they really shine when dealing with long-running jobs - and this is indeed the case with the vast majority of bioinformatics pipelines.

The problem is that, while the script was running, we couldn't do any more work on the terminal. Which means that, if the script didn't utilize all resources available, the instance was under-utilized. Once again, there are several approaches to addressing this issue, with `nohup` and `screen` being the most widely used.

For this lesson, we will focus on using `screen`. `Screen` is a full-screen software program that can be used to multiplex a physical console between several processes (typically interactive shells). It allows a user to open several separate terminal instances inside a one single terminal window manager.

First of all, we'll need to install the tool:

```
sudo apt-get install -y screen
```

If you just type `screen` on the command line, it will show a message, but afterwards you won't notice any difference:

```
Screen version 4.03.01 (GNU) 28-Jun-15

Copyright (c) 2010 Juergen Weigert, Sadrul Habib Chowdhury
Copyright (c) 2008, 2009 Juergen Weigert, Michael Schroeder, Micah Cowan, Sadrul Habib
Chowdhury
Copyright (c) 1993-2002, 2003, 2005, 2006, 2007 Juergen Weigert, Michael Schroeder
Copyright (c) 1987 Oliver Laumann

This program is free software; you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation; either
version 3, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this
program (see the file COPYING); if not, see http://www.gnu.org/licenses/, or contact
Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111-1301

                      [Press Space for next page; Return to end.]
dibtiger@js-168-134:~/work3$
```

When you enter the screen, you can do all your work as you are in the normal command line environment. But since the screen is an application, so it have command or parameters.

Type `Ctrl-A` and `?` at the same time. Then you will see all commands or parameters on screen.

```
Screen key bindings, page 1 of 3.

                           Command key:  ^A   Literal ^A:  a

    break       ^B b           license     ,              removebuf   =
    clear       C              lockscreen  ^X x           reset       Z
    colon       :              log         H              screen      ^C c
    copy        ^[ [           login       L              select      '
    detach      ^D d           meta        a              silence     _
    digraph     ^V             monitor     M              split       S
    displays    *              next        ^@ ^N sp n     suspend     ^Z z
    dumptermcap .              number      N              time        ^T t
    fit         F              only        Q              title       A
    flow        ^F f           other       ^A             vbell       ^G
    focus       ^I             pow_break   B              version     v
    hardcopy    h              pow_detach  D              width       W
    help        ?              prev        ^H ^P p ^?     windows     ^W w
    history     { }            quit        \              wrap        ^R r
    info        i              readbuf     <              writebuf    >
    kill        K k            redisplay   ^L l           xoff        ^S s

                      [Press Space for next page; Return to end.]
```

We won't cover all these parameters, but we will focus instead on two major functions: attaching and detaching a screen.

**Detaching a Screen**

One of the advantages of screen that is you can detach it. Then, you can restore it without losing anything you have done on the screen. Here’s the sample scenario:

You are in the middle of SSH-on your server. Let’s say that you are downloading a 15GB data file for your analysis using wget command.

The download process is estimated to take 2 hours long. If you disconnect the SSH session, or suddenly the connection lost by accident, then the download process will stop. You have to start from the beginning again. To avoid that, we can use screen and detach it.

In our case, let's just clear the work3 directory and re-run our `awesomeScript.sh` as follows:

```
rm -r ecoli*
rm -r prokka*
./awesomeScript.sh
```

While downloading in progress, you can press `Ctrl-A` and `d`. You will not see anything when you press those buttons. The output will be like this:

```
[detached from 20911.pts-3.js-168-134]
fpsom@js-168-134:~/work3$
```

What happened now is that the entire process that was currently in execution, got pushed in the background so that the terminal is free for our use (however, bear in mind that the resources *are* being used, so you need to be mindful of what's actually available).

**Re-attach a Screen**

After you detach the screen, let say you are disconnecting your SSH session and going home. In your home, you start to SSH again to your server and you want to see the progress of your download process. To do that, you need to restore the screen. You can run this command:

```
screen -r
```

This will **r** estore your session, so you will be able to see exactly the output so far, as well as monitor the progress.

Several times you may have multiple screen sessions in the background. In this case, we can *list* the sessions available:

```
screen -ls
```

The output of this command lists all the screen session IDs:

```
There is a screen on:
        20911.pts-3.js-168-134  (06/30/2017 01:11:41 PM)        (Detached)
1 Socket in /var/run/screen/S-fpsom.
```

If we have the ID, we can restore a session by invoking its ID as follows:

```
screen -r 20911
```

Given that we now have just a single screen session, this is equivalent to just `screen -r`, but usually we will have more than one running at a given time.

**Destroying a screen**

There are 2 (two) ways to leaving the screen. The easiest one (and the most commonly used) is to use the `exit` command in the screen, which will lead to the following output:

```
[screen is terminating]
fpsom@js-168-134:~/work3$
```

Alternatively, you can use `Ctrl-A` and `K` to kill the screen (nuclear option).


**Logging a Screen**

Finally, sometimes we will use a screen session to run multiple scripts and/or commands. In this case, it would be useful to keep log of whatever we input in screen. This is done when first running the screen tool by adding the parameter `-L`. This will create a log file named `screenlog.0` in the working directory, that will capture everything you do while you are in the screen window.

Let's start a new screen to re-run our script, which will include also the initial "clean-up" commands:

- Start the screen with logging

```
screen -L
```

- Run the clean-up commands and the script.

```
rm -r ecoli*
rm -r prokka*
./awesomeScript.sh
```

- Detach the screen with `Ctrl-A` and `d`.

When the screen finishes (check by restoring), you can exit the screen (`exit`) and have a look at the log file (`more screenlog.0`). It should contain something similar to the following:

```
fpsom@js-168-134:~/work3$ rm -r ecoli*
fpsom@js-168-134:~/work3$ rm -r prokka*
fpsom@js-168-134:~/work3$ ./awesomeScript.sh
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
 30  412M   30  126M    0     0  39.7M      0  0:00:13  0:00:02  0:00:11 29.7M
 77  412M   77  320M    0     0  42.6M      0  0:00:09  0:00:07  0:00:04 49.9M
100  412M  100  412M    0     0  44.7M      0  0:00:09  0:00:09 --:--:-- 50.7M
15.671Gb memory in total.
Using: 14.104Gb.
MEGAHIT v1.1.1-2-g02102e1
--- [Fri Jun 30 14:18:28 2017] Start assembly. Number of CPU threads 6 ---
--- [Fri Jun 30 14:18:28 2017] Available memory: 16826556416, used: 15143900774
--- [Fri Jun 30 14:18:28 2017] Converting reads to binaries ---
b'    [read_lib_functions-inl.h  : 209]     Lib 0 (ecoli_ref-5m.fastq.gz): interleaved, 5
000000 reads, 100 max length'
b'    [utils.h                   : 126]     Real: 9.1392\tuser: 3.1520\tsys: 0.3400\tmaxr
ss: 158660'
--- [Fri Jun 30 14:18:37 2017] k-max reset to: 119 ---

...
...
...
...
```

An alternative way of doing this is by `Ctrl-A` and `H`. At the bottom left of the screen, there will be a notification that tells you like: Creating logfile “screenlog.0“.  To close screen to log running activity, press `Ctrl-A` and `H` again.

> **Note** _Please be careful, we use capital `H` letter. Using non capital `h`, will only create a screenshot of screen in another file named hardcopy._


## Scripts: 'bash' shell scripts, R scripts, and Python scripts

The script file that we have been using so far, is essentially a `'bash'` shell script. The same approach can be used for any type of scripting language, with the variation being only in the way they are ultimately executed.

For example, if we have an R script named `rrules.R`, it can be executed at the command line level as follows:

```
Rscript rrules.R
```

Similarly, if we have a Python script named `snakes.py`, the execution will be:

```
python snakes.py
```

We've used this in various places (the `gather-counts.py` Python script and
the `yeast.salmon.R` R script in [RNAseq](counting.html)) as a way to keep
you from having to switch too many times from the shell to RStudio and
back.  But it is also just a convenient way to track what you're doing
and communicate it to yourself and others!

## Final thoughts

Automation like this is good for a few reasons.

1. Efficiency. You don't have to remember and/or mistype commands any more!  You just need to remember where you put the script!

2. Efficiency x 2: new projects! When you start a new project you can modify your old scripts to work in the new way!

3. Efficiency x 3: revisions! When reviewers *finally* get back to you
   and tell you to make all your red solid lines into blue dashed
   lines, you can edit and then re-run the R script that created the
   figure.

4. Repeatability and communication: graduate and your advisor wants to pass
   your methods on to the next student? Send them the scripts! Collaborator
   wants to rerun things with different parameters? Send them the script!
   etc.
   
5. Tracking and versioning. We didn't show you how to do this *yet*,
   but with script files it's easier to track what you did last time
   and compare it to how you're doing it now.  When you find yourself
   naming things "myscript.mar2017.sh" and "myscript.jun2017.sh", then
   you should talk to someone about something called "version control."
   We'll show you more on Monday!

...and probably some other things, too.
