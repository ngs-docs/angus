# Automating workflows using bash

## Objectives

+ Write a shell script

We have now worked through two workflows, RNA-seq analysis and variant calling.
In both cases, we pasted all of our commands into the terminal, and used for 
loops to run each step of the workflow over our six samples. This is a great way
to understand how each tool works and to establish a workflow you want to use.

Let's imagine now that we decided to download all 672 samples from this dataset
instead of only working with our six samples. We already know how each tool in
our workflow works, and we know how to run them. We can now use bash scripting
to automate our workflow. 

We are going to automate the quality control portion of our workflow. To write 
a script to run our FastQC analysis, we just need to take each of the commands we entered 
to run FastQC and process the output files, and then put them into a single file with 
a `.sh` extension. The `.sh` is not essential, but serves as a reminder to 
ourselves and to the computer that this is a shell script.

But that's it! A bash script (shell script) is just like when we type things one command at a time at the command line, except they are saved in a file!

## Accessing our JetStream intances

You should still have your jetstream instance running, you can follow the instructions [here](jetstream/boot.html) to log in to JetStream and find your instance. Then `ssh` into it following the instructions [here](jetstream/boot.html#ssh-secure-login).

Next, we're going to launch RStudio again from our instances. RStudio will be convenient for being able to edit a text file alongside an active terminal window (command-line window). We'll see this in action in a second. As earlier, run the following link to produce
a url that navigates to RStudio when entered in your browser.

```
echo http://$(hostname -i):8787/
```

## Starting with automation
Previously, we've been writting R code in the RStudio text-editor at the top left, but we can also write any kind of files there. So that will be where we build up our bash script. To start with a clean new document, select "File" from the top menu bar inside RStudio, then "New File", then "Text File". 

At the bottom left window is our "console" right now. This is our R console. But if we select "Tools" from the top menu bar within RStudio, then "Terminal", and then "New Terminal", we will open a command-line terminal in that bottom left window that looks something like this:

<center><img src="_static/rstudio-terminal.png" width="90%"></center>
<br>

This is the same as our terminal that we used to `ssh` into and then get the link to the RStudio environment, but we are just accessing it through RStudio here. Moving around in there, running commands like `cd`, `ls`, and others all work the same as if we were in a terminal outside of RStudio. If something isn't working there, double-check it is on the "Terminal" tab, and not the "Console" tab (the "Console" tab is the R environment, while the "Terminal" tab is our command-line environment).


## Our first bash script!
First, in the Terminal of RStudio, let's run these commands:

```bash
ls
echo "I'm here to show something random being done"
pwd  
```

Those 3 commands just execute as we are used to. Now let's copy those 3 commands, and put them into our text editor window at the top-left:

<center><img src="_static/bash-script-01.png" width="90%"></center>
<br>

And save the file in our home directory (which should be where it starts) as "my-first-bash-script.sh":

<center><img src="_static/bash-script-01a.png" width="90%"></center>
<br>

Now if we check at the command line, (our terminal window at the bottom left), we can see this file is there: 

```
ls *.sh
```

And we can run it by providing the command `bash`, and then the file as a positional argument like so:

```
bash my-first-bash-script.sh
```

And all three lines are executed as though we entered them one at a time! This is all a bash/shell script is 🙂

Next let's build up a larger script!

## Constructing a bash script
When we write a bash script, we need to add *all* commands that we ran in our
workflow. This includes making and changing directories, moving files, and 
running analysis programs like `fastqc`. 

Let's start by changing into the right directory and running `fastqc`. 

**Don't enter these commands in your Terminal, instead copy and paste them into a new text document in the top left panel (File -> New File -> Text File). Then save the file as "qc.sh".**


```
cd ~/data/

fastqc *fastq.gz
```

<center><img src="_static/bash-script-02.png" width="90%"></center>
<br>

Now we have a bash script that automates `fastqc`! 

**Let's run the bash script from our terminal in RStudio:**

```
bash qc.sh
```

Let's add more to our script. Next we'll organize our output files. **Add this text onto our file in the text-editor panel where our script is building, not in the Terminal window.**

```
cd ~/data/

fastqc *.fastq.gz

mkdir -p ~/fastqc_untrimmed

mv *.zip ~/fastqc_untrimmed/
mv *.html ~/fastqc_untrimmed/
``` 

We can also add a few pointers that
will print to our screen to let us know what point in the workflow we are in.


```
cd ~/data/

echo "Running FastQC ..."
fastqc *.fastq.gz

mkdir -p ~/fastqc_untrimmed

echo "Moving FastQC results..."
mv *.zip ~/fastqc_untrimmed/
mv *.html ~/fastqc_untrimmed/
```

And don't forget to save the file, the filename should turn black after being saved and look like this: 

<center><img src="_static/bash-script-03.png" width="90%"></center>
<br>


Now, let's run this again. 

```
bash qc.sh
```

We now see that our echo messages print to the terminal and tell us where we
are in the workflow. After it finishes, if we run `ls` on the "fastqc_untrimmed" directory, we see our `fastqc` output files are in there:

```
ls ~/fastqc_untrimmed
```


<blockquote>
<center><b>PRACTICE!</b></center>

Add some commands to the end of our script to: 1) change directories into the "fastqc_untrimmed" directory; 2) print out to the screen that we are running "Multiqc"; and 3) then run <code>multiqc</code> (the way to run <code>multiqc</code> is to execute <code>multiqc .</code>).

<div class="toggle-header closed">
    <strong>Solution</strong>
</div>

<div class="toggle-content docutils container" style="width:100%">

Our entire file should now look something like this:

<div class="highlight-bash notranslate">
<div class="highlight">
<pre>
<span class="nb">cd ~/data/

echo "Running FastQC ..."
fastqc *.fastq.gz 

mkdir -p ~/fastqc_untrimmed

echo "Saving FastQC results..."
mv *.zip ~/fastqc_untrimmed/
mv *.html ~/fastqc_untrimmed/

cd ~/fastqc_untrimmed/

echo "Running multiqc..."
multiqc .</span>
</pre>
</div>
</div>

</div>
</blockquote>

Now let's run our full script one last time! Note we need to change back into our home directory again first, because the last script we ran did not have that part added yet. **This is in the Terminal at the bottom left panel again**

```
bash qc.sh
```

And after that finishes We should now see the `multiqc` output files in this directory.

```
ls ~/fastqc_untrimmed
```

Bash scripting can be very powerful!
