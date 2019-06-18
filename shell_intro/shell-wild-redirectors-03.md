# Redirectors and wildcards
Here we are going to start covering some of the components of the command line that make it so versatile and powerful! If you'd like to follow along, but need to pull up the proper working environment again, revisit [here](shell-getting-started-01.html#how-to-access-the-shell-for-now) and then come back 🙂 

---

<blockquote>

To be sure we are starting in the same place, let's run: 

```bash
cd ~/shell_intro
```

</blockquote>

---

## Redirectors
When we are talking about "redirectors" here, we are referring to things that change where the output of something is going. The first we're going to look at is called a "pipe" (**`|`**). 

> A pipe (**`|`**) is used to connect multiple commands. It takes the output from the previous command and "pipes" it into the input of the following command.

Let's look at an example. Remember we used **`wc -l`** to count how many lines were in a file:

```bash
wc -l example.txt
```

And that **`ls`** lists the files and directories in our current working directory:

```bash
ls
```

If we "pipe" (**`|`**) the **`ls`** command into the **`wc -l`** command, instead of printing the output from **`ls`** to the screen as usual, it will go into **`wc -l`** which will print out how many items there are:

```bash
ls | wc -l
```

For another example, let's look at what's in the subdirectory, "data/all_samples/":

```bash
ls data/all_samples/
```

That prints out a lot of stuff, let's see how many things are in that directory: 

```bash
ls data/all_samples/ | wc -l
```

We'll get back to making sense of that when we get to *wildcards* in the next section.  

> Another important character is the greater than sign, **`>`**. This tells the command line to *redirect* the output to a file, rather than just printing it to the screen as we've seen so far.

For an example of this we will write the output of **`ls`** to a new file called "directory_contents.txt":

```bash
ls
ls > directory_contents.txt
```

Notice that when we redirect the output with the **`>`**, nothing printed to the screen. And we've just created a file called "directory_contents.txt":

```bash
ls
head directory_contents.txt
```

**It's important to remember that the `>` redirector will overwrite the file you are pointing to if it already exists.** 

```bash
ls experiment/ > directory_contents.txt
head directory_contents.txt
```

If we want to append an output to a file, rather than overwrite it like that, we can use two of them instead, `>>`:

```bash
ls >> directory_contents.txt
head directory_contents.txt
```

## Wildcards
> Wildcards as used at the command line are special characters that enable us to specify multiple items very easily. The **`*`** and **`?`** are probably the most commonly used, so let's try them out! 

<h3>The asterisk (<b>*</b>)</h3>

As we've seen, **`ls`** lists the contents of the current working directory, and by default it assumes you want everything: 

```bash
ls
```

But we can be more specific about what we're interested in by giving it a positional argument that narrows things down. Let's say we only want to look for files that end with the extension ".txt". The **`*`** wildcard can help us with that.

Here's an example:

```bash
ls *.txt
```

What this is saying is that no matter what comes before, if it ends with ".txt" we want it. 

> At the command line, the **`*`** means any character, any number of times (including 0 times). 

For a more practical example, let's change directories into that messy subdirectory we saw earlier: 

```bash
cd data/all_samples/
ls
ls | wc -l
```

So there are 900 files here, and it looks like there are 3 different extensions: ".txt"; ".tsv", and ".fq" (a common extension for the "fastq" format, which holds sequences and their quality information). 

<blockquote>
<center><b>PRACTICE!</b></center>

With 900 files and 3 file types (".txt", ".tsv", and ".fq"), we might expect there to be 300 of each type, but let's make sure. Using what we've seen above, how can we count how many files of each type there are in this directory?

<div class="toggle-header closed">
    <strong>Solution</strong>
</div>

<div class="toggle-content docutils container" style="width:100%">

<div class="highlight-bash notranslate">
<div class="highlight">
<pre>
<span class="nb">ls *.txt | wc -l</span>
<span class="nb">ls *.tsv | wc -l</span>
<span class="nb">ls *.fq | wc -l</span>
</pre>
</div>
</div>

Ah good, it's nice when things make sense 🙂
</div>
</blockquote>

So far we've just been using the **`*`** wildcard with the **`ls`** command. But wildcards can be used with many of the common shell commands we've seen so far. 

For example, we can use it with the **`mv`** command to move all 300 of the ".fq" files into their own directory at once:

```bash
ls | wc -l

mkdir fastq_files
ls fastq_files/

ls *.fq
mv *.fq fastq_files/

ls fastq_files/

ls | wc -l
```

<blockquote>
<center><b>QUICK QUESTION!</b></center>

Why does this say 601 instead of 600?

<div class="toggle-header closed">
    <strong>Solution</strong>
</div>

<div class="toggle-content docutils container" style="width:100%">

It's also counting the new directory we created 🙂
</div>
</blockquote>


> **Note:** When using wildcards, running **`ls`** first like done in the above example (**`ls *.fq`**) is good practice before actually running a command. It is a way of checking that we are specifying exactly what we think we are specifying. 

### BONUS ROUND: History!

The shell also keeps track of our previous commands for us. There are a few different ways we can take advantage of this, one is by using the **`history`** command. But that alone will print all of it to the screen. It can be more practical to "pipe" (**`|`**) that into something else like **`tail`** to see the last few commands:

```bash
history | tail
```

Or **`less`** so we can scroll through our previous commands:

```bash
history | less
```

To get out of **`less`**, press the <kbd>q</kbd> key. 

We can also use the up and down arrows at the command line to scroll through previous commands. This is useful for finding commands, but it's also useful for making sure we are acting on the files we want to act on when using wildcards. As mentioned above, we can check first with **`ls *.fq`**, press <kbd>return</kbd> to see we are acting on the files we want, and then press the up arrow to bring up the previous command, and change it to what we want without altering the "*.fq" part of the command – as we already know it's correct. Any time we can remove the chance of human error, we should 🙂


<blockquote>
<center><b>PRACTICE!</b></center>

We've already moved all the ".fq" files into their own directory. Create separate directories for the ".txt" files and the ".tsv" files too, and then try to move those files into their appropriate directories. 

<div class="toggle-header closed">
    <strong>Solution</strong>
</div>

<div class="toggle-content docutils container" style="width:100%">

<div class="highlight-bash notranslate">
<div class="highlight">
<pre>
<span class="nb">mkdir text_files</span>
<span class="nb">ls *.txt</span>
<span class="nb">mv *.txt text_files</span><br>
<span class="nb">mkdir tsv_files</span>
<span class="nb">ls *.tsv</span>
<span class="nb">mv *.tsv tsv_files</span><br>
<span class="nb">ls</span>
</pre>
</div>
</div>

It doesn't matter what the directories are named, but at the end they should be the only 3 things in the working directory 🙂
</div>
</blockquote>


<h3>The question mark (<b>?</b>)</h3>

> At the command line, the **`?`** wildcard represents *any* character that appears *only one time*. 

To see how this can be needed at times when the **`*`** won't do, let's change into the "fastq_file" subdirectory:

```bash
cd fastq_files/
```

And let's say we wanted only the ".fq" files for samples 10-19. If we tried to grab those with the **`*`**, we'd get more than we wanted:

```bash
ls sample_1*.fq
```

Because the **`*`** allows for any character *any number* of times, it is also grabbing those in the 100s. But if we use the **`?`** wildcard, which only allows any character *one time*, we get only the samples we want:

```bash
ls sample_1?.fq
```

## Summary
They may seem a little abstract at first, but redirectors and wildcards are two fundamental concepts of working at the command line that help make it a very powerful environment to work in. Just knowing they exist and generally what they do means that you can learn more about them when needed.

Next we're going to move on to [six glorious commands](shell-six-glorious-commands-04.html) that are worth knowing about 🙂

<center>

<h4><i>Special characters introduced:</i></h4>

|| | ||
||:----------:|:------------------:||
||**Characters**     |  **Meaning**  ||
|| **`|`** | a "pipe" allows stringing together multiple commands ||
|| **`>`** | sends output to a file (**overwrites** target file) ||
|| **`>>`** | sends output to a file (appends to target file) ||
|| **`*`** | represents any character appearing any number of times ||
|| **`?`** | represents any character appearing only once ||
|| | ||

</center>

---

<a href="shell-working-02.html" style="float: left"><b>Previous:</b> Working with files and directories</a>
<a href="shell-six-glorious-commands-04.html" style="float: right"><b>Next:</b> Six glorious commands</a><br>

<div style="text-align: center">
	<a href="../index.html">ANGUS Home</a>
</div>
