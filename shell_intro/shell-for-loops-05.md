# Variables and For loops

<blockquote>
<center><b>Welcome to the wonderful world of loops!!</b></center>
</blockquote>

Loops are extremely powerful in all programming languages. They are what let us write out a command or operation once, and have it run on all of our samples or files or whatever we want to act on. Not only is this powerful, but it also helps with keeping our code more concise and readable, and it helps elmininate some more of our mortal enemy (human error). There are multiple types of loops, but here we are going to cover what is probably the most common type: the for loop. First, let's quickly go a tiny bit about variables. 

If you'd like to follow along, but need to pull up the proper working environment again, revisit [here](shell-getting-started-01.html#how-to-access-the-shell-for-now) and then come back 🙂 


---

<blockquote>

To be sure we are starting in the same place, let's run: 

```bash
cd ~/shell_intro
```

</blockquote>

---


In order to use For loops, we first need a little understanding of variables.

## Variables
We can think of a variable as a placeholder for a value that will change with every iteration of our loop. To set a variable at the command line, we need to provide the variable name we want, an equals sign, and then the value we want the variable to hold (with no spaces in between any of that). Let's try it:

```bash
my_var=ANGUS
```

Nothing prints out when a variable is set, but the value "ANGUS" has been stored in the variable "my_var". 

To use what's stored in a variable, the variable name needs to be preceded by a **`$`** so the shell knows to *evaluate* what follows, rather than just treat it as generic characters. 

To see this in action, we'll use the **`echo`** command. **`echo`** is a command that prints out whatever is provided to it (it turns out this is actually really useful in some places – like in a program, for example, to report information to the user). Here we'll use it to check what's being stored in our variable:

```bash
echo $my_var
```

Note that if we don't put the **`$`** in front, **`echo`** just prints out the text we gave it:

```bash
echo my_var
```

Recall that spaces are special characters on the command line. If we wanted to set a variable that contained spaces, we could surround it in quotations to tell bash it should be considered as one thing:

```bash
my_new_var="ANGUS is awesome."

echo $my_new_var
```

Great, that's really all we need to know about variables for now. Let's get to the good stuff 🙂

## For loops
Let's make a new directory to work in: 

```bash
mkdir for_loops
cd for_loops/
```
### The 4 magic words
There are 4 special words in the syntax of a For Loop in the shell: `for`, `in`, `do`, and `done`. 

|| | ||
||:----------:|:------------------:||
||**Magic word**     |  **Purpose**  ||
|| **`for`** | set the loop variable name ||
|| **`in`** | specify whatever it is we are looping over ||
|| **`do`** | specify what we want to do with each item ||
|| **`done`** | tell the shell we are done telling it what to do with each item ||
|| | ||

Let's see what this looks like in practice. Here we are going to: name the variable "item" (we can name this whatever we want); loop over 3 words (cat, dog, and ukulele); and we're going to just `echo` each item, which will print each word to the terminal. 

```bash
for item in cat dog ukulele
do
  echo $item
done
```

> **Note:** Notice the prompt is different while we are within the loop syntax. On the jetstream instances we're working on, it shows `and...> ` until we finish the loop with **`done`**. This might look different if you are on a different system, but it will be something distinct from the normal prompt. If you get stuck with that alternate prompt and you want to get rid of it, you can press <kbd>ctrl + c</kbd> to cancel it. 

Just to note, we don't need to put these on separate lines, and we don't need to indent over the "body" of the loop like we did above (the `echo $item` part), but both can help with readability so we will continue doing that moving forward. As an example though, we could also enter it like this on one line, separating the major blocks with semicolons:

```bash
for word in cat dog ukulele; do echo $word; done
```

We can also do multiple things within the body of the loop (the lines between the special words **`do`** and **`done`**). Here we'll add another line that also writes the words into a file we'll call "words.txt":

```bash
for item in cat dog ukulele
do
  echo $item
  echo $item >> words.txt
done
```

Now we created a new file that holds these words:

```bash
ls
head words.txt
```

<blockquote>
<center><b>PRACTICE!</b></center>

Notice that we used <b>>></b> as the redirector inside the loop, and not <b>></b>. Why do you think this is? Try running the loop with the <b>></b> redirector instead and writing out to a new file (instead of "words.txt", call it anything else).

<div class="toggle-header closed">
    <strong>Solution</strong>
</div>

<div class="toggle-content docutils container" style="width:100%">

<div class="highlight-bash notranslate">
<div class="highlight">
<pre>
<span class="nb">

for item in cat dog ukulele
do
  echo $item
  echo $item > test.txt
done
</span>
<span class="nb">head test.txt</span>
</pre>
</div>
</div>

Since <b>></b> overwrites a file, each time we go through the loop it would overwrite the file and at the end we'd be left with just the last iteration, and we'd have a file holding only "ukulele".

</div>
</blockquote>

Usually we won't want to type out the items we're looping over, that was just to demonstrate what's happening. Often we will want to loop through items in a file, like a list of samples or genomes. 

### Looping through lines of a file
Instead of typing out the elements we want to loop over, we can execute a command in such a way that the output of that command becomes the list of things we are looping over. 

We're going to use the **`cat`** command to help us do this (which comes from con**cat**enate). **`cat`** is kind of like **`head`**, except that instead of just printing the first lines in a file, it prints the whole thing:

```bash
cat words.txt
```

Here we'll use **`cat`** to pull the items we want to loop over from the file, instead of us needing to type them out like we did above. The syntax of how to do this may seem a little odd at first, but let's look at it and then break it down. Here is an example with our "words.txt" file we just made:

```bash
for item in $(cat words.txt)
do
  echo $item
done
```

Here, where we say **`$(cat words.txt)`**, the command line is performing that operation first (it's evaluting what's inside the parentheses, similar to what the dollar sign does when put in front of our variable name, "item"), and then puts the output in its place. We can use **`echo`** to see this has the same result as when we typed the items out:

```bash
echo $(cat words.txt)
```

For a more practical example, let's pull multiple specific sequences we want from a file!

#### BONUS ROUND: interleaving files with **paste**
A pretty neat use of **`paste`** is to interleave two files. What **`paste`** is doing is sticking two files together, line-by-line, with some delimiter (separating character) in between them. This delimiter by default is a <kbd>tab</kbd> character, but we can set it to other things too, including a *newline* character. To demonstrate this, let's make a fasta-formatted sequence file from our genes in the previous lesson. 

>"Fasta" is a common format for holding sequence information. In it, each sequence entry takes up two lines: the first is the name of the sequence and needs to be preceded by a <kbd><b>></b></kbd> character; and the second line is the sequence. It looks like this:
>
>```
>>Seq_1
>ATGCGACC
>>Seq_2
>TCCGACTT
>```

To start, let's copy over our table that holds the gene IDs, lengths, and sequences (remember the **`.`** says to copy it to our current location and keep the same name):

```bash
cp ~/shell_intro/six_commands/genes_and_seqs.tsv .
```

This file holds the gene IDs in the first column and the sequences in the third:

```bash
head -n 1 genes_and_seqs.tsv
```

Let's get them into their own files. Note the use of **`-n +2`** in the **`tail`** command here. This takes everything in the file *except* the first line, which we don't want here because it is the header of the table:

```bash
cut -f 1 genes_and_seqs.tsv | tail -n +2 > ids.tmp
cut -f 3 genes_and_seqs.tsv | tail -n +2 > seqs.tmp

head ids.tmp
head seqs.tmp
```

We also need to add the <kbd><b>></b></kbd> character in front of our IDs though because that is part of the fasta format. We can do that with **`sed`** and using a special character that represents the start of every line (**`^`**):

```bash
sed 's/^/>/' ids.tmp > fasta_ids.tmp

head fasta_ids.tmp
```

This **`sed`** command is searching for the start of every line (**`^`**), and then adding in the <kbd><b>></b></kbd> character.

Now for the interleaving, to think about what's happening here, remember that **`paste`** is normally just sticking things together with a <kbd>tab</kbd> in between them:

```bash
paste fasta_ids.tmp seqs.tmp | head -n 2
```

But we can tell **`paste`** to use a different delimiter by providing it to the **`-d`** argument. Here is if we wanted the delimiter to be a dash:

```bash
paste -d "-" fasta_ids.tmp seqs.tmp | head -n 2
```


And we can also tell it to combine them with a newline character in between (which is represented by **`\n`**):

```bash
paste -d "\n" fasta_ids.tmp seqs.tmp | head -n 4
```

And that's our fasta-formatted file! So let's write it to a new file and get rid of the temporary files we made along the way:

```bash
paste -d "\n" fasta_ids.tmp seqs.tmp > genes.faa

head genes.faa

rm *.tmp
```

Now imagine we want to pull out all of the sequences that were annotated with that function we looked at before, epoxyqueuosine reductase, which we figured out had the KO identifier "K18979". We can get the gene IDs using **`grep`** like we did previously and then using **`cut`** to just keep the first column (note that we are providing the *relative* path to this file, starting from our current location):

```bash
grep "K18979" ../six_commands/gene_annotations.tsv | cut -f 1
```

And let's write them to a file:

```bash
grep "K18979" ../six_commands/gene_annotations.tsv | cut -f 1 > target_gene_ids.txt

ls
head target_gene_ids.txt
```

For pulling a few sequences out of a fasta file, **`grep`** can be very convenient. But remember the format of fasta is each entry takes two lines, and if we use **`grep`** with default settings to find a gene ID, we will only get the line with the gene ID:

```bash
grep "99" genes.faa
```

Fortunately, **`grep`** has a handy parameter that let's you pull out lines following your matched text also (in addition to just the line with the matched text), it's the **`-A`** parameter. So we can tell **`grep`** to pull out the line that matches *and* the following line like so: 

```bash
grep -A 1 "99" genes.faa
```

Cool! There's one more nuance we need to address though, and that is whether **`grep`** is looking for exact matches only or not. For example, trying to grab gene "9" does not do what we want:

```bash
grep -A 1 "9" genes.faa
```

It grabs everything that has a "9" in it. But we can tell **`grep`** to only take exact matches, meaning it needs to be the full word, if we provide the **`-w`** flag (for **w**ord). Here, the "word" can't have any whitespace inside it (spaces, tabs, and newline characters count as whitespace). We then also just need to add the leading <kbd><b>></b></kbd> character in front of the sequence ID we want:

```bash
grep -w -A 1 ">9" genes.faa
```

Great! Back to our target genes, we usually won't want to do that for every individual sequence we want (even though we only have 2 in our example here). So let's loop through our "target_gene_ids.txt" file! Here's just with **`echo`** like we did above, to see how we can add the <kbd><b>></b></kbd> in front of the variable:

```bash
for gene in $(cat target_gene_ids.txt)
do
  echo $gene
  echo ">$gene"  
done
```

>Note that since the **`>`** character is special at the command line (it redirects output), we need to put what we want to give to **`echo`** in quotes so that the whole thing goes to it (here ">$gene").

Now let's put the **`grep`** command we made above in the loop, and change what we're looking for to be the <kbd><b>></b></kbd> character followed by the **`$`** and variable name (just like it was provided to **`echo`**), and write the output to a new file:

```bash
for gene in $(cat target_gene_ids.txt)
do
  grep -w -A 1 ">$gene" genes.faa
done > target_genes.faa

ls
head target_genes.faa
```

## Summary
Even though loops can get much more complicated as needed, getting these foundational skills down is all that's needed to start harnessing their awesome power 🙂

---

<a href="shell-six-glorious-commands-04.html" style="float: left"><b>Previous:</b> Six glorious commands</a>

<a href="index.html" style="float: right"><b>Back to:</b> Shell crash course home</a><br>

<div style="text-align: center">
	<a href="../index.html">ANGUS Home</a>
</div>
