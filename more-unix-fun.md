# Some more practical use of Unix

In many of our lessons we have been working with RNASeq data from *Saccharomyces cerevisiae*. In the [transcriptome assembly lesson](transcriptome-assembly.md), which we didn't run through this year (yet, at least), we would have assembled these RNA reads into contigs using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki), a popular transcriptome assembler. Here we will compare the results of that assembly to the [reference](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/) RNA sequences from NCBI while getting some more practice at the command line! 

---

## Objectives

* Explore some more bash commands in practical usage examples (e.g. **`grep`**, **`sed`**, **`cut`**, **`comm`**, **`cat`**)
* Practice stringing together multiple commands with pipes (**`|`**) step-by-step
* Use blast to compare an assembled transcriptome to its reference transcriptome
* Run a "remote" blast from the command line

## Accessing our JetStream instances
You should still have your jetstream instance running, you can follow the instructions [here](jetstream/boot.html) to log in to JetStream and find your instance. Then `ssh` into it following the instructions [here](jetstream/boot.html#ssh-secure-login).

## Setting up a new conda environment
We will be using blast here, and need one more package for being able to BLAST remotely (send sequences to their servers from the command line). So let's make a new environment, install those, and enter our new environment:

```bash
conda create -y -n blast_env -c conda-forge -c bioconda blast gnutls
conda activate blast_env
```

## Setting up our working directory
Let's set up a new working directory and download the files we'll be working with:

```bash
cd ~
mkdir more-unix-fun/
cd more-unix-fun/

curl -L https://ndownloader.figshare.com/files/16200176 -o our-transcriptome-assembly.fa
curl -L https://ndownloader.figshare.com/files/16218335 -o ref-transcripts.fa
```

Out of curiousity, let's see how many contigs we assembled vs how many open-reading frames there are in the reference we have:

```bash
grep -c ">" our-transcriptome-assembly.fa
grep -c ">" ref-transcripts.fa
```

<blockquote>
<center><b>QUICK QUESTION!</b></center>

What might be some of the reasons for the large difference between the number of transcripts we assembled and the number of reference transcripts? 

<div class="toggle-header closed">
    <strong>Possible Causes</strong>
</div>

<div class="toggle-content docutils container" style="width:100%">

<ul>
  <li>not all transcripts may have been expressed at the time of sampling</li>
  <li>even expressed transcripts at low relative abundance may not have been amplified and sequenced</li>
  <li>not all transcripts that were expressed and sequenced may have assembled successfully</li>
</ul>
(This is not an exhaustive list of possibilities.)
</div>
</blockquote>

Even though we assembled fewer transcripts than the reference holds, let's use BLAST and the command line to try to find out if all of what we *did* assemble can actually be found in our reference!

## Making a blast database of our reference
Here we are going to make our blast database using our original reference fasta. Then we are going to blast the assembled transcripts against it.

```bash
makeblastdb -dbtype nucl -in ref-transcripts.fa -out ref-transcripts.blastdb
```

>**CODE BREAKDOWN**
>
> - **`makeblastdb`** - this is our command to make a blast database out of our reference fasta file
>   - **`-dbtype nucl`** - here we are specifying there are nucleotides
>   - **`-in`** - specifying the input file
>   - **`-out`** - naming the prefix of the output database files

## Blasting our assembly against the reference
Now we're going to try to align our assembled contigs against the reference. There are a lot of options for `blastn` that can be seen by running `blastn -help`, and there's an explanation of what we've done here in the following code breakdown block. 

```bash
blastn -query our-transcriptome-assembly.fa -db ref-transcripts.blastdb \
       -max_target_seqs 1 -max_hsps 1 \
       -out assembly-to-ref-blastout.tsv \
       -outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore"
```

>**CODE BREAKDOWN** 
>
> - **`blastn`** - this is our base command, to blast nucleotide queries against a nucleotide database (that's what a blast**n** is)
>   - **`-query`** - specifies our input fasta
>   - **`-db`** - specifies our database we just made
>   - **`-max_target_seqs 1`** - sets the maximum number of reference sequences returned to 1
>   - **`-max_hsps 1`** - sets maximum **h**ighest-**s**coring **p**airs returned (since BLAST is a local alignment, you can get more than one alignment returned between one query and one reference)
>   - **`-out`** - specifies name of output file
>   - **`-outfmt`** - specifies the output format

Now let's take a look (remember we can exit `less` by pressing the `q` key): 

```bash
column -t assembly-to-ref-blastout.tsv | less -NS
```

> **Code breakdown:**  
> * `column` – the `column` command with the `-t` flag will attempt to keep columns together based on the tabs. This can be convenient sometimes. 
> * `less` – since there are a lot of lines here, we are piping (`|`) the output into `less`
>   * the `-NS` flags we are providing add line numbers and prevent lines from softwrapping at the edge of the terminal.

Let's see how many of our assembled contigs successfully aligned to the orf reference fasta.

```bash
wc -l assembly-to-ref-blastout.tsv
```

Remember that `wc -l` tells us how many lines there are in the file. 3,262 of our contigs successfully aligned to the orf reference fasta. But we had a total of 3,323 contigs from our assembly as we saw with `grep -c ">" our-transcriptome-assembly.fa`, so some didn't successfully align to the reference.

Here is one way we can do a quick calculation at the command line:

```bash
echo "3323-3262" | bc
```

**This tells us we have 61 contigs from our assembly that did *not* successfully align to the reference. Let's find out what they are!**

## What didn't align?

Here are the steps we'll take to get the sequences that didn't align, and then try to find out what they are: 

1. **get all the names of the contigs from our assembly**  
2. **get all the names of the contigs from our assembly that were reported as "hits" in the blast output**  
3. **compare these to figure out which contigs from our assembly are not in the list of contigs reported as successfully aligning to the reference**  
4. **use this new list of contig names to pull out their sequences in fasta format**  
5. **blast the sequences against NCBI's nr database "remotely" (from the command line, sending our sequences to the NCBI servers)**  

---

### 1. Get all the names of the contigs from our assembly into a file

```bash
grep ">" our-transcriptome-assembly.fa | tr -d ">" | cut -f1 -d " " | sort > all-assembly-contig-IDs.txt
```

> **Code breakdown:** 
> * `grep` - Just like we used `grep` above to count how many sequences we had by providing the `-c` flag, here we are leaving off the `-c` flag so that it will pull out the lines that match. Since fasta format dictates the only place the ">" character appears is in front of sequence headers, we can pull out all head lines with that.
> * `tr` - We then "pipe" (`|`) that output into a command called `tr`, which is useful for manipulating indvidual characters. Here we are using to delete all of the ">" in the file (so our names are just names without the ">" in front of them, like they are in the blast output file).
> * `cut` - We then pipe that into the `cut` command, which is good for manipulating columns. Here, we're going to use it to cut the first column (`-f1`) setting the delimiter to a space (`-d " "`)
>   * this is because blast automatically cut the trailing space off of our sequence headers, and we want to match their format
>   * you can see this by running `head assembly_to_orf_coding_blastnout_with_header.tsv` and `head yeast-transcriptome-assembly.fa`
> * `sort` - We use sort here because the command we're going to use later to compare our two lists of headers needs them to be sorted in the same fashion, and running this on both will ensure that
> * `>` – We then redirect the output to a new file called "all-assembly-contigs.txt"

### 2. Get all the names of the contigs from our assembly that were reported as "hits" in the blast output

```bash
cut -f1 assembly-to-ref-blastout.tsv | sort > all-assembly-contig-IDs-with-hits.txt
```

> **Code breakdown:** 
> * `cut` – Here we are using `cut` to cut the first column (the sequence name) from the blast output file.
>   * Note that here we didn't need to provide the `-d` flag to set the delimiter. This is because by default `cut` sets the delimiter to tabs.
> * `sort` – We then pipe the output from `cut` into `sort`
>   * as mentioned above, to compare the two lists of names the way we are going to below requires that the lists be sorted in the same fashion, and running `sort` on both of them ensures that
> * `>` – We then redirect the output to a new file called "all-assembly-contig-hits.txt"

### 3. Compare these to figure out which contigs from our assembly are not in the list of contigs reported as successfully aligning to the reference

We can use the `comm` command (compare with an extra "m" for some reason...) to quickly find out which sequences we assembled that didn't successfully align to the reference transcriptome.

```bash
comm -23 all-assembly-contig-IDs.txt all-assembly-contig-IDs-with-hits.txt > all-assembly-contig-IDs-that-did-not-hit-ref.txt
```

> **Code breakdown:** 
> * `comm` – this command compares the lines in two files and by default returns 3 columns:
>   1. The lines unique to file 1 (first file given)
>   2. The lines unique to file 2 (second file given)
>   3. The lines common to both files
>     * Looking at the manual for it with `man comm` we can see you can suppress columns by provide them as arguments
>     * So by providing the flags `-23` we are saying to keep those, all we want are the lines (contig names) in file 1 that are not in file 2. This gives us all the contig names that we assembled but that did not successfully align to the reference.
> * `>` – We then redirect the output to a new file called "all-assembly-contigs-that-did-not-hit-ref.txt"

And this should hold the 65 contig names that we're looking for (sanity check):

```bash
wc -l all-assembly-contig-IDs-that-did-not-hit-ref.txt
```

Stellar 🙂

### 4. Use this new list of contig names to pull out their sequences in fasta format

Now let's use `grep` to pull out the sequences! 

> **NOTE:** This requires that the fasta file is in "true" fasta form – meaning that each entry (sequence) takes up two lines. Sometimes, fasta files are formatted with what are known as softwraps. There is a one-liner at the end of this page to convert them if that is the case with your data. 

If we wanted to pull one sequence, we could do it like this:

```bash
grep -A1 "TRINITY_DN501_c0_g1_i1" our-transcriptome-assembly.fa
```

```bash
for header in $(cat all-assembly-contig-IDs-that-did-not-hit-ref.txt)
do
  grep -A1 "$header" our-transcriptome-assembly.fa
done
```

> **Code breakdown:** 
> The special words of the loop are `for`, `in`, `do`, and `done`.
> * `for` – here we are specifying the variable name we want to use ("header"), this could be anything
> * `in` – here we are specifying what we are going to loop through, in this case it is every line of the "all-assembly-contigs-that-did-not-hit-ref.txt" file
>   * `$(cat ...)` – this is a special type of notation in shell. The operation within the parentheses here will be performed and the output of that replaces the whole thing. It would be the same as if we typed out each line in the file with a space in between them
> * `do` – indicates we are about to state the command(s) we want to be performed on each item we are looping through
>   * `grep` – just like with the individual example above, for each sequence name in our file we are pulling out the header line and the following line containing the sequence
> * `done` – tells the loop it is finished


This just printed to the screen for right now, but let's grab one and take it to [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and do a web blast to see what we get.

```
TRINITY_DN939_c0_g1_i1 len=900 path=[0:0-899]
TTCGACTATGATACTAATTCGCTGTGCACCCAGTGGGAGGATGTTATTTCGCTTACATGCCATATGCTGTATATTCTTGCGTTAACGGTTTCCGCGTGCGATCTCTCTTGTGTTCAGACGAGGCCCAATTGAGCACCATCCCCCTCGGGTAGTTTCCCGATCAAACTGGAAGATAGGCGTCTTTACTTACGCGCTCCTCTTGACCGAGACCCCCAATGCGCGATGTATCGAACCTTCACTAACCCTAGAAATTAGTGGTGGGAATCAGCGAAGTTACAATGTGGGGTTGGACCCAGGATGTTAGCCTGCAAGCTATACAATTCTCTTAGATTAGACGAGAACGGAGAATTTAACCCCTGCAGCATTGGAGGTATGGTCTTGGGCATACCCGATACATGCAACGCAGCTCGGGATGTTCATGGTAGCACCTAACTGTATGGCATAGTTATGCAGAAGTGCGCTGCTTAAGAGCGATACCCCATAAAGAACGATTTTGGTGGTATTGCCCAAAGATAATGTCCCACGTTATCATCTGGTCAACGATGAGGTGGGTTGTTTTGTGATTGTTTGAGATGCTGAGTGCTGTTTAATGCGGGACATAAGGAAGGATATTAGTAGGGAGAAACGCTTGATGCCGGAAATATCCTTGCCTGGTTAACTGCTCGAAGTTAATCTGCGACGCTCGCCCTCATTCGGATGCATCGAAGGGCTCCCCTGCAGTTGCAAAGTCTTTGTTCTGCGAACTCGTAAAGTCGTAATGCCGTTGGTGGACCGTGCTTGTTAGGGATATTAAATGTTTCCTGGCCTTTAAAGCTATTGGCACGGCGGTTTAGATGGGACACCCTATCTCGTTTTCTACTTGCGCTTCAAGCGTCCCAACGAAACGAAATTGCGGACCGG
```

Now let's write these all to a file instead of printing them to a screen. We're just adding the `>` redirector and a file name at the end of the loop to save the output.

```bash
for header in $(cat all-assembly-contig-IDs-that-did-not-hit-ref.txt)
do
  grep -A1 "$header" our-transcriptome-assembly.fa
done > contigs-not-in-ref.fa
```

### 5. Blast the sequences against NCBI's nr database "remotely" (from the command line, sending our sequences to the NCBI servers)

We can send our sequences to NCBI through the command line to blast them instead of doing it at their website. But in order to include taxonomic information with the output, we need to have another small blast database that holds that information.  

We can download and unpack that database with these commands:

```bash
curl -L https://ndownloader.figshare.com/files/16219079 -o taxdb.tar.gz
tar -xzvf taxdb.tar.gz
```

Now we have the taxonomy database too so we can add taxonomy to our blast output. Let's blast them!

```
nohup blastn -query contigs-not-in-ref.fa -remote -db nr \
       -max_target_seqs 1 -max_hsps 1 \
       -out contigs-not-in-ref-blast-to-nr.tsv \
       -outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore ssciname" &
```

> **CODE BREAKDOWN**
> 
> The new things added here different from the previous blast are:
> 
> - **`nohup`** – this tells the computer not to stop the process if we disconnect from the server (very handy)
> - **`-db`** – here we are specifying "nr" which is for NCBI's "non-redundant" database (nucleotide in this case, which it knows because we used `blastn`
> - **`-remote`** – this tells the computer to send the BLAST job to the NCBI servers
> - **`&`** - this at the end tells the computer to run the process in the background so we can still work in our terminal. We can check on the job with the command `jobs`. This will finish when it does and it won't be interrupted disconnect or sign off before then. 

This may take a few minutes (the amount of traffic the blast servers deal with fluctuates). So we can download a results file as follows and move on while that continues in the background:

```bash
curl -L https://ndownloader.figshare.com/files/16219238 -o contigs-not-in-ref-blast-to-nr-DL.tsv
```

## What did we get?
The 9th column contains the taxonomy information, so we can use the `cut` command to look at just that.

```bash
cut -f 9 contigs-not-in-ref-blast-to-nr-DL.tsv | sort | uniq -c
```

> **CODE BREAKDOWN** 
> 
> - **`cut`** – this crops out the specified columns, here we're cutting column 9 (taxonomy column, labeled as "ssciname" by blast)
> - **`sort`** – sorting is required in order for the following `uniq` command to function properly
> - **`uniq`** – this removes all duplicates, leaving only single copies of the originals
>   - the **`-c`** flag provided to `uniq` also tells us how many of each were found

```
      1 Bacillus subtilis subsp. subtilis str. NCIB 3610
      2 Cloning vector pBAD-FLP
      1 Methanocaldococcus jannaschii
     33 Saccharomyces cerevisiae
      1 Yeast expression vector pSCUDT2MFE6H
     23 synthetic construct
```

We see there are quite a few synthetic constructs in here, we might want to remove them and look further into these other sequences if this were a real project. And there are 33 *S. cerevisiae* sequences that didn't successfully align to our reference cDNA file.


<blockquote>
<center><b>QUICK QUESTION!</b></center>

Why might there by <i>S. cerevisiae</i> sequences assembled from our RNA that are not aligning to the reference cDNA file?

<div class="toggle-header closed">
    <strong>Possible Causes</strong>
</div>

<div class="toggle-content docutils container" style="width:100%">

It's possible some of these might be non-coding RNAs, and therefore aren't in the cDNA reference. There may (probably are) other reasons too.

</div>
</blockquote>

**Isn't Unix grand 🙂**
