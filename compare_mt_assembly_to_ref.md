# Comparing our de novo transcriptome assembly to the reference

In the [transcriptome assembly lesson](https://angus.readthedocs.io/en/2018/transcriptome-assembly.html), we assembled short RNA reads into contigs. Here we will compare this to the reference transcriptome we worked with in the [mapping and variant calling](https://angus.readthedocs.io/en/2018/mapping-variant-calling.html) lesson, and get some more practice at the command line and working with blast!  

---

Learning objectives:

	* Explore some more bash commands in a practical way (e.g. `grep`, `sed`, `cut`, `comm`, `cat`)
	* Practice stringing together multiple commands with pipes (`|`) step-by-step
	* Use blast to compare two transcriptomes to find sequences that are different
	* Run a "remote" blast from the command line

## Boot up a Jetstream
[Boot an m1.medium Jetstream instance](https://angus.readthedocs.io/en/2018/jetstream/boot.html) and log in.

## Software needed
We will be using blast here, so if you've run through the [using blast at the command line](https://angus.readthedocs.io/en/2018/running-command-line-blast.html) lesson, you're already good to go. You can check that you have blast by entering:

```bash
makeblastdb -h
```

If that returns some help information, blast is already installed. If it returns a message with "command not found", then you can install blast with:

```bash
conda install -y blast
```

## Setting up our working directory

Let's set up a new working directory to compare our assembly to the orf reference and copy over the reference we worked with and the result from our de novo transcriptome assembly (they're both small enough files that copying isn't a big deal).

> NOTE: ".fa", ".fasta", and ".fna" are all common extensions for a nucleotide fasta file. Amino acid fasta files typically have ".faa".

```bash
cd ~
mkdir comparing_assembly_to_ref
cd comparing_assembly_to_ref/
cp ../mapping/orf_coding.fasta .
cp ../assembly/yeast-transcriptome-assembly.fa .
```

> If either of those returned an error message (file not found), you can download the reference by running `curl -O https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz` followed by `gunzip orf_coding.fasta.gz` to download and unzip the reference, and/or `curl -O https://ndownloader.figshare.com/files/12323351` to download a result of the transcriptome assembly (this one is not gzipped).

Out of curiousity, let's see how many contigs we assembled vs how many are in the orf reference file we have?

```bash
grep -c ">" yeast-transcriptome-assembly.fa
grep -c ">" orf_coding.fasta
```

> **Code breakdown:** This is a common, quick way to check how many sequences are in a fasta file. 
> * `grep` by default will pull out all lines containing the pattern we give it, from the file we provide. Here we are adding the `-c` flag telling it to just count the lines for us, and the pattern we are searching for is just the ">" character. 
>   * This is because if you remember in true fasta format those should only be found once for each sequence directly in front of the header line.
>   * You can peak at one of the files to remind yourself of the format with `head yeast-transcriptome-assembly.fa`. 

**How can we check if the contigs we assembled are all in the orf reference? Or find which ones aren't and might be spurious? Blast is one way.**

### Making a blast database of our reference
Here we are going to make our blast database using our original reference fasta. Then we are going to blast the assembled transcripts against it.

```bash
makeblastdb -dbtype nucl -in orf_coding.fasta -out orf_coding_blastdb
```

> Note that we are providing a name to the `-out` flag. This is going to be the "basename" of the files that are created. Running `ls orf_coding_blastdb*` after it's done will show them.

### Blasting our assembly against the reference
Now we're going to try to align our assembled contigs against the reference. You can see all the options for `blastn` by running `blastn -help`, and there's an explanation of what we've done here in the following code breakdown block. 

> **NOTE:** Ignore the "\\" at the end of lines within code blocks like below. They are there to tell the terminal to ignore the newline characters that were added here so that these can be copied and pasted and still run properly. Feel free to copy and paste. Typing out is good practice too sometimes, but thinking about what the code is doing is probably more useful here.

```bash
blastn -query yeast-transcriptome-assembly.fa -db orf_coding_blastdb \
-max_target_seqs 1 -max_hsps 1 \
-out assembly_to_orf_coding_blastnout.tsv -outfmt "6 qseqid qlen sseqid slen pident \
length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore score"
```

> **Code breakdown:** We've added a couple new ones here, namely `-max_target_seqs 1` and `-max_hsps 1`. Both of these are here to ensure we only get one entry back for each of our query sequences, in theory, the "best hit". The first tells the program to only record one "target" sequence. In blast terminology, the "query" is what you're blasting (our assembly contigs), and the "target" is what you're blasting against (our reference fasta in this case). The second, dealing with "hsps" is for "highest scoring pairs". Since blast is a "local" alignment (**B**asic **L**ocal **A**lignment **S**earch **T**ool), you can get more than one alignment between the same two sequences in different places.

Sometimes it's useful to add headers to tabular blast output files right at the command line. This can be done in `nano` (be sure to put tabs between the column names if you do it this way), but sometimes files are large and `nano` can be sluggish too. So here's a little command-line detour showing another way.

```bash
cat <(printf "qseqid\tqlen\tsseqid\tslen\tpident\tlength\tqcovhsp\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tscore\n") \
assembly_to_orf_coding_blastnout.tsv > assembly_to_orf_coding_blastnout_with_header.tsv
```

There's quite a bit going on here in these two lines so here's a full breakdown. 

> **Code breakdown:**  
> In the previous code block, we're using the `cat` command to stick our header on top of our tabular blast output file.
> * the first argument we're passing to `cat` is using a new syntax "<(...)"
>   * this will evaulate what's inside the parentheses, then pass it the `cat` command
>   * within these parentheses we are using the `printf` command because it nicely translates the "\t" characters into tabs
>     * copy and paste just that part of the command to see what it does: `printf "qseqid\tqlen\tsseqid\tslen\tpident\tlength\tqcovhsp\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tscore\n"`
>     * we also add a newline character ("\n") at the end of the `printf` command so the blast output file starts at line 2
> * the second argument to `cat` is our tabular blast output file
> * last we redirect (">") the output to a new file rather than printing it to the screen

And don't worry if this still just looks like a jumbled mess. It's not so much about remembering commands and arguments or the oddities of using `sed`. Spending time understanding the fundamental rules of running commands and understanding the code is very different than memorizing it. Knowing the base rules and the capabilities is enough to help us google the rest :)  
Keep in mind that we didn't overwrite the initial file (with no header), we just made a new one for our viewing.

Now let's take a look! Remember you can exit `less` by pressing the "q" key.

```bash
column -t assembly_to_orf_coding_blastnout_with_header.tsv | less -NS
```

> **Code breakdown:**  
> * `column` – the `column` command with the `-t` flag will attempt to keep columns together based on the tabs. This can be convenient sometimes. 
> * `less` – since there are a lot of lines here, we are piping (`|`) the output into `less`
>   * the `-NS` flags we are providing add line numbers and prevent lines from softwrapping at the edge of the terminal.

Let's see how many of our assembled contigs successfully aligned to the orf reference fasta.

```bash
wc -l assembly_to_orf_coding_blastnout.tsv 
```

`wc -l` tells us how many lines there are in the file. 3,258 of our contigs successfully aligned to the orf reference fasta. But we had a total of 3,323 contigs from our assembly as we saw with `grep -c ">" yeast-transcriptome-assembly.fa`, so some didn't successfully align to the reference.

This is one way you can do a quick calculation at the command line:

```bash
echo "3323-3258" | bc # look at the man page (`man bc`) or google "bc command examples in unix" if interested
```

**This tells us we have 65 contigs from our assembly that did not successfully align to the reference. Let's find out what they are!**

---

## Here's the path we'll take to get the sequences that didn't align and try to find out what they are: 
**1. get all the names of the contigs from our assembly**  
**2. get all the names of the contigs from our assembly that were reported as "hits" in the blast output**  
**3. compare these to figure out which contigs from our assembly are not in the list of contigs reported as successfully aligning to the reference**  
**4. use this new list of contig names to pull out their sequences in fasta format**  
**5. blast the sequences against NCBI's nr database "remotely" (from the command line, sending our sequences to the NCBI servers)**  

---

**1. Get all the names of the contigs from our assembly into a file**

```bash
grep ">" yeast-transcriptome-assembly.fa | tr -d ">" | cut -f1 -d " " | sort > all_assembly_contigs.txt
```

> **Code breakdown:** 
> * `grep` - Just like we used `grep` above to count how many sequences we had by providing the `-c` flag, here we are leaving off the `-c` flag so that it will pull out the lines that match. Since fasta format dictates the only place the ">" character appears is in front of sequence headers, we can pull out all head lines with that.
> * `tr` - We then "pipe" (`|`) that output into a command called `tr`, which is useful for manipulating indvidual characters. Here we are using to delete all of the ">" in the file (so our names are just names without the ">" in front of them, like they are in the blast output file).
> * `cut` - We then pipe that into the `cut` command, which is good for manipulating columns. Here, we're going to use it to cut the first column (`-f1`) setting the delimiter to a space (`-d " "`)
>   * this is because blast automatically cut the trailing space off of our sequence headers, and we want to match their format
>   * you can see this by running `head assembly_to_orf_coding_blastnout_with_header.tsv` and `head yeast-transcriptome-assembly.fa`
> * `sort` - We use sort here because the command we're going to use later to compare our two lists of headers needs them to be sorted in the same fashion, and running this on both will ensure that
> * `>` – We then redirect the output to a new file called "all_assembly_contigs.txt"

**2. Get all the names of the contigs from our assembly that were reported as "hits" in the blast output**  

```bash
cut -f1 assembly_to_orf_coding_blastnout.tsv | sort > all_assembly_contig_hits.txt
```

> **Code breakdown:** 
> * `cut` – Here we are using `cut` to cut the first column (the sequence name) from the blast output file.
>   * Note that here we didn't need to provide the `-d` flag to set the delimiter. This is because by default `cut` sets the delimiter to tabs.
> * `sort` – We then pipe the output from `cut` into `sort`
>   * as mentioned above, to compare the two lists of names the way we are going to below requires that the lists be sorted in the same fashion, and running `sort` on both of them ensures that
> * `>` – We then redirect the output to a new file called "all_assembly_contig_hits.txt"

**3. Compare these to figure out which contigs from our assembly are not in the list of contigs reported as successfully aligning to the reference**  

We can use the `comm` command (compare with an extra "m" for some reason...) to quickly find out which sequences we assembled that didn't successfully align to the reference transcriptome.

```bash
comm -23 all_assembly_contigs.txt all_assembly_contig_hits.txt > all_assembly_contigs_that_did_not_hit_ref.txt
```

> **Code breakdown:** 
> * `comm` – this command compares the lines in two files and by default returns 3 columns:
>   1. The lines unique to file 1 (first file given)
>   2. The lines unique to file 2 (second file given)
>   3. The lines common to both files
>     * Looking at the manual for it with `man comm` we can see you can suppress columns by provide them as arguments
>     * So by providing the flags `-23` we are saying to keep those, all we want are the lines (contig names) in file 1 that are not in file 2. This gives us all the contig names that we assembled but that did not successfully align to the reference.
> * `>` – We then redirect the output to a new file called "all_assembly_contigs_that_did_not_hit_ref.txt"

And this should hold the 65 contig names that we're looking for (sanity check):

```bash
wc -l all_assembly_contigs_that_did_not_hit_ref.txt
```

Good news? It was when putting this together :)

**4. Use this new list of contig names to pull out their sequences in fasta format**  

Now let's use `grep` to pull out the sequences! 

> NOTE: This requires that the fasta file is in "true" fasta form – meaning that each entry (sequence) takes up two lines: a header on one line preceded by a ">", and a sequence on one line. Both of our fasta files here are formatted this way. You can check by using `head`. If you're sequences are not formatted like this, and have sequences spread across multiple lines in order to keep column width narrow, there is a one-liner at the bottom of this page to convert them.

If we wanted to pull one sequence, we could do it like this:

```bash
grep -A1 "TRINITY_DN501_c0_g1_i1" yeast-transcriptome-assembly.fa
```

Here we're using `grep` to grab the line that contains the contig name we want, from the file we provide at the end. By default this would just grab the header line, but we also want the sequence. The `-A1` flag tells `grep` to also give us back the line following the matching name – which gives us the header and the sequence for each entry we want. 

Now we're going to loop through the file of names we want, and do this for each one:

```bash
for header in $(cat all_assembly_contigs_that_did_not_hit_ref.txt)
do
grep -A1 "$header" yeast-transcriptome-assembly.fa
done
```

> **Code breakdown:** 
> The special words of the loop are `for`, `in`, `do`, and `done`.
> * `for` – here we are specifying the variable name we want to use ("header"), this could be anything
> * `in` – here we are specifying what we are going to loop through, in this case it is every line of the "all_assembly_contigs_that_did_not_hit_ref.txt" file
>   * `$(cat ...)` – this is a special type of notation in shell. The operation within the parentheses here will be performed and the output of that replaces the whole thing. It would be the same as if we typed out each line in the file with a space in between them
> * `do` – indicates we are about to state the command(s) we want to be performed on each item we are looping through
>   * `grep` – just like with the individual example above, for each sequence name in our file we are pulling out the header line and the following line containing the sequence
> * `done` – tells the loop it is finished


This just printed to the screen for right now, but grab one and take it to [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and do a web blast to see what you get!

Maybe try this one...

>TRINITY_DN939_c0_g1_i1 len=900 path=[0:0-899]
TTCGACTATGATACTAATTCGCTGTGCACCCAGTGGGAGGATGTTATTTCGCTTACATGCCATATGCTGTATATTCTTGCGTTAACGGTTTCCGCGTGCGATCTCTCTTGTGTTCAGACGAGGCCCAATTGAGCACCATCCCCCTCGGGTAGTTTCCCGATCAAACTGGAAGATAGGCGTCTTTACTTACGCGCTCCTCTTGACCGAGACCCCCAATGCGCGATGTATCGAACCTTCACTAACCCTAGAAATTAGTGGTGGGAATCAGCGAAGTTACAATGTGGGGTTGGACCCAGGATGTTAGCCTGCAAGCTATACAATTCTCTTAGATTAGACGAGAACGGAGAATTTAACCCCTGCAGCATTGGAGGTATGGTCTTGGGCATACCCGATACATGCAACGCAGCTCGGGATGTTCATGGTAGCACCTAACTGTATGGCATAGTTATGCAGAAGTGCGCTGCTTAAGAGCGATACCCCATAAAGAACGATTTTGGTGGTATTGCCCAAAGATAATGTCCCACGTTATCATCTGGTCAACGATGAGGTGGGTTGTTTTGTGATTGTTTGAGATGCTGAGTGCTGTTTAATGCGGGACATAAGGAAGGATATTAGTAGGGAGAAACGCTTGATGCCGGAAATATCCTTGCCTGGTTAACTGCTCGAAGTTAATCTGCGACGCTCGCCCTCATTCGGATGCATCGAAGGGCTCCCCTGCAGTTGCAAAGTCTTTGTTCTGCGAACTCGTAAAGTCGTAATGCCGTTGGTGGACCGTGCTTGTTAGGGATATTAAATGTTTCCTGGCCTTTAAAGCTATTGGCACGGCGGTTTAGATGGGACACCCTATCTCGTTTTCTACTTGCGCTTCAAGCGTCCCAACGAAACGAAATTGCGGACCGG

**What is the source of this contig we assembled?**

Now let's write these all to a file instead of printing them to a screen. We're just adding the `>` redirector and a file name at the end of the loop to save the output.

```bash
for header in $(cat all_assembly_contigs_that_did_not_hit_ref.txt)
do
grep -A1 "$header" yeast-transcriptome-assembly.fa
done > contigs_not_in_ref.fa
```

**5. Blast the sequences against NCBI's nr database "remotely" (from the command line, sending our sequences to the NCBI servers)**  
We can send our sequences to NCBI through the command line to blast them, we don't need their big sequence database for that. But in order to include taxonomic information with the output, we need to have another small blast database that holds that information.  

We can get this taxonomy database very easily using a script that came with the blast installation we did. 

```bash
update_blastdb.pl taxdb
```

This downloaded with 2 files: the tar.gz and an md5 file. This md5 file allows us to make sure the file wasn't corrupted when we downloaded it. You can check that like this:

```bash
md5sum -c taxdb.tar.gz.md5 # should print out file name and "OK"
```

Now unzipping and untarring:

```bash
tar -xz -f taxdb.tar.gz
```

> **Code breakdown:** 
> * `tar` – This is a program that can uncompress directories with multiple folders.
>   * `-xz -f`– the "x" flag tells it to decompress the target object ("taxdb.tar.gz" here), this is splitting it into multiple files; the "z" flag tells it to unzip things after uncompressing; and the "f" is for the input file

Now we have the taxonomy database too so we can add taxonomy to our blast output!
Let's blast them all!

```bash
nohup blastn -db nr -query contigs_not_in_ref.fa -remote -max_target_seqs 1 -max_hsps 1 \
-out contigs_not_in_ref_to_nr_blastout.tsv \
-outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score taxid ssciname" &
```

> **Code breakdown:** 
> * `nohup` – This tells the computer not to stop the process if you disconnect from the server. (very handy)
> * `blastn` – This is our base command here, telling the computer we're using the `blastn` program
>   * `-db` – This is specifying which database we want to use.
>   * `-query` – Our input sequences, here we are providing the contigs we assembled that did not successfully align to the reference.
>   * `-remote` – This tells the computer to send the blast job to the NCBI servers and use their database.
>   * `-max_target_seqs` and `-max_hsps` are the same as described above in our first blast.
>   * `-out` specifies the output file, and `-outfmt` specifies how we want it – tabular and which columns (see `blastn -help`)

This may take a few minutes (the amount of traffic the blast servers deal with fluctuates. If you'd like to just move forward, you can download a results file as follows and move on while that continues in the background:

```bash
curl -L https://ndownloader.figshare.com/files/12323003 -o contigs_not_in_ref_to_nr_blastout_dl.tsv
```

> **Code breakdown:** 
> * `curl` – "Client for URLs", this command downloads files from the web.
>   * the `-L` tells the program to follow links if needed to get the file
>   * the `-o` flag is where you can specify what you want the file to be saved as
>     * note the "_dl" added to the downloaded file here. That will be used in the following command, so modify if required to fit the remote blast results. 

**What did we get?**
The 16th column contains the taxonomy information, so we can use the `cut` command to look at just that.

```bash
cut -f16 contigs_not_in_ref_to_nr_blastout_dl.tsv | sort | uniq -c
```

> **Code breakdown:** 
> * `cut` – this crops out the specified columns, here we're cutting column 16 (taxonomy column)
> * `sort` – sorting is required in order for the following `uniq` command to function properly
> * `uniq` – this removes all duplicates, leaving on single copies of the originals
>   * the `-c` flag provided to `uniq` also tells us how many of each were found


```
      1 Bacillus subtilis
      3 Methanocaldococcus jannaschii
     39 Saccharomyces cerevisiae
      1 Saccharomyces cerevisiae YJM1592
     21 synthetic construct
```

And we see there are mostly the yeast we expect in here (*Saccharomyces*), but there are also quite a few "synthetic constructs". 

**Using the skills we saw above, see if you can get the sequence that blasted to *Bacillus subtilis*. Then run a web blast again to see what else it has closest hits too – remember that here we only kept the "best hit".**

