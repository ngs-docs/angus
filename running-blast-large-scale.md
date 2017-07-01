# Running large and long command line jobs - using shmlast!

Our goal for this tutorial is for you become more familiar with running longer programs
on the command line. You'll be introduced to
[shmlast](http://joss.theoj.org/papers/3cde54de7dfbcada7c0fc04f569b36c7), which is implements
an algorithm for discovering potential orthologs between an RNA-seq assembly and a protein database.

## Installing shmlast

Install base packages:

```
sudo apt-get -y update && \
sudo apt-get install -y python3.5-dev python3.5-venv make \
    libc6-dev g++ zlib1g-dev last-align parallel
```

Then create a Python environment with `virtualenv`, which will isolate your python packages:
```
python3.5 -m venv ~/py3
. ~/py3/bin/activate
pip install -U pip
```

And now install shmlast 1.2:
```
pip install shmlast>=1.2
```

This downloads and installs the latest version of shmlast.

## Download some data

Next we need some data! Here we're going to grab one of the three
mouse RNA data sets,

```
curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.1.rna.fna.gz
```

and all 8 of the cow protein data sets.

```
for i in 1 2 3 4 5 6 7 8
do
   curl -O ftp://ftp.ncbi.nih.gov/refseq/B_taurus/mRNA_Prot/cow.$i.protein.faa.gz
done
```

shmlast wants one query database (here, we'll use mouse) and one
database to be searched (here, cow) - but first we have to combine
all of the databases into one:

```
gunzip -c cow.*.faa.gz > cow.faa
```

## Run shmlast!

Now run shmlast:
```
shmlast crbl -q mouse.1.rna.fna.gz -d cow.faa --n_threads=6
```
this will take 16 minutes (!!) and produce some large files.

## Digression: What is shmlast doing?

shmlast is going to compute putative orthologs between mouse
transcripts and cow proteins.  Orthologs are genes that duplicated
from speciation (i.e. are the "same" gene in cow and mouse) and
are presumed to have the same function, although that is a computational
inference that needs to be treated with care.

As we'll see on Friday, the computation of orthologs (or homologs more
generally) is a core step in annotating genomes and transcriptomes.
The reason is that you don't automatically get gene assignments when
you build a new genome or transcriptome - you just get unidentified
DNA or RNA sequence! And then you have to name each transcript or gene,
and generally people want that name to be the same across species.
And that involves computing orthologs.

One of the most common ways to compute ortholog assignments is to use
reciprocal-best-hit BLAST, in which you use BLAST to find the two
sequences that match each other best in the database.  However,
reciprocal best hit has a few problems in the face of complicated
evolutionary scenarios or deep RNA sequencing; from the supp. material
of [Aubry et al., 2014](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004365),

```
[... reciprocal best hit ] perform[s] well when sequences are present
as single copy genes in the datasets being compared and perform[s] less
well when trying to distinguish highly similar sequence groups (such
as multi-copy genes) from each other.  The issue of multiple-copy
genes and near identical gene-groups is particularly relevant for the
analysis of transcriptome data. It is to be expected that following de
novo assemblies of RNAseq data, most gene loci will be represented by
multiple assembled transcript variants.
```

basically this is saying that in many realistic scenarios (most especially
multi-copy genes, but also multiple isoforms) reciprocal best hit is
too conservative and will ignore real orthologs.

So Aubry et al. invent *conditional* reciprocal best hit, which tries
to find close homolog *groupings* that can deal with multi-copy genes.
shmlast is a reimplementation of that, done by our very own Camille Scott.

### Why does shmlast take "so long"?

Well, we're calculating all pairwise matches between 36,000 mouse transcripts
and 64,000 cow proteins!  So frankly it's amazing it works so fast in the
first place!

----

## Looking at the output

Like most bioinformatics software, shmlast produces *a lot* of output.
How do we explore the results?

The main output is `mouse.1.rna.fna.gz.x.cow.faa.crbl.csv`, which is a
Comma-Separated Value file that you can load into any spreadsheet
program.  If we look at the file by typing `head
mouse.1.rna.fna.gz.x.cow.faa.crbl.csv` we should see something like this:

```
E,EG2,E_scaled,ID,bitscore,q_aln_len,q_frame,q_len,q_name,q_start,q_strand,s_aln_len,s_len,s_name,s_start,s_strand,score
6.6e-24,9.8e-16,23.18045606445813,641897,109.65804469295703,89,1,390,"ref|NM_001013372.2| Mus musculus neural regeneration protein (Nrp), mRNA",64,+,89,389,ref|XP_005212262.1| PREDICTED: DNA oxidative demethylase ALKBH1 isoform X1 [Bos taurus],0,+,241.0
5.4e-194,4.4e-165,193.26760624017703,719314,605.7589445367834,313,0,331,"ref|NM_207235.1| Mus musculus olfactory receptor 358 (Olfr358), mRNA",0,+,313,313,ref|XP_607965.3| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1365.0
2.8e-188,5e-160,187.5528419686578,423289,588.9868500580775,307,0,323,"ref|NM_146368.1| Mus musculus olfactory receptor 361 (Olfr361), mRNA",0,+,307,313,ref|XP_607965.3| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1327.0
6.6e-183,5.6e-155,182.18045606445813,725159,572.2147555793716,307,0,318,"ref|NM_146622.1| Mus musculus olfactory receptor 360 (Olfr360), mRNA",0,+,307,313,ref|XP_607965.3| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1289.0
5.4e-194,4.4e-165,193.26760624017703,719315,605.7589445367834,313,0,331,"ref|NM_207235.1| Mus musculus olfactory receptor 358 (Olfr358), mRNA",0,+,313,313,ref|XP_002691614.1| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1365.0
2.8e-188,5e-160,187.5528419686578,423290,588.9868500580775,307,0,323,"ref|NM_146368.1| Mus musculus olfactory receptor 361 (Olfr361), mRNA",0,+,307,313,ref|XP_002691614.1| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1327.0
6.6e-183,5.6e-155,182.18045606445813,725160,572.2147555793716,307,0,318,"ref|NM_146622.1| Mus musculus olfactory receptor 360 (Olfr360), mRNA",0,+,307,313,ref|XP_002691614.1| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1289.0
4.8e-183,5.6e-155,182.3187587626244,373474,572.2147555793716,266,0,310,"ref|XR_001782298.1| PREDICTED: Mus musculus predicted gene 4786 (Gm4786), misc_RNA",29,+,266,266,ref|NP_001035610.1| 60S ribosomal protein L7a [Bos taurus],0,+,1289.0
3.2e-153,3.1e-138,152.4948500216801,643504,516.6020212552417,246,1,659,"ref|NR_003628.1| Mus musculus predicted gene 5766 (Gm5766), non-coding RNA",357,+,246,266,ref|NP_001035610.1| 60S ribosomal protein L7a [Bos taurus],0,+,1163.0
```
(note that you can scroll to the right within the text to see all the output.)

Here the columns are helpfully labeled, but it's still kind of a mess
to look at - we'll look at in more detail in R, instead of using the
command line.  The key bits are the `q_name`_ and `s_name` column, which
tell you which *query* and which *subject* sequences match each other.

How big is this file? Big!  You can calculate how many lines are
present by using the `wc` command, which will tell you how many lines,
words, and characters are in the file:

```
  132901  2918423 39291181 mouse.1.rna.fna.gz.x.cow.faa.crbl.csv
```


(A precomputed version of this file is available [for download here]((https://github.com/ngs-docs/angus/raw/17a0ba3b1d915de90a5b8bd1fbc1027eba47baf8/_static/shmlast/mouse.1.rna.fna.gz.x.cow.faa.crbl.csv.gz).)

Next we're going to do some exploratory data analysis of this file
[using RStudio](visualizing-blast-scores-with-RStudio.html).

## Some points for discussion

* What evolutionary scenarios are there that complicate orthology and
  homology assignment?
