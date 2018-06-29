# Visualizing BLAST score distributions in RStudio

Learning objectives:

* learn the basics of the RStudio interface.

* explore plotting in RStudio.

* explore some characteristics of the data resulting from your BLAST search.

## Getting started

Connect to RStudio by setting your password (note, password will not
be visible on the screen):

```
sudo passwd $USER
```

figuring out your username:

```
echo My username is $USER
```

and finding YOUR RStudio server interface Web address:

```
echo http://$(hostname):8787/
```

Now go to that Web address in your Web browser, and log in with the username
and password from above.

## Enter some R commands

(Enter the below commands into RStudio, not the command line.)

Load the data you created with BLAST:

```
blast_out <- read.table('blast/mm-second.x.zebrafish.tsv', sep='\t')
```

If you run `View(blast_out)`, you'll see the same information as in
the previous section, *but* loaded into R for plotting and manipulation.

The only problem is that the column names are kind of opaque - what does V1
mean!?  To fix this, we can reset the column names like so, using the
information from [the BLAST outfmt documentation](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6):

```
colnames(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
View(blast_out)
```

`blast_out` is called a dataframe, which is a sort of R-ish version of a
spreadsheet with named columns.  `View` can be used to present it nicely,
and `head(blast_out)` can be used to look at just the first few rows.

Another useful command is `dim` which will tell you the DIMENSIONS of this
data frame:

```
dim(blast_out)
```

That's a big data frame! 14,720 rows (and 12 columns!)

Let's do some data visualization to get a handle on what our blast output looked like. First, let's look at the evalue:

```
hist(blast_out$evalue)
```

This is telling us that MOST of the values in the `evalue` column are
quite low.  What does this mean? How do we figure out what this is?

(You can also try plotting the distribution of `-log(blast_out$evalue)` - why
is this more informative?)

So these are a lot of *low* e-values.  Is that good or bad?  Should we
be happy or concerned?

We can take a look at some more stats -- let's look at the `bitscore` column:

```
hist(blast_out$bitscore) 
```

what are we looking for here? (And how would we know?)

(Hint: longer bitscores are better, but even bitscores of ~200 mean a
nucleotide alignment of 200 bp - which is pretty good, no? Here we really
want to rescale the x axis to look at the distribution of bitscores in the
100-300 range.)

Another question - if 'bitscore' is a score of the match, and 'pident'
is the percent identity - is there a relationship between bitscore and
pident?

Well, we can ask this directly with `plot`:

```
plot(blast_out$pident, blast_out$bitscore)
```

why does this plot look the way it does?  (This may take a minute to show
up, note!)

The answer is that bitscores are only *somewhat* related to pident; they
take into account not only the percent identity but the length.  You can
get a napkin sketch estimate of this by doing the following:

```
plot(blast_out$pident  * (blast_out$qend - blast_out$qstart), blast_out$bitscore)
```

which constructs a *new* variable, the percent identity times the
length of the match, and then plots it against bitscore; *this*
correlation looks much better.

### Summary points

This is an example of initial *exploratory data analysis*, in which we poke
around with data to see roughly what it looks like.  This is opposed to
other approaches where we might be trying to do statistical analysis to
confirm a hypothesis.

Typically with small replicate sizes (n < 5) it is hard to do confirmatory
data analysis or hypothesis testing, so a lot of NGS work is done for
*hypothesis generation* and then confirmed via additional experimental
work.

## Some questions for discussion/points to make:

* Why are we using R for this instead of the UNIX command line, or Excel?

  One important thing to note here is that we're looking at a pretty large
  data set - with ease.  It would be much slower to do this in Excel.
  
* What other things could we look at?
