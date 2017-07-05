# Visualizing BLAST score distributions in RStudio

If we want to take a look at the results of
[computing orthologs between mouse and cow](running-blast-large-scale.html),
we can't page through that big text file and understand it usefully.
We have to use another program to look at it.  Here, we're going to use
R and RStudio!

...but first we have to install it.

## Installing and running RStudio on Jetstream

### Install RStudio!

The following commands install the prerequisites for RStudio Web,
download the latest version, and then install it. 

```
sudo apt-get update && sudo apt-get -y install gdebi-core r-base
```

After that finishes, download and install RStudio:
```
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi -n rstudio-server-1.0.143-amd64.deb 
```

You should see now see text indicating an RStudio server has started:

```
Jun 27 11:33:40 js-17-66.jetstream-cloud.org systemd[1]: Starting RStudio Server...
Jun 27 11:33:40 js-17-66.jetstream-cloud.org systemd[1]: Started RStudio Server.
```


### Figure out the Web site to connect to.

Because we're using the cloud to run things, everyone will have a different
computer that they're running RStudio on.  To find out what Web site to
connect to to access YOUR server, run:

```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

After running this, copy/paste the URL into your Web browser; you should
see login page. Enter the XSEDE username and password you were given
(should be `tx160085` username, with associated password).

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

Once R is up and running, we'll give you a quick tour of the RStudio
Web interface for those of you who haven't seen it.

----

## Enter some R commands

(Enter the below commands into RStudio, not the command line.)

Download the precomputed data, and give it a better name; you could also
use your precomputed shmlast results --

```
download.file("https://github.com/ngs-docs/angus/raw/17a0ba3b1d915de90a5b8bd1fbc1027eba47baf8/_static/shmlast/mouse.1.rna.fna.gz.x.cow.faa.crbl.csv.gz", "shmlast_mouse.rna.fna.gz.x.cow.faa.crbl.csv")
```

Next, read the object in to R, and name it something that you might remember

```
shmlast_out <- read.csv("shmlast_mouse.rna.fna.gz.x.cow.faa.crbl.csv")
```

Now we can take a look at the data in a slightly nicer way!

This is called a dataframe, which is a sort of R-ish version of a
spreadsheet with named columns. Use the `head()` command to take a
look at the beginning of it all:

```
head(shmlast_out)
```

this is a generic R command that acts very much like the UNIX command
line `head` command but gives you a slightly nicer more structured view.
In RStudio you can also use the `view` command,

```
View(shmlast_out)
```
which is much nicer altogether!

Another useful command is `dim` which will tell you the DIMENSIONS of this
data frame:

```
dim(shmlast_out)
```

That's a big data frame! 132,900 rows (and 17 columns!)

Let's do some data visualization to get a handle on what our blast output looked like: first, let's look at the `E_scaled` column.

```
hist(shmlast_out$E_scaled)
```

This is telling us that MOST of the values in the `E_scaled` column are
quite high.  What does this mean? How do we figure out what this is?

(Hint: the [shmlast documentation](https://github.com/camillescott/shmlast) should tell us! Go to that page and search for "To fit the model".)

So these are a lot of *low* e-values.  Is that good or bad?  Should we
be happy or concerned that 

We can take a look at some more stats -- let's look at the `bitscore` column:

```
hist(shmlast_out$bitscore) 
```

what are we looking for here? (And how would we know?)

(Hint: longer bitscores are better, but even bitscores of ~1000 mean a
nucleotide alignment of 1000 bp - which is pretty good, no? Here we really
want to rescale the x axis to look at the distribution of bitscores in the
100-300 range.)

We can also look at the length of the queries, which are the mouse sequences
in this case. 

```
hist(shmlast_out$q_len)
```
 Compare this to the bitscores... do things match? Are
most mouse sequences getting matched by something of equivalent length?

Well, we can ask this directly with `plot`:

```
plot(shmlast_out$q_len, shmlast_out$bitscore)
```

why does this plot look the way it does?  (This may take a minute to show
up, note!)

(The bitscores are limited by the length of the sequences! You can't get a
*longer* bitscore than you have bases to align.)

### Summary points

This is an example of initial *exploratory data analysis*, in which we poke
around with data to see roughly what it looks like.  This is opposed to
other approaches where we might be trying to do statistical analysis to
confirm a hypothesis.

Typically with small sample sizes (n < 5) it is hard to do confirmatory
data analysis or hypothesis testing, so a lot of NGS work is done for
*hypothesis generation* and then confirmed via additional experimental
work.

## Some questions for discussion/points to make:

* Why are we using R for this instead of the UNIX command line, or Excel?

  (We'll talk more about this later, too!)
  
  One important thing to note here is that we're looking at a pretty large
  data set - with ease.  It would be much slower to do this in Excel.
  
* What other things could we look at?

* Have we done a basic check of just *looking* at the data? Go back and
  look at the data frame!  Do the gene name assignments look right?
  
  How might we do this a bit more systematically, while still "looking"
  at things? Try googling 'choose random rows from data frame' and then
  run
  
  ```
  shmlast_sub = shmlast_out[sample(nrow(shmlast_out), 10),]
  View(shmlast_sub)
  ```
  
* What does the following code do?

  ```
  tmp <- subset(shmlast_out, q_len >= 8000 & q_len <= 11000 & bitscore <=2000)    
  functions <- tmp[, c("q_name", "s_name")]
  ```

  and do you notice anything interesting about the names?  (They're
  all predicted/inferred genes.) What does this suggest about that
  "quadrant" of the data?
