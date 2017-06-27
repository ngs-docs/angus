# Visualizing BLAST score distributions in RStudio

## Installing and running RStudio on Jetstream

### Install RStudio!

```
sudo apt-get update && sudo apt-get -y install gdebi-core r-base

wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb

sudo gdebi -n rstudio-server-1.0.143-amd64.deb 
```

### Figure out the Web site to connect to.

Because we're using the cloud to run things, everyone will have a different
computer that they're running RStudio on.  To find out what Web site to
connect to to access YOUR server, run:

```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

After running this, copy/paste the URL into your Web browser; you should
see login page. Enter your Jetstream username and password.


If login is unsuccessful, return to the terminal and run:

```
sudo passwd <XSEDE username>
```
You will be prompted to enter a new password:

```
Enter new UNIX password: 
Retype new UNIX password:
```
Return to the browser login page and enter your new password. Note this will not change your XSEDE login info (only affects this instance).


Once R is up and running:

Download the data, and give it a better name
```
download.file("https://github.com/ngs-docs/angus/raw/17a0ba3b1d915de90a5b8bd1fbc1027eba47baf8/_static/shmlast/mouse.1.rna.fna.gz.x.cow.faa.crbl.csv.gz", "shmlast_mouse.rna.fna.gz.x.cow.faa.crbl.csv")
```

Read the object in to R, and name it something that you might remember
```
shmlast_out <- read.csv("shmlast_mouse.rna.fna.gz.x.cow.faa.crbl.csv")
```
Take a look at the data. This is called a dataframe. Use the `head()` command, or the structure `str()` command, or the dimensions command `dim()`.
```
head(shmlast_out)
str(shmlast_out)
dim(surveys)
```
That's a big data frame!
Let's do some data visualization to get a handle on what our blast output looked like
```
hist(shmlast_out$E_scaled)
hist(shmlast_out$bitscore) 
hist(shmlast_out$q_len)
```
