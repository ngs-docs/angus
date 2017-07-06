# Differential expression analysis with DESeq2

(Mike Love and Rob Patro)

## Upgrade R to the very latest (3.4.x)

```
sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -

sudo apt-get update && sudo apt-get install -y r-base r-base-dev gdebi-core
```

## Make sure you're running RStudio

For this, we will again be working exclusively in RStudio!  Try to connect to a
running RStudio Web server instance -- you can get the Web address by
running this command:

```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

## Install RStudio Web server

If you cannot connect, you'll need to install it:

```
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi -n rstudio-server-1.0.143-amd64.deb
```

And, finally, change the password to something you can remember:
```
sudo passwd tx160085
```

## Install the `DESeq2` prereqs

```
sudo apt-get install -y libxml2 libxml2-dev libcurl4-gnutls-dev libssl-dev
```

and then install DESeq2:

```
curl -O -L https://github.com/ngs-docs/angus/raw/2017/_static/install-deseq2.R

sudo Rscript --no-save install-deseq2.R
```

## Learn!

From this point on, we will be using the lesson from
[Mike Love's asthma repository](https://github.com/mikelove/asthma).  To
get a copy of it, execute:

```
git clone https://github.com/mikelove/asthma
```

and then open `asthma/scripts/asthma.Rmd` in RStudio.
