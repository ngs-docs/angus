# Visualizing BLAST score distributions in RStudio

## Installing and running RStudio on Jetstream

### Install RStudio!

```
sudo apt-get update && sudo apt-get -y install gdebi-core r-base

wget https://download2.rstudio.org/rstudio-server-1.0.136-amd64.deb

sudo gdebi -n rstudio-server-1.0.136-amd64.deb
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
