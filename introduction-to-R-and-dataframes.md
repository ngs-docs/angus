# An Introduction to R and Data Analysis

(Guest lecturer for the Tigers - Tracy Teal, Executive Director of Data Carpentry!)

----

For this, we will be working exclusively in RStudio!  Try to connect to a
running RStudio Web server instance -- you can get the Web address by
running this command:

```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

## Install RStudio Web server

If you cannot connect, you'll need to install the prerequisite software:

```
sudo apt-get update && sudo apt-get install -y gdebi-core r-base r-base-dev
```

After that finishes, download and install RStudio itself.
```
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi -n rstudio-server-1.0.143-amd64.deb 
```

And, finally, change the password to something you can remember:
```
sudo passwd tx160085
```

## Install the `tidyverse` packages

As per [these installation instructions](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-16-04-2), we can install the
so-called `tidyverse` packages like so:

```
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
```

and then upgrade `r-base` and install `r-base-dev`:

```
sudo apt-get install -y --allow-unauthenticated r-base r-base-dev \
    libxml2-dev libcurl4-openssl-dev
```

Now we'll want to install `tidyverse` at the command line:

```
cd ~/
cat > install.R <<EOF
install.packages('tidyverse', repos='http://cran.us.r-project.org')
EOF
sudo Rscript install.R
echo 'done!'
```

this will take a long time to run - while it's running you can switch to
the RStudio Web tab in your browser and start working with R!

## Learn!

From this point on, we will be using the standard lesson from
Data Carpentry on R: [R for Data Analysis](https://data-lessons.github.io/R-genomics/).
