# Exploratory RNAseq data anlysis using RMarkdown

During this lesson, you’ll learn how to use RMarkdown for reproducible data analysis.  We will work with the RNAseq data from the yeast `mut` and `wt` dataset from last week.  

## Getting started

[Start up an m1.medium instance running Ubuntu 16.04 on Jetstream.](jetstream/boot.html)


## Make sure R & RStudio are installed:

```
sudo apt-get update && sudo apt-get -y install gdebi-core r-base
```

Try to connect to a running RStudio Web server instance – you can get the Web address by running this command:
```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

If you cannot connect, download and install RStudio.

```
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi -n rstudio-server-1.0.143-amd64.deb 
```
And, finally, change the password to something you can remember. If your username is different than the one below (i.e. `diblions` or `dibtiger`), you'll need to change that.
```
sudo passwd your_username
```      

## Introduction to RMarkdown 

1. Why RMarkdown?  
	- YAML Header --> renders pretty docs (we will use HTML)   
2. Markdown  
3. Code Chunks  
    - Knitr 
    - Figures 
    - Tables 
4. Inline Code  
5. Breifly  
    - Can do citations & bibliography 
    - Share easily on Rpubs/github  

## Exploratory data analysis with Yeast RNAseq data  

1. Wild Type vs Mutant plots  
	- Facet the plots  
	- Run linear regressions on correlations between the replicates  
2. More Plots  
3. Heatmap
