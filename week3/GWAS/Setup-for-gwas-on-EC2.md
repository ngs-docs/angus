### Launch a new Amazon EC2 instance
instructions: http://angus.readthedocs.org/en/2015/amazon/index.html

### Change prompt
PS1='$ '

### Install programs
```
echo 'deb http://cran.mtu.edu/bin/linux/ubuntu trusty/' | sudo tee -a /etc/apt/sources.list.d/r.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update
sudo apt-get install -y r-base zlib1g-dev vcftools putty-tools plink \
     gcc gcc-c++ libstdc++ gcc-gfortran glibc glibc-devel make blas-devel lapack lapack-devel atlas-devel \
     libatlas-dev libatlas-base-dev
```

and then install plink:
```
wget https://www.cog-genomics.org/static/bin/plink150805/plink_linux_x86_64.zip
unzip plink_linux_x86_64.zip
```

```
### Get data
curl -L https://github.com/ttimbers/SKAT_NGS-2015/blob/master/data.zip?raw=true > data.zip
unzip data.zip
```

### Download and Install packages
```
wget https://cran.r-project.org/src/contrib/Archive/plyr/plyr_1.7.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/stringr/stringr_0.5.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_0.4.1.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/SKAT/SKAT_1.0.7.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/fdrtool/fdrtool_1.2.14.tar.gz
wget https://cran.r-project.org/src/contrib/assertthat_0.1.tar.gz
wget https://cran.r-project.org/src/contrib/R6_2.1.1.tar.gz
wget https://cran.r-project.org/src/contrib/Rcpp_0.12.0.tar.gz
wget https://cran.r-project.org/src/contrib/magrittr_1.5.tar.gz
wget https://cran.r-project.org/src/contrib/lazyeval_0.1.10.tar.gz
wget https://cran.r-project.org/src/contrib/DBI_0.3.1.tar.gz
wget https://cran.r-project.org/src/contrib/BH_1.58.0-1.tar.gz


# Go into R
R

setwd("/home/ubuntu")
install.packages("plyr_1.7.tar.gz", type="source", repos=NULL)	
install.packages("stringr_0.5.tar.gz", type="source", repos=NULL)
install.packages("SKAT_1.0.7.tar.gz", type="source", repos=NULL)
install.packages("fdrtool_1.2.14.tar.gz", type="source", repos=NULL)
install.packages("assertthat_0.1.tar.gz", type="source", repos=NULL)
install.packages("R6_2.1.1.tar.gz", type="source", repos=NULL)
install.packages("Rcpp_0.12.0.tar.gz", type="source", repos=NULL)
install.packages("magrittr_1.5.tar.gz", type="source", repos=NULL)
install.packages("lazyeval_0.1.10.tar.gz", type="source", repos=NULL)
install.packages("DBI_0.3.1.tar.gz", type="source", repos=NULL)
install.packages("BH_1.58.0-1.tar.gz", type="source", repos=NULL)
install.packages("dplyr_0.4.1.tar.gz", type="source", repos=NULL)
```

