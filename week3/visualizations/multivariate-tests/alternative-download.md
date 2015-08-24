```
download.file("https://raw.githubusercontent.com/ryanjw/ngs-3rdweek/master/datasets/fly_data.txt",destfile="fly_data.txt",method="curl")
dataset<-read.table("fly_data.txt",header=T,sep="\t",check.names=F)
```