## Multivariate Visualizations with NGS data

We previously explored multivarate tests with the fly data [(lesson here)](https://github.com/ryanjw/ngs-3rdweek/blob/master/multivariate-tests/tests.md), and now we are going to work on different ways of visualizing the data. 

Here we are going to use some multivariate visualizations with a particular goal in mind.  We assume that in addition to the fly data we had previously, 40 more fly transcriptomes were sequenced.  Unfortunately, someone forgot to label the samples.  We know that two different flies (combination of `fly` and `type`) were sequenced.  

Let's download the new data
```R
URL<-("https://raw.githubusercontent.com/ryanjw/ngs-3rdweek/master/multivariate-viz/fly_data_with_unknowns.txt")
dataset<-read.table(textConnection(getURL(URL)),header=T,check.names=F,sep="\t")
```
If you can't use `getURL` modify the URL in the code [here](https://github.com/ryanjw/ngs-3rdweek/blob/master/multivariate-tests/alternative-download.md)

Let's look at the first and last few rows of the data

```R
head(dataset[,1:10])
tail(dataset[,1:10])
```

We can see that we have the addition of 40 unknown flies.  [In reality, these data were simulated from the previous samples, and the rough code can be seen here.](https://github.com/ryanjw/ngs-3rdweek/blob/master/multivariate-viz/data-sim.md)

The first visualization we will perform is based on nonmetric multidimensional scaling.  In short, this algorithmic approach samples are placed in N-dimensional space based on a distance matrix.  We then view a 2D representation of these differences and are given a statistic called *stress* which, when minimized, the visualization is easiest to interpret (Smaller Stress, Better Visualization)

Let's generate an NMDS using the `metaMDS` function from the `vegan` package.  We will make an object called `sc` that has the scores (i.e. coordinates in space) for each sample along with metadata from our samples
```R
 mds<-metaMDS(dataset[,-c(1:4)],distance="bray",autotransform=F,k=3 )
 sc<-data.frame(scores(mds),dataset[,1:3])
 head(sc)
 ```

 Let's visualize these relationships via a scatterplot in ggplot
 ```R
 ggplot(sc)+geom_point(aes(x=NMDS1,y=NMDS2))
 ```
 ![alt text](https://raw.githubusercontent.com/ryanjw/ngs-3rdweek/master/multivariate-viz/plain-nmds.jpg)

 While we can see the points in space, colors and shapes would be nice to help us understand the indentity of each point.

 ##Challenge
 Add colors and shapes based on the metadata within our dataset.  Do this by passing a variable name like `fly` or `type` to the arguments `shape=...` or `colour=...` within the `aes()` function. [*See solution here*](https://github.com/ryanjw/ngs-3rdweek/blob/master/multivariate-viz/nmds-soln.md)

 # PCoA and/or PCA

 Another popular method of visualization is Principle Components Analysis (PCA), and often Principle Components Analysis is used.  What's the difference?  PCA and PCoA are performed differently, but PCoA produces identical results as PCA when the euclidian distance is used.  Therefore we will focus on PCoA as it can work better for non-normally distributed data. 

 Let's produce a PCA first
 ```R
pca<-capscale(dataset[,-c(1:4)]~1)
pca
```
```R
Call: capscale(formula = dataset[, -c(1:4)] ~ 1)

               Inertia Rank
Total         9.59e+08     
Unconstrained 9.59e+08   51
Inertia is mean squared Euclidean distance 

Eigenvalues for unconstrained axes:
     MDS1      MDS2      MDS3      MDS4      MDS5      MDS6      MDS7      MDS8 
728233494  55344575  30956887  26646256  20802970  11395822   9047108   8339889 
(Showed only 8 of all 51 unconstrained eigenvalues)
```
```R
eigs<-eigenvals(pca)
eigs/sum(eigs)
```
```R
     MDS1      MDS2      MDS3      MDS4      MDS5      MDS6      MDS7      MDS8      MDS9     MDS10     MDS11     MDS12     MDS13     MDS14     MDS15     MDS16 
0.7593702 0.0577109 0.0322805 0.0277856 0.0216924 0.0118831 0.0094339 0.0086965 0.0073555 0.0068737 0.0055017 0.0043177 0.0038954 0.0038647 0.0029479 0.0029354 
    MDS17     MDS18     MDS19     MDS20     MDS21     MDS22     MDS23     MDS24     MDS25     MDS26     MDS27     MDS28     MDS29     MDS30     MDS31     MDS32 
0.0025209 0.0023627 0.0021998 0.0019720 0.0018604 0.0017608 0.0017294 0.0014893 0.0013986 0.0012682 0.0011862 0.0010923 0.0010440 0.0009666 0.0008388 0.0007844 
    MDS33     MDS34     MDS35     MDS36     MDS37     MDS38     MDS39     MDS40     MDS41     MDS42     MDS43     MDS44     MDS45     MDS46     MDS47     MDS48 
0.0007641 0.0007490 0.0007321 0.0006646 0.0006423 0.0006083 0.0005679 0.0005434 0.0004940 0.0004746 0.0004725 0.0004220 0.0003919 0.0003593 0.0003352 0.0003105 
    MDS49     MDS50     MDS51 
0.0002960 0.0000866 0.0000660 
```
```R
sc<-data.frame(scores(pca)$sites,dataset[,1:4])
ggplot(sc)+geom_point(aes(x=MDS1,y=MDS2,colour=info,shape=type))+labs(x="MDS1 (75.9% of variation explained)",y="MDS1 (5.8% of variation explained)")
```
![alt text](https://github.com/ryanjw/ngs-3rdweek/raw/master/multivariate-viz/pca.jpg) 

Mirroring the tests we performed earlier, let's do the same thing but with analysis based on composition
```R  
pcoa<-capscale(decostand(dataset[,-c(1:4)],"total")~1,distance="bray")
pcoa
eigs<-eigenvals(pcoa)
eigs/sum(eigs)
sc<-data.frame(scores(pcoa)$sites,dataset[,1:4])
ggplot(sc)+geom_point(aes(x=MDS1,y=MDS2,colour=info,shape=type))+labs(x="MDS1 (33.0% of variation explained)",y="MDS1 (8.2% of variation explained)")
```
![alt text](https://raw.githubusercontent.com/ryanjw/ngs-3rdweek/master/multivariate-viz/pcoa.jpg)

We now see some nice separation that can help us interpret what the samples originate from, any ideas?

We can also constrain the variation within the plot to get a better idea of how things are falling out.  Now we assign variables within the capscale function to do this.  
```R
pcoa<-capscale(decostand(dataset[,-c(1:4)],"total")~dataset$fly,distance="bray")
scores(pcoa)$centroids
ggplot(sc)+geom_point(aes(x=MDS1,y=MDS2,colour=info,shape=type))+labs(x="MDS1 (33.0% of variation explained)",y="MDS1 (8.2% of variation explained)")+annotate("text",x=c(-.157,-.077,.17),y=c(-.004,.670,.156),label=c("HYB","ORE","SAM"))
```
![alt text](https://raw.githubusercontent.com/ryanjw/ngs-3rdweek/master/multivariate-viz/contrained-fly-pcoa.jpg))

Can we guess where our unknown flies originated from?

##Challenge

Determine the fly `type` using the same methodology as above.  [*See solution here*](https://github.com/ryanjw/ngs-3rdweek/blob/master/multivariate-viz/capscale-soln.md)


