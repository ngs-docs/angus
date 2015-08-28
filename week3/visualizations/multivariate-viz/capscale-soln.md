```R
pcoa<-capscale(decostand(dataset[,-c(1:4)],"total")~dataset$type,distance="bray")
scores(pcoa)$centroids
ggplot(sc)+geom_point(aes(x=MDS1,y=MDS2,colour=info,shape=type))+labs(x="MDS1 (33.0% of variation explained)",y="MDS1 (8.2% of variation explained)")+annotate("text",x=c(-.06,-.23,.25,.02),y=c(.195,-.111,-.10,.521),label=c("sdE3","Unknown 1","Unknown2","wt"))

```
![alt text](https://raw.githubusercontent.com/ryanjw/ngs-3rdweek/master/multivariate-viz/capscale-soln-fig.jpg)