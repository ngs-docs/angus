```R
ggplot(sc)+geom_point(aes(x=NMDS1,y=NMDS2,colour=info,shape=type))+theme_bw(base_size=15)+theme(aspect.ratio=1)
```
![alt text](https://raw.githubusercontent.com/ryanjw/ngs-3rdweek/master/multivariate-viz/nmds-soln.jpg)