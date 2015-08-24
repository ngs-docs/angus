```R
df<-read.table("fly_data.txt",header=T,check.names=F,sep="\t")
head(df)
df_melt<-melt(df, id=c("fly","type","file"))
library(plyr)
fly_stats<-ddply(df_melt, .(fly,type,variable),summarise,.progress="text",
means=mean(value),
sds=sd(value)
)

head(df_melt)


head(fly_stats)
df_melt<-melt(df, id=c("fly","type","file"))
simmer1<-subset(fly_stats, fly=="HYB" & type=="sdE3")
simmer2<-subset(fly_stats, fly=="SAM" & type=="sdE3")
head(simmer1)
head(subset(simmer2, value > 0))
head(subset(simmer1, means > 0))
simmer1$means<-round(simmer1$means)
for(i in 1:20){
	simmer1<-subset(fly_stats, fly=="HYB" & type=="sdE3")
	simmer1$value<-abs(rnorm(dim(simmer1)[1], simmer1$means,simmer1$sds/4))
	simmer1$fly<-paste("Unknown",i,sep="")
	simmer1$type<-"Unknown 1"
	simmer2<-subset(fly_stats, fly=="SAM" & type=="sdE3")
	simmer2$value<-abs(rnorm(dim(simmer2)[1], simmer2$means,simmer2$sds/4))
	simmer2$fly<-paste("Unknown",i+20,sep="")
	simmer2$type<-"Unknown 2"
	simmer1$file<-"unk"
	simmer2$file<-"unk"
	head(df_melt)
	head(simmer1)
	simmer1<-simmer1[,c(1,2,7,3,6)]
	simmer2<-simmer2[,c(1,2,7,3,6)]
	df_melt<-rbind(df_melt,simmer1,simmer2)
	print(i/20)
}
```
