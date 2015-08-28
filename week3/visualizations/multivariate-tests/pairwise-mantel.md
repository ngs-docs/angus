```R
HYB_subset<-subset(dataset, fly=="HYB")[,-c(1:3)]
ORE_subset<-subset(dataset, fly=="ORE")[,-c(1:3)]
SAM_subset<-subset(dataset, fly=="SAM")[,-c(1:3)]

HYB_dist<-vegdist(decostand(HYB_subset,"pa"),method="jaccard")
ORE_dist<-vegdist(decostand(ORE_subset,"pa"),method="jaccard")
SAM_dist<-vegdist(decostand(SAM_subset,"pa"),method="jaccard")

mantel(ORE_dist,HYB_dist,method="spearman",permutations=9999)
mantel(ORE_dist,SAM_dist,method="spearman",permutations=9999)
mantel(HYB_dist,SAM_dist,method="spearman",permutations=9999)
```