library("edgeR")

files <- c(
"ERR458493.fastq.gz.quant.counts",
"ERR458494.fastq.gz.quant.counts",
"ERR458495.fastq.gz.quant.counts",
"ERR458500.fastq.gz.quant.counts",
"ERR458501.fastq.gz.quant.counts",
"ERR458502.fastq.gz.quant.counts"
)

labels=c("A", "B", "E", "C", "D", "F")

data <- readDGE(files)

print(data)

###

group <- c(rep("mut", 3), rep("wt", 3))

dge = DGEList(counts=data, group=group)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# make an MA-plot

et <- exactTest(dge, pair=c("wt", "mut"))
etp <- topTags(et, n=100000)
etp$table$logFC = -etp$table$logFC
pdf("yeast-edgeR-MA-plot.pdf")
plot(
  etp$table$logCPM,
  etp$table$logFC,
  xlim=c(-3, 20), ylim=c(-12, 12), pch=20, cex=.3,
  col = ifelse( etp$table$FDR < .2, "red", "black" ) )
dev.off()

# plot MDS
pdf("yeast-edgeR-MDS.pdf")
plotMDS(dge, labels=labels)
dev.off()

# output CSV for 0-6 hr
write.csv(etp$table, "yeast-edgeR.csv")
