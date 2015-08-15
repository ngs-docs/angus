# reads kmer abundance data output from python script
hista <- read.table("aahist.txt")
histb <- read.table("bbhist.txt")
histc <- read.table("cchist.txt")
histd <- read.table("ddhist.txt")

# opens pdf file for writing:
pdf("kmer10hists.pdf")
# clears graph space
par(mfrow=c(2,2))
# generates 4 histogram graphs on one page
hist(hista[,1],breaks=100,ylim=c(0,20))
hist(histb[,1],breaks=100,ylim=c(0,20))
hist(histc[,1],breaks=100,ylim=c(0,20))
hist(histd[,1],breaks=100,ylim=c(0,40))
# close pdf file
dev.off()
