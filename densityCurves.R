args=(commandArgs(TRUE))
label=args[1]
table1=paste("known",label,sep=".")
table2=paste("novel",label,sep=".")


library(ggplot2)

dataA=read.table(table1)
dataB=read.table(table2)

dataA$V1="known"
dataB$V1="novel"

allData <- rbind(dataA, dataB)

outputPDF <- paste(label,"pdf",sep=".")
pdf(outputPDF)
ggplot(allData, aes(allData$V2, fill = V1)) + geom_density(alpha = 0.2)
dev.off()
