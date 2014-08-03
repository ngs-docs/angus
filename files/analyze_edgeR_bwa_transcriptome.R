transcriptome.data <- read.table("bwa_output_transcriptome_mapped/bwa_transcriptome_counts.txt", header=FALSE, col.names=c("id", "start", "len", "w_1", "w_2", "vg_1", "vg_2"))

library(edgeR)

edger.counts <- transcriptome.data[,c("w_1", "w_2", "vg_1", "vg_2")]
edger.data <- DGEList(counts=edger.counts, group=c(1,1,2,2))

edger.data.norm <- calcNormFactors(edger.data)
edger.data.norm <- estimateCommonDisp(edger.data.norm)
edger.data.norm <- estimateTagwiseDisp(edger.data.norm)
edger.results <- exactTest(edger.data.norm)

topTags(edger.results)

library(qvalue)
p.vals <- edger.results$table$PValue
q.vals <- qvalue(p.vals)$qvalues

transcriptome.data$p <- p.vals
transcriptome.data$q <- q.vals

transcriptome.data[transcriptome.data$q < 0.01,]

fdr.bh <- p.adjust(p.vals, method="BH")