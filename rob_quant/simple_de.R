
#from https://github.com/rob-p/angus/blob/2015/files/drosophila_deseq2.R as a template
#Load the DESeq2, tximport and readr libraries
library(DESeq2)
library(tximport)
library(readr)

quant_files <- file.path("quants", list.files("quants"), "quant.sf")
samples <- c("ORE_sdE3_rep1", "ORE_sdE3_rep2", "ORE_wt_rep1","ORE_wt_rep2", "HYB_sdE3_rep1", "HYB_wt_rep1", "SAM_sdE3_rep1","SAM_sdE3_rep2", "SAM_wt_rep1","SAM_wt_rep2", "HYB_sdE3_rep2", "HYB_wt_rep2")
names(quant_files) <- samples
tx2gene <- read.csv("txp_to_gene.csv", col.names=c("TXNAME", "GENEID"))

#Read the input
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
#Get gene-level summary
#txi.sum <- summarizeToGene(txi, tx2gene)

backgrounds <- c(rep("ORE", 4), rep("HYB", 2), rep("SAM", 4), rep("HYB", 2))
genotypes <- c("sdE3", "sdE3", "wt", "wt", "sdE3", "wt", "sdE3", "sdE3", "wt", "wt", "sdE3", "wt")
rna.design <- data.frame(sample=samples, file=quant_files, background=backgrounds, genotype=genotypes)

load.model <- formula(~ genotype)

all.data <- DESeqDataSetFromTximport(txi, rna.design, design=load.model)
all.data <- DESeq(all.data)

plotDispEsts(all.data)

#Now let's look at differentially expressed genes
genotype.results <- as.data.frame(results(all.data, alpha=0.05))

#Viewing the top of the results data frame reveals that there's a lot of genes
# with missing data [no mapped reads]
#That makes sense for this sample dataset, because we kept only reads derived from
# X-linked genes [to speed up run times for the workshop]
print(head(genotype.results))

#Let's pull out the genes that are significant at a false discovery rate of 0.05
sig.genotype.results <- genotype.results[(genotype.results$padj <= 0.05) & !is.na(genotype.results$padj),]

#Sort them
sig.genotype.results <- sig.genotype.results[order(sig.genotype.results$padj, decreasing=F),]

#View them
print(head(sig.genotype.results, n=15))

plot(x=genotype.results$log2FoldChange, y=-log10(genotype.results$padj), pch=16, col=ifelse(genotype.results$padj <= 0.05, "red", "black"))