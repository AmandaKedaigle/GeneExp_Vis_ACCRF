#Load HPA per sample data
#Trimmed and aligned and summarizedExperiment-ed as with PDX data
load("after_se_HPA.RData")
se_norm <- se
load("after_se_acc.RData")
se_acc <- se
colData(se_acc)$Disease <- as.factor(rep("ACC",8))
colData(se_norm)$Disease <- as.factor(rep("Cont",6))
se <- cbind(se_acc, se_norm)

#DESeq2
#biocLite("DESeq2")
library("DESeq2")

dds <- DESeqDataSet(se, design = ~Disease)
dds <- dds[ rowSums(counts(dds)) > 1, ]

#Run differential expression (using predefined experimental design)
dds <- DESeq(dds)
(res <- results(dds))
(summary(res))

#plot counts for top gene
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("Disease"))
#MA plot
plotMA(res)

#load annotation information to change gene names
library("AnnotationDbi")
library("org.Hs.eg.db")
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENTREZID",
                     multiVals="first")
resOrdered <- res[order(res$padj),]
resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file="RNAseq_vHPA_results.csv")

#plot gene counts nicer
library(ggplot2)
library(reshape2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
(ebg <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene"))
lens <- sum(width(reduce(ebg)))
lens <- lens/1000

rpk <- assay(se)/lens
tpm <- as.data.frame(round(t(t(rpk)/colSums(rpk))*1000000,1))
tpm$Gene <- mapIds(org.Hs.eg.db,
                        keys=row.names(tpm),
                        column="SYMBOL",
                        keytype="ENTREZID",
                        multiVals="first")
tpm <- tpm[complete.cases(tpm$Gene),]
tpm <- melt(tpm)
tpm = within(tpm, {
  sample = ifelse(variable %in% c("ERR315325.bam","ERR315382.bam","ERR315418.bam","ERR315420.bam","ERR315449.bam","ERR315459.bam"), "Cont", "ACC")
})
gene="NOTCH1"
ggplot(tpm[tpm$Gene==gene,], aes(variable,value)) + geom_point() +
  xlab(NULL) + ylab("Transcripts per Million") + ggtitle(gene) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
save.image("after_deseq.RData")
