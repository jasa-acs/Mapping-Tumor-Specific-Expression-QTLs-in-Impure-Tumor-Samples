
setwd("~/research/TCGA/BRCA/data/")

cts = read.table("gene_counts_all.txt", sep="\t")
dim(cts)
cts[1:2,1:5]

total = colSums(cts)

ct0 = read.table("../info/samples_count.txt", header=TRUE, sep="\t")
dim(ct0)
ct0[1:2,]

all(names(cts) == ct0$sample)

w1 = which(ct0$TReC_after_QC/2 < total/1e6)
ct0[w1,]

summary(ct0$ratio)
ct1 = ct0[ct0$ratio < 0.7,]
ct1[order(ct1$ratio),]

png("../figures/proportion_reads_passed_QC.png", width=5, height=4, units="in", res=300)
par(mar=c(5,4,1,1))
hist(ct0$ratio, breaks=20, main="", xlab="Proportion of peads passed QC")
dev.off()

png("../figures/nreads_genome_vs_exon.png", width=5, height=5, unit="in", res=300)
par(mar=c(5,4,1,1))

plot(ct0$TReC_after_QC, total/1e6, , bty="n", 
  xlab="genome-wide # of reads (million)", cex=0.7, 
  ylab="# of paired-end reads mapped to exons")
abline(0, 0.5,   lwd=2, col="red", lty=1)
abline(0, 0.45,  lwd=2, col="green", lty=2)
abline(0, 0.4,   lwd=2, col="orange", lty=3)

legend("topleft", c("y=0.5x", "y=0.45x", "y=0.4x"), bty="n", 
  lty=c(1,2,3), col=c("red", "green", "orange"), lwd=c(2,2,2))
dev.off()
