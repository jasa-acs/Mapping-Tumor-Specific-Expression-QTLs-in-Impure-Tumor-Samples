
setwd("/Volumes/Moon/TCGA_BRCA/data/")

datAA = read.table(file = "gene_counts_AA.txt", sep = "\t", header=TRUE, as.is=TRUE)
dim(datAA)
datAA[1:2,1:5]

datAA = data.matrix(datAA)
tot   = colSums(datAA)
s75   = apply(datAA, 2, quantile, prob=0.75)

cor(tot, s75)

pdf("../figures/AA_total_vs_75_percentile.pdf", width=4, height=4)
par(mar=c(5,4,1,1), bty="n")
plot(tot/1e6, s75/1000, xlab="total reads (million)", 
  ylab="75 percentile (thousand)", cex=0.5)
dev.off()

nDat = t(log10(t((datAA + 1))/s75))

dim(nDat)

pDat = data.frame(id=rownames(nDat), nDat)
dim(pDat)
pDat[1:2,1:5]

write.table(pDat, file = "../data/log_TReC_AA.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# Run PCA only use the AA TCGA samples 
# ------------------------------------------------------------

datR14Pr = nDat - rowMeans(nDat, na.rm=TRUE)

datR14Pr[is.na(datR14Pr)] = 0 
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]
PC3 =  prdatR1$vectors[,3]

setwd("/Volumes/Moon/TCGA_BRCA/figures/")

pdf("PCs_log_TReC_AA.pdf", width=6, height=6)
par(mar=c(5,4,1,1), mfrow=c(2,2))
barplot(prdatR1$values[1:20], main="", 
xlab="Index", ylab="Eigen-value")

par(mar=c(5,4,1,1))
plot(PC1, PC2,  bty="n")
plot(PC1, PC3,  bty="n")
plot(PC2, PC3,  bty="n")

dev.off()

PCs = t(prdatR1$vectors)
PCs = data.frame(id=paste("PC", 1:nrow(PCs), sep=""), PCs)
names(PCs) = names(pDat)

dim(PCs)
PCs[1:2,1:5]

write.table(PCs, file = "../data/PCs_log_TReC_AA.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

