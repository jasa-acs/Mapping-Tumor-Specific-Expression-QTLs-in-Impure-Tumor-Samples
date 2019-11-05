
setwd("~/research/TCGA/BRCA/data")

# ------------------------------------------------------------
# read data
# ------------------------------------------------------------

dat = read.table("gene_counts_EA.txt", header=TRUE, sep="\t")
dim(dat)
dat[1:2,1:5]

datAs1 = read.table("gene_asCounts_EA_hap1.txt", header=TRUE, sep="\t")
dim(datAs1)
datAs1[1:2,1:5]

datAs2 = read.table("gene_asCounts_EA_hap2.txt", header=TRUE, sep="\t")
dim(datAs2)
datAs2[1:2,1:5]

all(rownames(dat) %in% rownames(datAs1))
all(rownames(dat) %in% rownames(datAs2))

mat1 = match(rownames(dat), rownames(datAs1))
mat2 = match(rownames(dat), rownames(datAs2))
all(mat1 == mat2)

# ------------------------------------------------------------
# find a cutoff to filter out low expressed genes
# ------------------------------------------------------------

q1s  = apply(datAs1, 1, quantile, probs=c(0.75,0.90))
q2s  = apply(datAs2, 1, quantile, probs=c(0.75,0.90))

qs = t(rbind(q1s, q2s))
colnames(qs) = c("hap1_75", "hap_90", "hap2_75", "hap2_90")
round(cor(qs),2)

pdf("../figures/asCts_summary.pdf", width=6, height=6)
par(mfrow=c(2,2), mar=c(5,4,1,1), bty="n")
hist(log10(1+qs[,1]),  xlab="log10(hap1 75 percentile + 1)", main="")
hist(log10(1+qs[,2]),  xlab="log10(hap1 90 percentile + 1)", main="")

hist(log10(1+qs[,3]),  xlab="log10(hap2 75 percentile + 1)", main="")
hist(log10(1+qs[,4]),  xlab="log10(hap2 90 percentile + 1)", main="")
dev.off()

# ------------------------------------------------------------
# match gene IDs
# ------------------------------------------------------------

datAs1 = datAs1[mat1,]
datAs2 = datAs2[mat2,]

dim(datAs1)
dim(datAs2)

mat3 = match(names(dat), names(datAs1))
mat4 = match(names(dat), names(datAs2))
all(mat3==mat4)
all(names(datAs1)[mat3] == names(dat))

datAs1 = datAs1[,mat3]
datAs2 = datAs2[,mat4]

# ------------------------------------------------------------
# plot it
# ------------------------------------------------------------

TReC = colSums(dat)/1e6
ASE  = (colSums(datAs1) + colSums(datAs2))/1e6

l1 = lm(ASE ~ -1 + TReC)
w1 = which(ASE/TReC < 0.025)
w1

l2 = lm(ASE[-w1] ~ -1 + TReC[-w1])

pdf("../figures/total_counts_vs_as_counts.pdf", width=6, height=3)
par(mfrow=c(1,2), mar=c(5,4,1,1), bty="n")

plot(TReC, ASE, xlab="Total # of reads per sample", 
ylab="Total # of AS reads", cex=0.6)
mtext(side=1, line=4, "(million)")
abline(l1, lwd=2, col="orange")
legend("topleft", legend=sprintf("y = %.3f x", l1$coef), bty="n")
points(TReC[w1], ASE[w1], pch=19, col="red", cex=0.8)

plot(TReC[-w1], ASE[-w1], xlab="Total # of reads per sample", 
ylab="Total # of AS reads", cex=0.6)
mtext(side=1, line=4, "(million)")
abline(l2, lwd=2, col="orange")
legend("topleft", legend=sprintf("y = %.3f x", l2$coef), bty="n")

dev.off()

# ------------------------------------------------------------
# output gene counts for EA or AA population
# ------------------------------------------------------------

write.table(dat[,-w1], file = "gene_counts_EA_filtered.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(datAs1[,-w1], file = "gene_asCounts_EA_hap1_filtered.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(datAs2[,-w1], file = "gene_asCounts_EA_hap2_filtered.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
