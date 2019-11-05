
setwd("/lustre/scr/w/e/weisun/TCGA/bam/")

# ------------------------------------------------------------
# read data
# ------------------------------------------------------------

idx1 = seq(1,601, by=50)
idx2 = c(seq(50,600, by=50), 685)

i = 1
dat = read.table(sprintf("gene_counts_%d_%d.txt",idx1[i], idx2[i]), 
header=TRUE, sep="\t")

for(i in 2:length(idx1)){
  ct = read.table(sprintf("gene_counts_%d_%d.txt",idx1[i], idx2[i]), 
                  header=TRUE, sep="\t")
  dat = cbind(dat, ct)
}

dim(dat)
dat[1:2,1:5]

write.table(dat, file = "../data/gene_counts_all.txt", append = FALSE, 
  quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# ------------------------------------------------------------
# find a cutoff to filter out low expressed genes
# ------------------------------------------------------------

rMin = apply(dat, 1, min)
rMed = apply(dat, 1, median)
r75  = apply(dat, 1, quantile, probs=0.75)
r90  = apply(dat, 1, quantile, probs=0.90)

summary(rMin)
summary(rMed)
summary(r75)
summary(r90)

cor(rMin, rMed)
cor(r75, rMed)
cor(r90, rMed)

pdf("../figures/cts_summary.pdf", width=6, height=6)
par(mfrow=c(2,2), mar=c(5,4,1,1), bty="n")
hist(log10(1+rMin), xlab="log10(min + 1)", main="")
hist(log10(1+rMed), xlab="log10(median + 1)", main="")
hist(log10(1+r75),  xlab="log10(75 percentile + 1)", main="")
hist(log10(1+r90),  xlab="log10(90 percentile + 1)", main="")
dev.off()

summary(rMin[rMed >=10])
summary(rMed[r75 >=20])

table(rMed >=10)
table(r75 >=20)

dim(dat)
dat = dat[which(r75 >=20),]
dim(dat)

# ------------------------------------------------------------
# read in sample list
# ------------------------------------------------------------

aa = read.table("../info/brca_samples2use_after_qc_female_aa.txt", 
  header=TRUE, sep="\t", as.is=TRUE)
dim(aa)
aa$RNA_patientID[1:5]

ea = read.table("../info/brca_samples2use_after_qc_female_caucasian.txt", 
header=TRUE, sep="\t", as.is=TRUE)
dim(ea)
ea$RNA_patientID[1:5]

# ------------------------------------------------------------
# output gene counts for EA or AA population
# ------------------------------------------------------------

waa = match(aa$RNA_patientID, names(dat))
datAA = dat[,waa]
dim(datAA)
datAA[1:2,1:5]

wea = match(ea$RNA_patientID, names(dat))
datEA = dat[,wea]
dim(datEA)
datEA[1:2,1:5]

summary(apply(datAA, 1, median))
summary(apply(datEA, 1, median))

summary(apply(datAA, 1, max))
summary(apply(datEA, 1, max))

datAA[apply(datAA, 1, max)>1e7,]

write.table(datEA, file = "../data/gene_counts_EA.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(datAA, file = "../data/gene_counts_AA.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

