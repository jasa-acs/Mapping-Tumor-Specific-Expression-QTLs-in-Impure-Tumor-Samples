
# ------------------------------------------------------------
# read data
# ------------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/data/")

dat = read.table("genotype_4_PCA_forward.txt", header=TRUE, 
	na.strings="NN", as.is=TRUE)

dim(dat)
dat[1:2,1:5]

table(dat$chrom)

# ------------------------------------------------------------
# read CEU data
# ------------------------------------------------------------

chrs = c(1:22, "X", "Y")

setwd("~/research/data/HapMap/genotypes_CEU/")

for(chr1 in chrs){
  message(chr1, " ", date())
  f1   = sprintf("genotypes_chr%s_CEU_r27_nr.b36_fwd.txt", chr1)
  
  ceuD = read.table(f1, as.is=TRUE, comment.char="", 
									header=TRUE, na.strings="NN")
  
  if(chr1 == "1"){
    datCEU = matrix("", nrow=nrow(dat), ncol=ncol(ceuD))
    colnames(datCEU) = names(ceuD)
  }
  
  snp2use = intersect(ceuD$rs, dat$rs)
  mat1    = match(snp2use, dat$rs)
  mat2    = match(snp2use, ceuD$rs)
  ceuD    = as.matrix(ceuD)

  datCEU[mat1,] = ceuD[mat2,]
}

dim(datCEU)
datCEU[1:2,1:13]
table(datCEU[,3])

nNA = colSums(is.na(datCEU))
w2kp = which(nNA/nrow(datCEU) < 0.1)
length(w2kp)

datCEU = datCEU[,w2kp]

# ------------------------------------------------------------
# read YRI data
# ------------------------------------------------------------

setwd("~/research/data/HapMap/genotypes_YRI/")

for(chr1 in chrs){
  message(chr1, " ", date())
  f1   = sprintf("genotypes_chr%s_YRI_r27_nr.b36_fwd.txt", chr1)
  
  yriD = read.table(f1, as.is=TRUE, comment.char="", 
                    header=TRUE, na.strings="NN")
  
  if(chr1 == "1"){
    datYRI = matrix("", nrow=nrow(dat), ncol=ncol(yriD))
    colnames(datYRI) = names(yriD)
  }
  
  snp2use = intersect(yriD$rs, dat$rs)
  mat1    = match(snp2use, dat$rs)
  mat2    = match(snp2use, yriD$rs)
  yriD    = as.matrix(yriD)
  
  datYRI[mat1,] = yriD[mat2,]
}

dim(datYRI)
datYRI[1:2,1:13]
table(datYRI[,3])


nNA = colSums(is.na(datYRI))
w2kp = which(nNA/nrow(datYRI) < 0.1)
length(w2kp)

table(datYRI[,12],useNA="ifany")
datYRI = datYRI[,w2kp]

# ------------------------------------------------------------
# read CHB data
# ------------------------------------------------------------

setwd("~/research/data/HapMap/genotypes_CHB/")

for(chr1 in chrs){
  message(chr1, " ", date())
  f1   = sprintf("genotypes_chr%s_CHB_r27_nr.b36_fwd.txt", chr1)
  
  chbD = read.table(f1, as.is=TRUE, comment.char="", 
                    header=TRUE, na.strings="NN")
  
  if(chr1 == "1"){
    datCHB = matrix("", nrow=nrow(dat), ncol=ncol(chbD))
    colnames(datCHB) = names(chbD)
  }
  
  snp2use = intersect(chbD$rs, dat$rs)
  mat1    = match(snp2use, dat$rs)
  mat2    = match(snp2use, chbD$rs)
  chbD    = as.matrix(chbD)
  
  datCHB[mat1,] = chbD[mat2,]
}

dim(datCHB)
datCHB[1:2,1:13]
table(dat$chrom[datCHB[,3]==""])

nNA = colSums(is.na(datCHB))
w2kp = which(nNA/nrow(datCHB) < 0.1)
length(w2kp)

table(datCHB[,16],useNA="ifany")
datCHB = datCHB[,w2kp]

# ------------------------------------------------------------
# find intersection
# ------------------------------------------------------------

dim(dat)
snp2use = intersect(dat$rs, datCEU[,1])
length(snp2use)

snp2use = intersect(snp2use, datYRI[,1])
length(snp2use)

snp2use = intersect(snp2use, datCHB[,1])
length(snp2use)

which(snp2use == "")

# ------------------------------------------------------------
# merge them into one dataset
# ------------------------------------------------------------

table(datCEU[,5])
table(datYRI[,5])
table(datCHB[,5])

dTCGA = dat[match(snp2use, dat$rs),]
dCEU  = datCEU[match(snp2use, datCEU[,1]),]
dYRI  = datYRI[match(snp2use, datYRI[,1]),]
dCHB  = datCHB[match(snp2use, datCHB[,1]),]

dim(dTCGA)
dim(dCEU)
dim(dYRI)
dim(dCHB)

all(dCEU[,1] == dTCGA$rs)
all(dCEU[,1] == dYRI[,1])
all(dCEU[,1] == dCHB[,1])

table(dTCGA$alleles, dCEU[,2])
table(dCEU[,2], dYRI[,2])
table(dCEU[,2], dCHB[,2])

mD = cbind(dTCGA[,-(1:4)], dCEU[,-(1:11)], dYRI[,-(1:11)], dCHB[,-(1:11)])
mD = as.matrix(mD)
dim(mD)

lab = c(rep("TCGA",ncol(dTCGA)-4), rep("CEU", ncol(dCEU)-11))
lab = c(lab, rep("YRI",ncol(dYRI)-11), rep("CHB", ncol(dCHB)-11))
table(lab)

nNA = colSums(is.na(mD))
summary(nNA)

# ------------------------------------------------------------
# row by row, fix it
# ------------------------------------------------------------

alleles = as.character(dTCGA$alleles)  
alleles = strsplit(alleles, "/")
n.allel = sapply(alleles, length)
table(n.allel)

## numerical version of genotype data
nDat = matrix(NA, nrow=nrow(mD), ncol=ncol(mD))
colnames(nDat) = colnames(mD)

for(i in 1:nrow(mD)){
	
  if(i %% 10000 == 0){
    cat(i, date(), "\n")	
  }

  pdrow = mD[i,]

  a1 = alleles[[i]][1]
  a2 = alleles[[i]][2]

  g1 = paste(a1, a1, sep="")
  g2 = paste(a1, a2, sep="")
  g3 = paste(a2, a2, sep="")
  g4 = paste(a2, a1, sep="")
        
  which.g1 = which(pdrow==g1)
  which.g2 = which(pdrow==g2 | pdrow==g4)
  which.g3 = which(pdrow==g3)
  which.na = which(is.na(pdrow))

  len.g1 = length(which.g1)
  len.g2 = length(which.g2)
  len.g3 = length(which.g3)
  len.na = length(which.na)	

  if(len.g1 + len.g2 + len.g3 + len.na != ncol(mD)){
    next
  }

  nDat[i, which.g2] = 1

  if(len.g1 >= len.g3){
    nDat[i, which.g1] = 0
    nDat[i, which.g3] = 2
  }else{
    nDat[i, which.g1] = 2
    nDat[i, which.g3] = 0
  }

}

nNA = rowSums(is.na(nDat))
summary(nNA)

w2kp = which(nNA < ncol(nDat)*0.05 )
length(w2kp)

dim(nDat)
nDat = nDat[w2kp,]
dim(nDat)

# ------------------------------------------------------------
# PCA together with HapMap samples
# ------------------------------------------------------------

dim(nDat)
nDat[1:2,1:5]

table(lab)
datR14Pr = nDat - rowMeans(nDat, na.rm=TRUE)

datR14Pr[is.na(datR14Pr)] = 0 
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:10]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]

setwd("/Volumes/Moon/TCGA_BRCA/figures/")

pdf("PCs_withHapMap.pdf", width=6, height=3)
par(mar=c(5,4,1,1), mfrow=c(1,2))
barplot(prdatR1$values[1:10], main="", 
  xlab="Index", ylab="Eigen-value")

w0 = which(lab == "TCGA")
w1 = which(lab == "CEU")
w2 = which(lab == "YRI")
w3 = which(lab == "CHB")

par(mar=c(5,4,1,1))
plot(PC1, PC2, type="n", bty="n")
points(PC1[w1], PC2[w1], col="darkred",   pch=19, cex=0.6)
points(PC1[w2], PC2[w2], col="darkgreen", pch=0,  cex=0.6)
points(PC1[w3], PC2[w3], col="darkblue",  pch=1,  cex=0.6)
points(PC1[w0], PC2[w0], col="black",     pch=3,  cex=0.4)
legend("topright", c("TCGA", "CEU", "YRI", "CHB"), pch=c(3,19,0,1), bty="n", 
	col=c("black", "darkred", "darkgreen", "darkblue"))

abline(h=0.00, col="grey")
abline(v=0.00, col="grey")

dev.off()

table(lab[PC1 > 0.08])
table(lab[PC2 > 0.08])

table(PC1 < 0.01, PC2 < 0.01)
table(lab[PC1 < 0.01 & PC2 < 0.01])

table(PC1 < 0.0, PC2 < 0.0)
table(lab[PC1 < 0.0 & PC2 < 0.0])

# ------------------------------------------------------------
# write out PCA results
# ------------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/data/")

colnames(prdatR1$vectors) = paste("PC", 1:ncol(prdatR1$vectors), sep="")

pcDat = cbind(colnames(nDat), lab, prdatR1$vectors[,1:10])
colnames(pcDat)[1] = "sample"
dim(pcDat)
pcDat[1:2,1:5]

write.table(pcDat, file = "PCs_with_hapmap_samples.txt", 
  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
  col.names = TRUE)

# ------------------------------------------------------------
# Run PCA only use the Caucasian TCGA samples 
# ------------------------------------------------------------

w2use   = (lab=="TCGA" & PC1 < 0.0 & PC2 < 0.0)
sam2use = colnames(nDat)[w2use]
length(sam2use)

datR14Pr = nDat[,w2use] - rowMeans(nDat[,w2use], na.rm=TRUE)

datR14Pr[is.na(datR14Pr)] = 0 
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]
PC3 =  prdatR1$vectors[,3]

setwd("/Volumes/Moon/TCGA_BRCA/figures/")

pdf("PCs_Caucasian_without_HapMap.pdf", width=6, height=6)
par(mar=c(5,4,1,1), mfrow=c(2,2))
barplot(prdatR1$values[1:20], main="", 
xlab="Index", ylab="Eigen-value")

par(mar=c(5,4,1,1))
plot(PC1, PC2,  bty="n")
plot(PC1, PC3,  bty="n")
plot(PC2, PC3,  bty="n")

dev.off()

# ------------------------------------------------------------
# read in covaraites 
# ------------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/info/")

sam = read.table("brca_samples2use_after_qc_female.txt", 
header=TRUE, sep="\t", as.is=TRUE)
dim(sam)
sam[1,]

sam = sam[match(sam2use, sam$DNAnorml_patientID),]
all(sam2use == sam$DNAnorml_patientID)


pvals = matrix(NA, nrow=20, ncol=4)

for(i in 1:20){
  PCi = prdatR1$vectors[,i]
  ai  = anova(lm(PCi ~  as.factor(sam$DNAnorml_plate) 
                      + as.factor(sam$DNAnorml_institution)
                      + as.factor(sam$DNAnorml_type)
                      + as.factor(sam$DNAnorml_portion)))
  pvals[i,] = ai$Pr[1:4]

}

cbind(signif(prdatR1$values[1:20]), signif(pvals,3))

table(sam$DNAnorml_plate)
table(sam$DNAnorml_institution)
table(sam$DNAnorml_plate, sam$DNAnorml_institution)

# ------------------------------------------------------------
# write the new sample list 
# ------------------------------------------------------------

write.table(sam, file = "brca_samples2use_after_qc_female_caucasian.txt", 
append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
col.names = TRUE)

