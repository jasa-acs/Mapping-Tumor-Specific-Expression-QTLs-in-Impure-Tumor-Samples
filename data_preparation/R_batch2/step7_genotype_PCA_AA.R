
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
# read in the samples to be used
# ------------------------------------------------------------

sam = read.table("../info/brca_samples2use_after_qc_female_aa.txt", 
header=TRUE, sep="\t", as.is=TRUE)
dim(sam)
sam[1,]

mat1 = match(sam$DNAnorml_patientID, names(dat))

pdat = as.matrix(dat[,mat1])
dim(pdat)

## numerical version of genotype data
nDat = matrix(NA, nrow=nrow(pdat), ncol=ncol(pdat))
colnames(nDat) = colnames(pdat)

for(i in 1:nrow(nDat)){
	
  if(i %% 10000 == 0){
    cat(i, date(), "\n")	
  }
  
  pdrow = pdat[i,]
  
  alleles = unlist(strsplit(pdrow, split=""))
  alleles = names(sort(table(alleles), decreasing=TRUE))
  
  if(length(alleles)==1){
    a1 = alleles[1]
    a2 = alleles[1]
  }else if (length(alleles)==2){
    a1 = alleles[1]
    a2 = alleles[2]
  }else{
    stop("more than 2 alleles\n") 
  }

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
  
  if(len.g1 + len.g2 + len.g3 + len.na != ncol(pdat)){
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

dim(nDat)
nDat[1:2,1:5]

# ------------------------------------------------------------
# remove the SNPs that have more than 5% missing values 
# ------------------------------------------------------------

nNA = rowSums(is.na(nDat))
summary(nNA)

w2kp = which(nNA < ncol(nDat)*0.05 )
length(w2kp)

dim(nDat)
nDat = nDat[w2kp,]
dim(nDat)

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

pdf("PCs_genotype_AA.pdf", width=6, height=6)
par(mar=c(5,4,1,1), mfrow=c(2,2))
barplot(prdatR1$values[1:20], main="", 
xlab="Index", ylab="Eigen-value")

par(mar=c(5,4,1,1))
plot(PC1, PC2,  bty="n")
plot(PC1, PC3,  bty="n")
plot(PC2, PC3,  bty="n")

dev.off()

setwd("/Volumes/Moon/TCGA_BRCA/data/")

PCs = t(prdatR1$vectors)
PCs = data.frame(id=paste("PC", 1:nrow(PCs), sep=""), PCs)
names(PCs)[-1] = colnames(nDat)

dim(PCs)
PCs[1:2,1:5]

write.table(PCs, file = "PCs_genotype_AA.txt", append = FALSE, 
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


