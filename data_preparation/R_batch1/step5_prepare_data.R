
setwd("/Volumes/Moon/TCGA_BRCA/genotype_normal/")

# ------------------------------------------------------------
# read data
# ------------------------------------------------------------

date()
dat = read.table("birdseed-v2.calls.txt", header=TRUE, na.strings="-1")
date()

dim(dat)
dat[1:2,1:5]

table(dat[,2])
table(is.na(dat[,2]))

# ------------------------------------------------------------
# read in the samples to be used
# ------------------------------------------------------------

sam = read.table("../info/brca_samples2use_after_qc_female.txt", 
  header=TRUE, sep="\t", as.is=TRUE)
dim(sam)
sam[1,]

mat1 = match(sam$DNAnorml_arrayFile, names(dat))

if(any(is.na(mat1))){
  stop("some sampls in the sample file are not in the data\n")
}

dat = dat[,c(1,mat1)]
dim(dat)
dat[1:2,1:5]

if(any(sam$DNAnorml_arrayFile != names(dat)[-1])){
  stop("sample name mismatch\n")
}

# ------------------------------------------------------------
# read in SNP annotation 
# ------------------------------------------------------------

path = "/Users/suninsky/research/data/Affy/Affy6/anno/"
info = read.csv(sprintf("%s/GenomeWideSNP_6.na32.annot.csv", path), 
                as.is=TRUE, comment.char="#")

dim(info)
names(info)
info[1:5,1:4]
table(info$Chromosome)

# ------------------------------------------------------------
# extract data 
# ------------------------------------------------------------

if(! all(dat$probeset_id %in% info$Probe.Set.ID)){
  stop("some probeset_id in data file are not recognized\n")
}

info = info[match(dat$probeset_id, info$Probe.Set.ID),]
dim(info)
all(dat$probeset_id == info$Probe.Set.ID)
table(info$Chromosome)

w2rm = which(info$Chromosome == "---")
length(w2rm)

if(length(w2rm) > 0){
  info = info[-w2rm,]
  dat  = dat[-w2rm,]
}

dim(info)
dim(dat)
table(info$Chromosome)

# ------------------------------------------------------------
# sort by SNP location 
# ------------------------------------------------------------

chrs = info$Chromosome
chrs[which(chrs=="X")]   = "23"
chrs[which(chrs=="Y")]   = "24"
chrs[which(chrs=="MT")]  = "25"

table(chrs)
chrs = as.numeric(chrs)
table(chrs)

pos  = as.numeric(info$Physical.Position)
summary(pos)

od   = order(chrs, pos)

info = info[od,]
dat  = dat[od,]

dim(info)
info[1:5,1:4]

dim(dat)
dat[1:2,1:2]

if(! all(dat$probeset_id == info$Probe.Set.ID)){
  stop("probeset_id mismatch\n")
}

# ------------------------------------------------------------
# remove the SNPs that have more than 5% missing values 
# ------------------------------------------------------------

pdata = data.matrix(dat[,-1])
dim(pdata)

nNA = rowSums(is.na(pdata))
summary(nNA)

w2rm = which(nNA > 0.05*ncol(pdata))
length(w2rm)

if(length(w2rm) > 0){
  dat   = dat[-w2rm,]
  info  = info[-w2rm,]
  pdata = pdata[-w2rm,]
}

dim(dat)
dim(info)
dim(pdata)

# ------------------------------------------------------------
# switch to forward strand
# ------------------------------------------------------------

table(info$Strand)
table(info$Strand.Versus.dbSNP)
table(info$Strand, info$Strand.Versus.dbSNP)

table(info$Allele.A)
table(info$Allele.B)

switchStrand <- function(xx){
	wA = which(xx == "A")
	wC = which(xx == "C")
	wG = which(xx == "G")
	wT = which(xx == "T")
  
  if(length(wA) > 0){ xx[wA] = "T" }
  if(length(wT) > 0){ xx[wT] = "A" }
  if(length(wC) > 0){ xx[wC] = "G" }
  if(length(wG) > 0){ xx[wG] = "C" }
  
  xx
}

wRev = which(info$Strand == "-")

info$Allele.A.Forward = info$Allele.A
info$Allele.A.Forward[wRev] = switchStrand(info$Allele.A[wRev])

info$Allele.B.Forward = info$Allele.B
info$Allele.B.Forward[wRev] = switchStrand(info$Allele.B[wRev])

table(info$Allele.A, info$Allele.A.Forward)
table(info$Allele.B, info$Allele.B.Forward)

# ------------------------------------------------------------
# prune SNPs (50, 5, 0.1)
# using sliding window of 50 SNPs, step size of 5 SNP, and 
# remove one of a pair of SNPs if the R2 is larger than 0.1 
# ------------------------------------------------------------

starts = seq(1, nrow(pdata)-49, 5)
ends   = seq(50, nrow(pdata), 5)

length(starts)
length(ends)

snp2drop = rep(0, nrow(pdata))

r2  = matrix(0, nrow=50, ncol=50)
lr2 = lower.tri(r2, diag = TRUE)

message("pruning SNPs, (50, 5, 0.1)")

for(i in 1:length(starts)){
  
  if(i %% 1000 == 0){
    message(i, " ", date())
  }
  
  p1 = starts[i]
  p2 = ends[i]
  r2 = cor(t(pdata[p1:p2,]), use="pair")
  r2 = r2*r2
  r2[lr2] = 0

  w2 = which(r2 > 0.1, arr.ind=TRUE)
  
  if(nrow(w2) > 0){
    w3 = (p1:p2)[w2]
    snp2drop[w3] = 1
  }
}

# ------------------------------------------------------------
# further check MAF and remove those SNPs with MAF=0
# ------------------------------------------------------------

MAF = 0.5*rowSums(pdata, na.rm=TRUE)/rowSums((!is.na(pdata)), na.rm=TRUE)
MAF = pmin(MAF, 1 - MAF)
summary(MAF)

xx = strsplit(info$Minor.Allele.Frequency, split=" // ", fixed=TRUE)
table(sapply(xx, length))
xx = matrix(unlist(xx), ncol=5, byrow=TRUE)
xx = as.numeric(xx[,1])

cor(MAF, xx)^2

png("../figures/MAF_data_vs_affy_anno.png", width=5, height=5, 
  res=200, units="in")
par(mar=c(5,4,1,1))
smoothScatter(xx, MAF, xlab="MAF @ Caucasian", ylab="MAF @ sample")
dev.off()

table(snp2drop)
table(info$In.Hapmap)
table(MAF > 0)
table(info$Strand != "---")

snp2kp = snp2drop==0
table(snp2kp)
snp2kp = snp2kp & MAF > 0 
table(snp2kp)
snp2kp = snp2kp & info$In.Hapmap=="YES" 
table(snp2kp)
snp2kp = snp2kp & info$Strand != "---"
table(snp2kp)

snp2kp = which(snp2kp)
length(snp2kp)

# ------------------------------------------------------------
# select the SNPs for PCA analysis
# ------------------------------------------------------------

dat4Pr = pdata[snp2kp,]
dim(dat4Pr)
dat4Pr[1:2,1]

inf4Pr = info[snp2kp, c(1:10,19,28,29)]
dim(inf4Pr)
inf4Pr[1:2,]

write.table(inf4Pr, file = "../data/snps_4_PCA.txt", append = FALSE, 
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# generate genotype data in terms of nucleotides for PCA
# ------------------------------------------------------------

dat4PrN = matrix("NN", nrow=nrow(dat4Pr), ncol=ncol(dat4Pr))
geno0   = paste(inf4Pr$Allele.A.Forward, inf4Pr$Allele.A.Forward, sep="")
geno1   = paste(inf4Pr$Allele.A.Forward, inf4Pr$Allele.B.Forward, sep="")
geno2   = paste(inf4Pr$Allele.B.Forward, inf4Pr$Allele.B.Forward, sep="")

for(i in 1:nrow(dat4Pr)){
  w0 = which(dat4Pr[i,] == 0)  
  w1 = which(dat4Pr[i,] == 1)  
  w2 = which(dat4Pr[i,] == 2)
  
  if(length(w0)>0){
    dat4PrN[i,w0] = geno0[i]
  }
  
  if(length(w1)>0){
    dat4PrN[i,w1] = geno1[i]
  }

  if(length(w2)>0){
    dat4PrN[i,w2] = geno2[i]
  }
}

colnames(dat4PrN) = sam$DNAnorml_patientID
dim(dat4PrN)
dat4PrN[1:3,1:4]

table(dat4PrN[,1])

alleles = paste(inf4Pr$Allele.A.Forward, inf4Pr$Allele.B.Forward,sep="/")
db      = data.frame(rs=inf4Pr$dbSNP.RS.ID, alleles=alleles, 
          chrom=inf4Pr$Chromosome, pos=inf4Pr$Physical.Position, 
          dat4PrN, stringsAsFactors=FALSE)
dim(db)
db[1:2,1:6]

write.table(db, file = "../data/genotype_4_PCA_forward.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



