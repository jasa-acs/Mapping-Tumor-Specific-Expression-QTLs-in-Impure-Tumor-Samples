
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

sam = read.table("../info/brca_samples2use_after_qc_female_caucasian.txt", 
header=TRUE, sep="\t", as.is=TRUE)
dim(sam)
sam[1,]

w1 = which(sam$DNAnorml_patientID == "A15R")
w1
sam = sam[-w1,]
dim(sam)

mat1 = match(sam$DNAnorml_arrayFile, names(dat))

if(any(is.na(mat1))){
  stop("some sampls in the sample file are not in the data\n")
}

datEA = dat[,c(1,mat1)]
dim(datEA)
datEA[1:2,1:5]

if(any(sam$DNAnorml_arrayFile != names(datEA)[-1])){
  stop("sample name mismatch\n")
}
names(datEA)[-1] = sam$DNAnorml_patientID
names(datEA)[1]  = "id"


dim(datEA)
datEA[1:2,1:5]

# ------------------------------------------------------------
# remove the SNPs that have more than 5% missing values 
# ------------------------------------------------------------

pdata = data.matrix(datEA[,-1])
dim(pdata)

nNA = rowSums(is.na(pdata))
summary(nNA)

w2rm = which(nNA > 0.05*ncol(pdata))
length(w2rm)

if(length(w2rm) > 0){
  datEA = datEA[-w2rm,]
  pdata = pdata[-w2rm,]
}

dim(datEA)
dim(pdata)

# ------------------------------------------------------------
# remove the SNPs with MAF smaller than 1% 
# ------------------------------------------------------------

MAF = 0.5*rowSums(pdata, na.rm=TRUE)/rowSums((!is.na(pdata)), na.rm=TRUE)
MAF = pmin(MAF, 1 - MAF)
summary(MAF)

w2rm = which(MAF <0.01)
length(w2rm)

if(length(w2rm) > 0){
  datEA = datEA[-w2rm,]
  pdata = pdata[-w2rm,]
}

dim(datEA)
dim(pdata)
datEA[1:2,1:8]

write.table(datEA, file = "datEA.txt", append = FALSE, quote = FALSE, 
  sep = "\t", row.names = FALSE, col.names = TRUE)


# ------------------------------------------------------------
# read in the samples to be used
# ------------------------------------------------------------

sam = read.table("../info/brca_samples2use_after_qc_female_aa.txt", 
header=TRUE, sep="\t", as.is=TRUE)
dim(sam)
sam[1,]

mat1 = match(sam$DNAnorml_arrayFile, names(dat))

if(any(is.na(mat1))){
  stop("some sampls in the sample file are not in the data\n")
}

datAA = dat[,c(1,mat1)]
dim(datAA)
datAA[1:2,1:5]

if(any(sam$DNAnorml_arrayFile != names(datAA)[-1])){
  stop("sample name mismatch\n")
}
names(datAA)[-1] = sam$DNAnorml_patientID
names(datAA)[1]  = "id"


dim(datAA)
datAA[1:2,1:5]

# ------------------------------------------------------------
# remove the SNPs that have more than 5% missing values 
# ------------------------------------------------------------

pdata = data.matrix(datAA[,-1])
dim(pdata)

nNA = rowSums(is.na(pdata))
summary(nNA)

w2rm = which(nNA > 0.05*ncol(pdata))
length(w2rm)

if(length(w2rm) > 0){
  datAA = datAA[-w2rm,]
  pdata = pdata[-w2rm,]
}

dim(datAA)
dim(pdata)

# ------------------------------------------------------------
# remove the SNPs with MAF smaller than 5% 
# ------------------------------------------------------------

MAF = 0.5*rowSums(pdata, na.rm=TRUE)/rowSums((!is.na(pdata)), na.rm=TRUE)
MAF = pmin(MAF, 1 - MAF)
summary(MAF)

w2rm = which(MAF <0.05)
length(w2rm)

if(length(w2rm) > 0){
  datAA = datAA[-w2rm,]
  pdata = pdata[-w2rm,]
}

dim(datAA)
dim(pdata)
datAA[1:2,1:8]

write.table(datAA, file = "datAA.txt", append = FALSE, quote = FALSE, 
sep = "\t", row.names = FALSE, col.names = TRUE)



