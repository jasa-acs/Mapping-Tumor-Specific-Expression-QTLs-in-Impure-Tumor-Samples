
library(pTReCASE)

# chromsome to work on
chrNum   = 21

#------------------------------------------------------------------#
# read in total read count and gene information
#------------------------------------------------------------------#

datT = read.table(file = "data/gene_counts_EA_filtered.txt", sep = "\t", 
                  header = TRUE, as.is=TRUE)

info = read.table(file = "data/ens_gene_loci.txt", sep = "\t", 
                  header = TRUE, as.is=TRUE)

dim(datT)
dim(info)

datT[1:2,1:5]
info[1:2,]

mati = match(rownames(datT), info$id)
table(is.na(mati))

info = info[mati,]
dim(info)
all(rownames(datT) == info$id)
datT[1:2, 1:5]
info[1:2,]

sams = names(datT)

#------------------------------------------------------------------#
# read in allele-specific read counts
#------------------------------------------------------------------#

datA1 = read.table(file = "data/gene_asCounts_EA_hap1_filtered.txt", 
                   sep = "\t", header = TRUE)
dim(datA1)

datA1[1:2, 1:5]
all(names(datA1) == sams)

datA2 = read.table(file = "data/gene_asCounts_EA_hap2_filtered.txt", 
                   sep = "\t", header = TRUE)
dim(datA2)

datA2[1:2, 1:8]
all(names(datA2) == sams)

all(rownames(datA1) == rownames(datT))
all(rownames(datA2) == rownames(datT))

#------------------------------------------------------------------#
# read in allele-specific read counts
#------------------------------------------------------------------#

genes2kp = which(info$chr==sprintf("chr%s",chrNum))
length(genes2kp)

Y  = t(data.matrix(datT[genes2kp,]))
Y1 = t(data.matrix(datA1[genes2kp,]))
Y2 = t(data.matrix(datA2[genes2kp,]))
info1  = info[genes2kp,]
dim(info1)
info1[1:2,]

#------------------------------------------------------------------#
# read in covariates
#------------------------------------------------------------------#

covs = read.table("data/covariates_EA.txt", sep = "\t", header = TRUE)
dim(covs)
covs[1:2,1:5]

X0 = t(data.matrix(covs[,-1]))
colnames(X0) = covs$id
all(rownames(X0) == sams)

X = X0[,c(1:3,5,6,8:14,16:38)]
X = cbind(X,(X0[,4]+X0[,7]+X0[,15]))
colnames(X)[36] = c("AC+AQ+GI")

dim(X)
X[1:2,1:5]

# read depth
total.reads = colSums(datT)
rd.depth    = apply(datT, 2, quantile, prob=0.75)
cor(total.reads, rd.depth)

X = cbind(log(rd.depth+1), X)

# make sure X is full rank
s1 = svd(X)
summary(s1$d)

#------------------------------------------------------------------#
# read in tumor purity
#------------------------------------------------------------------#

purity.est = read.table(file="data/ABSOLUTE_purity_and_ploidy_BRCA.txt",
                        sep="\t", header=TRUE, as.is=TRUE)
dim(purity.est)
purity.est[1:2,]

purity.est$samp_id = substr(purity.est$individual_id, start = 9,stop = 12)
table(colnames(datT) %in% purity.est$samp_id)

purity.est = purity.est[match(colnames(datT), purity.est$samp_id),]
dim(purity.est)
purity.est[1:2,]

all(colnames(datT)==purity.est$samp_id)

rhos = purity.est$absolute_extract_purity

#------------------------------------------------------------------#
# read in SNP genotype information
# here we use randomly generated genotype data for demonstration
#------------------------------------------------------------------#

gdata = readRDS("data/simulated_gdata.rds")
ginfo = readRDS("data/simulated_ginfo.rds")

dim(gdata)
gdata[1:2,1:5]

dim(ginfo)
ginfo[1:2,]

Z = data.matrix(gdata)

#------------------------------------------------------------------#
# run pTReCASE
#------------------------------------------------------------------#

dim(Y)
dim(Y1)
dim(Y2)
dim(Z)
dim(X)

summary(rhos)
w2kp = which(!is.na(rhos))

Y  = Y[w2kp,]
Y1 = Y1[w2kp,]
Y2 = Y2[w2kp,]
X  = X[w2kp,]
rhos = rhos[w2kp]

# first 300 SNPs
Z  = Z[w2kp,1:300]
ginfo  = ginfo[1:300,]

p1 = pTReCASE_multapp(Y, Y1, Y2, Z, X, rhos, F_out = "pTReCASE_ex_output.txt", 
                      geno_pos = ginfo$V2, gene_start = info1$start,
                      gene_end = info1$end, Perform_CT_co=1)

sessionInfo()

q(save="no")

