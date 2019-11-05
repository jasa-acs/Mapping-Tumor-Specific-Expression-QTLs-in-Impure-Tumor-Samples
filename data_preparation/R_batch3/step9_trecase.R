
library(asSeq, lib="/nas02/home/w/e/weisun/R/Rlibs/")

hinfo = help(package="asSeq")$info[[1]][4]
hinfo = unlist(strsplit(hinfo, split="\\s+", perl=TRUE))[2]
if(hinfo!="0.99.0"){ stop("wrong version") }

setwd("/lustre/scr/w/e/weisun/TCGA/data")

# -------------------------------------------------------------------------
# read the data and info
# -------------------------------------------------------------------------

datT  = read.table(file = "gene_counts_EA_filtered.txt", sep = "\t", 
header = TRUE, as.is=TRUE)

info = read.table(file = "~/research/data/human/ens_gene_loci.txt", 
sep = "\t", header = TRUE, as.is=TRUE)

dim(datT)
dim(info)

mati = match(rownames(datT), info$id)
table(is.na(mati))

info = info[mati,]

all(rownames(datT) == info$id)
datT[1:2, 1:8]
info[1:2,]

sams = names(datT)

# -------------------------------------------------------------------------
# The total number of reads per sample
# -------------------------------------------------------------------------

tot = colSums(datT)
s75 = apply(datT, 2, quantile, prob=0.75)
cor(tot, s75)

# -------------------------------------------------------------------------
# read the information covariates
# -------------------------------------------------------------------------

covs = read.table("covariates_EA.txt", sep = "\t", header = TRUE)
dim(covs)
covs[1:2,1:5]

X = t(data.matrix(covs[,-1]))
all(rownames(X) == sams)

X = cbind(log(s75+1), X)
dim(X)
X[1:4,1:5]

# -------------------------------------------------------------------------
# read AS expression data and info
# -------------------------------------------------------------------------

datA1 = read.table(file = "gene_asCounts_EA_hap1_filtered.txt", 
sep = "\t", header = TRUE)
dim(datA1)

datA1[1:2, 1:8]
all(names(datA1) == sams)

datA2 = read.table(file = "gene_asCounts_EA_hap2_filtered.txt", 
sep = "\t", header = TRUE)
dim(datA2)

datA2[1:2, 1:8]
all(names(datA2) == sams)

all(rownames(datA1) == rownames(datT))
all(rownames(datA2) == rownames(datT))

# ---------------------------------------------------------------------
# test for cis-eQTL
# ---------------------------------------------------------------------

setwd("/lustre/scr/w/e/weisun/TCGA/data/")

chr2 = sprintf("chr%d", chr1)

st1 = "\n---------------------------------------------------\n"
message(sprintf("%s chr%d, %d-%d %s", st1, chr1, id1, id2, st1))

# -------------------------------------------------------------------------
# Prepare gene expression data
# -------------------------------------------------------------------------

wchr = which(info$chr == chr2)

if(length(wchr) == 0){
  warninng("no genes at", chr2)
  next
}

Y  = t(data.matrix(datT[wchr,]))
Y1 = t(data.matrix(datA1[wchr,]))
Y2 = t(data.matrix(datA2[wchr,]))
einfo  = info[wchr,]

dim(Y)
Y[1:2,1:4]
dim(Y1)
dim(Y2)

dim(einfo)
einfo[1:2,]

# -------------------------------------------------------------------------
# Prepare genotype data
# -------------------------------------------------------------------------

datFile = sprintf("../phased_geno_EA/%s_MAF_0.02_SNP_%d_%d.txt", chr2, id1, id2)
Z       = read.table(datFile, comment.char = "", sep="\t", 
                     header = FALSE, as.is=TRUE, row.names=1)
Z       = data.matrix(Z)
dim(Z)
Z[1:2,1:5]

if(any(rownames(Y) != rownames(Z))){ stop("row names do not match\n") }

infFile = sprintf("../phased_geno_EA/%s_MAF_0.02_info.txt", chr2)
genoInf = read.table(infFile, comment.char = "", sep="\t", 
                     header = FALSE, as.is=TRUE)
dim(genoInf)
genoInf[1:2,]

genoInf = genoInf[id1:id2,]
dim(genoInf)
genoInf[1:2,]

# -------------------------------------------------------------------------
# The actual computation
# -------------------------------------------------------------------------

chrI = as.integer(chr1)
eChr = rep(chrI, ncol(Y1))
ePos = (einfo$start + einfo$end)/2

mChr = rep(chrI, ncol(Z))
mPos = genoInf$V2

ta  = trecase(Y, Y1, Y2, X, Z, p.cut=1e-3, 
              output.tag=sprintf("trecase_%s_SNP_%d-%d", chr2 ,id1, id2), 
              local.only=TRUE, eChr = eChr, ePos = ePos, 
              mChr = mChr, mPos = mPos, trace=0, maxit=100)


