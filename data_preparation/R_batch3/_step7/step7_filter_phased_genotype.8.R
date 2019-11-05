chr1 = 9

setwd("/lustre/scr/w/e/weisun/TCGA/phased_geno_EA/")

chr2 = sprintf("chr%d", chr1)

st1 = "\n---------------------------------------------------\n"
message(sprintf("%s chr: %d %s", st1, chr1, st1))

# -------------------------------------------------------------------------
# Prepare genotype data
# -------------------------------------------------------------------------

datFile = sprintf("%s.txt", chr2)
genoDat = read.table(datFile, comment.char = "", sep="\t", 
                     header = FALSE, as.is=TRUE, row.names=1)
dim(genoDat)
genoDat[1:2,1:5]

infFile = sprintf("../1000G/anno_out/%s.annotation.txt", chr2)
genoInf = read.table(infFile, comment.char = "", sep="\t", 
header = TRUE, as.is=TRUE)

dim(genoInf)
genoInf[1:2,]

if(ncol(genoDat) != nrow(genoInf)){
  stop("dimension does not match\n")
}

genoDat = data.matrix(genoDat)

w2rm  = which(rownames(genoDat) == "A15R")
w2rm

if(length(w2rm) > 0){
  genoDat = genoDat[-w2rm,]
}

nna     = colSums(is.na(genoDat))
summary(nna)
table(nna > 0)
if(any(nna > 0)){
  genoDat[is.na(genoDat)] = 0
}

tmpDat = genoDat
ww2    = tmpDat > 2
tmpDat[ww2] = genoDat[ww2] - 2

maf = 0.5*colSums(tmpDat)/nrow(tmpDat)
summary(maf)

rm(tmpDat)

# -------------------------------------------------------------------------
# keep those MAF >= 0.02
# -------------------------------------------------------------------------

wkp     = which( maf >= 0.02 & maf <= 0.98)
genoDat = genoDat[,wkp]
dim(genoDat)
genoDat[1:2,1:5]

genoInf = genoInf[wkp,]
dim(genoInf)
genoInf[1:2,]

summary(maf[wkp])

message(sprintf("%d out of %d SNPs with MAF >= 0.02 are used", 
                length(wkp), length(maf)))

write.table(genoDat, file = sprintf("%s_MAF_0.02.txt", chr2), 
            append = FALSE, quote = FALSE, sep = "\t", 
            row.names = TRUE, col.names = FALSE)

write.table(genoInf, file = sprintf("%s_MAF_0.02_info.txt", chr2), 
            append = FALSE, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)
