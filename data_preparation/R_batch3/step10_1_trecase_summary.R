
library(MASS)
library(limma)

# -------------------------------------------------------------------------
# read the gene expression information
# -------------------------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/data/")

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
# first check the frequency of eQTL
# -------------------------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/data/trecase")

freqSumTReC = freqSumASE = freqSumJoint = numeric(101)

ffs = list.files(pattern="_freq.txt")
length(ffs)
ffs[1:2]

for(ffi in ffs){
  fi = read.table(ffi, header=TRUE)
  
  f1 = as.numeric(fi[1,-1])
  freqSumTReC = freqSumTReC + f1

  f1 = as.numeric(fi[2,-1])
  freqSumASE = freqSumASE + f1

  f1 = as.numeric(fi[3,-1])
  freqSumJoint = freqSumJoint + f1
}

pdf("../trecase_figure/TReCASE_eQTL_freq.pdf", width=6, height=4)
par(mar=c(3,2,4,1))

barplot(freqSumTReC[-1], xlab="", main="TReC model")
mtext(sprintf("%d P-values", sum(freqSumTReC[-1])), side=1, line=0.5)

barplot(freqSumASE[-1], xlab="", main="ASE model")
mtext(sprintf("%d P-values", sum(freqSumASE[-1])), side=1, line=0.5)

barplot(freqSumJoint[-1], xlab="", main="TReCASE model")
mtext(sprintf("%d P-values", sum(freqSumJoint[-1])), side=1, line=0.5)

dev.off()

# -------------------------------------------------------------------------
# read in the eQTL results
# -------------------------------------------------------------------------

eqtls = NULL

ffs = list.files(pattern="_eqtl.txt")
pcutNew = 1e-3

kk = 0
for(ffi in ffs){
  kk   = kk + 1
  if(kk %% 100 == 0) { cat(kk, date(), "\n") }
  
  fnam = unlist(strsplit(ffi, split="_"))
  chri = gsub("chr", "", fnam[2])
  if(chri=="X"){ chri="23" }
  if(chri=="Y"){ chri="24" }
  if(chri=="M" || chri=="Mt"){ chri="25" }
  chri = as.numeric(chri)
  
  posi = as.numeric(unlist(strsplit(fnam[4], split="-")))
  ei   = read.table(ffi, header=TRUE)
  ei[["Chr"]] = rep(chri, nrow(ei))
  ei$MarkerRowID = ei$MarkerRowID + posi[1] - 1

  ei = ei[which(ei$TReC_Pvalue < pcutNew | ei$ASE_Pvalue < pcutNew | ei$final_Pvalue < pcutNew),]
  eqtls = rbind(eqtls, ei)
}

dim(eqtls)
eqtls[1:5,]

summary(eqtls$ASE_Pvalue)
summary(eqtls$TReC_Pvalue)
summary(eqtls$final_Pvalue)

eqtls = eqtls[which(!is.na(eqtls$final_Pvalue)),]

summary(eqtls$ASE_Pvalue)
summary(eqtls$TReC_Pvalue)
summary(eqtls$final_Pvalue)

# -------------------------------------------------------------------------
# looping
# -------------------------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/data/")

uchrs = c(1:22)

eqtls[["gene"]]  = rep(NA, nrow(eqtls))
eqtls[["start"]] = rep(NA, nrow(eqtls))
eqtls[["end"]]   = rep(NA, nrow(eqtls))
eqtls[["snp"]]   = rep(NA, nrow(eqtls))

for(chr1 in uchrs){
  
  chr2 = sprintf("chr%d", chr1)
  
  st1 = "\n---------------------------------------------------\n"
  message(sprintf("%s %s, %s %s", st1, chr2, date(), st1))  
  
  # -------------------------------------------------------------------------
  # Prepare gene expression info
  # -------------------------------------------------------------------------
  
  wchr = which(info$chr == chr2)
  
  if(length(wchr) == 0){
    warninng("no genes at", chr2)
    next
  }
  
  einfo  = info[wchr,]
  dim(einfo)
  einfo[1:2,]
  
  # -------------------------------------------------------------------
  # read in SNP info
  # -------------------------------------------------------------------
  
  datF  = sprintf("../phased_geno_EA/%s_MAF_0.02.txt", chr2)
  infF  = sprintf("../phased_geno_EA/%s_MAF_0.02_info.txt", chr2)
  ginfo = read.table(infF, comment.char = "", sep="\t", 
                     header = FALSE, as.is=TRUE)
  dim(ginfo)
  ginfo[1:2,]
  
  wchr = which(eqtls$Chr == chr1)
  geneRowIDs = as.numeric(eqtls[["GeneRowID"]][wchr])
  
  eqtls[["gene"]][wchr]  = einfo$id[geneRowIDs]
  eqtls[["start"]][wchr] = einfo$start[geneRowIDs]
  eqtls[["end"]][wchr]   = einfo$end[geneRowIDs]
  eqtls[["snp"]][wchr]   = ginfo$V5[eqtls$MarkerRowID[wchr]]  
}

dim(eqtls)
eqtls[1:2,]

write.table(eqtls, file = "eqtls_1e-3.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# -------------------------------------------------------------------------
# trec p-value vs. final p-value
# -------------------------------------------------------------------------

which.min(eqtls$TReC_Pvalue)
eqtls[which.min(eqtls$TReC_Pvalue),]

w1 = which(eqtls$trans_Pvalue < 0.05)
w2 = which(eqtls$trans_Pvalue >= 0.05)

dim(eqtls)
length(w1)
length(w2)

png("trecase_figure/trecP_vs_finalP.png", width=6.5, height=3.5, 
  units="in", res=300)
par(mfrow=c(1,2), mar=c(5,4,4,1), bty="n")

plot(-log10(eqtls$TReC_Pvalue[w1]), -log10(eqtls$final_Pvalue[w1]), 
  main="trans-eQTL", xlab="TReC p-value", ylab="Final p-value", cex=0.5, 
  col=rgb(0, 100/256, 0, 0.3))
abline(0, 1, lty=2, col="red")

plot(-log10(eqtls$TReC_Pvalue[w2]), -log10(eqtls$final_Pvalue[w2]),
  main="cis-eQTL", xlab="TReC p-value", ylab="Final p-value", cex=0.5,
  col=rgb(0, 100/256, 0, 0.3))
abline(0, 1, lty=2, col="red")
dev.off()

# -------------------------------------------------------------------------
# keep the most significant eQTL per gene
# -------------------------------------------------------------------------

geneIDs = paste(eqtls$Chr, eqtls$GeneRowID, sep=":") 
length(geneIDs)

eqtls   = eqtls[order(geneIDs),]
geneIDs = paste(eqtls$Chr, eqtls$GeneRowID, sep=":") 
ugeneID = unique(geneIDs)
length(ugeneID)
length(geneIDs)

tmin = tapply(eqtls$final_Pvalue, geneIDs, which.min)
tt1  = table(geneIDs)
csum = cumsum(tt1)

if(any(ugeneID != names(tt1))){ stop("ugeneID != names(tt1)\n") }
if(any(ugeneID != names(tmin))){ stop("ugeneID == names(tmin)\n") }
if(any(ugeneID != names(csum))){ stop("ugeneID == names(csum)\n") }

tt1[1:5]
csum[1:5]
tmin[1:5]
ugeneID[1:5]

w2kp = tmin + c(0, csum[-length(csum)])
geneFinal = eqtls[w2kp,]
geneFinal = geneFinal[order(geneFinal$Chr, geneFinal$GeneRowID),]
geneFinal[1:5,]

write.table(geneFinal, file = "eqtls_1e-3_per_gene.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# -------------------------------------------------------------------------
# for different p-value cutoffs of pvalTReC, check the proportion of 
# cis-eQTL
# -------------------------------------------------------------------------

pcuts = c(1e-5, 1e-7, 1e-10, 1e-15)

nAll = nCis = nTrans = rep(NA, length(pcuts))

for(i in 1:length(pcuts)){
  pcut1     = pcuts[i]
  nAll[i]   = length(which(geneFinal$TReC_Pvalue < pcut1))
  nCis[i]   = length(which(geneFinal$TReC_Pvalue < pcut1 & geneFinal$trans_Pvalue >= 0.05))
  nTrans[i] = length(which(geneFinal$TReC_Pvalue < pcut1 & geneFinal$trans_Pvalue < 0.05))
}

nTest   = nCis+nTrans
cisProp = round(nCis/nTest,3)
cisP    = rbind(pcuts, nAll, nTest, nCis, cisProp)

write.table(cisP, sep=" & ", quote=FALSE, col.names=FALSE)

# -------------------------------------------------------------------------
# try another cutoff of transPval 
# -------------------------------------------------------------------------

nAll = nCis = nTrans = rep(NA, length(pcuts))

for(i in 1:length(pcuts)){
  pcut1     = pcuts[i]
  nAll[i]   = length(which(geneFinal$TReC_Pvalue < pcut1))
  nCis[i]   = length(which(geneFinal$TReC_Pvalue < pcut1 & geneFinal$trans_Pvalue >= 0.001))
  nTrans[i] = length(which(geneFinal$TReC_Pvalue < pcut1 & geneFinal$trans_Pvalue <  0.001))
}

nTest   = nCis+nTrans
cisProp = round(nCis/nTest,3)
cisP    = rbind(pcuts, nAll, nTest, nCis, cisProp)

write.table(cisP, sep=" & ", quote=FALSE, col.names=FALSE)

# -------------------------------------------------------------------------
# plot Venn Diagram
# -------------------------------------------------------------------------

geneIDs = paste(eqtls$Chr, eqtls$GeneRowID, sep=":") 

pminTReC = tapply(eqtls$TReC_Pvalue,  geneIDs, min)
pminASE  = tapply(eqtls$ASE_Pvalue,   geneIDs, min)
pminJoin = tapply(eqtls$final_Pvalue, geneIDs, min)

is.numeric(pminTReC)
is.numeric(pminASE)
is.numeric(pminJoin)

length(pminTReC)
length(pminASE)
length(pminJoin)

png("trecase_figure/vennDiag.png", width=6, height=6, units="in", res=200)
par(mfrow=c(2,2))

pvals = cbind(pminTReC, pminASE, pminJoin)
pvals[is.na(pvals)] = 1.0
colnames(pvals) = c("TReC", "ASE", "Joint")

a1 = data.matrix(pvals <= 1e-6)
b1 = vennCounts(a1)
vennDiagram(b1, mar=rep(0,4), cex=1)
mtext("P-value <= 1e-6", line=2, font=2)

a1 = data.matrix(pvals <= 1e-8)
b1 = vennCounts(a1)
vennDiagram(b1, mar=rep(0,4), cex=1)
mtext("P-value <= 1e-8", line=2, font=2)

a1 = data.matrix(pvals <= 1e-10)
b1 = vennCounts(a1)
vennDiagram(b1, mar=rep(0,4), cex=1)
mtext("P-value <= 1e-10", line=2, font=2)

a1 = data.matrix(pvals <= 1e-15)
b1 = vennCounts(a1)
vennDiagram(b1, mar=rep(0,4), cex=1)
mtext("P-value <= 1e-15", line=2, font=2)

dev.off()
