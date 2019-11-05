
library(MASS)
library(limma)

source("~/research/R/asSeq/_test/tests.R")
source("~/research/R/asSeq/_test/aseR.R")
source("~/research/R/asSeq/_test/trecR.R")
source("~/research/R/asSeq/_test/trecaseR.R")
source("~/research/R/asSeq/_test/plots.R")

# -------------------------------------------------------------------------
# read the data and info
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

# -------------------------------------------------------------------------
# read in the eQTL results
# -------------------------------------------------------------------------

eqtls = read.table("eqtls_5e-5_per_gene.txt", header=TRUE, sep="\t", as.is=TRUE)
dim(eqtls)
eqtls[1:2,]

# -------------------------------------------------------------------------
# looping
# -------------------------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/data/")

uchrs = c(1:22)

for(chr1 in uchrs){
  
  chr2 = sprintf("chr%d", chr1)
  
  st1 = "\n---------------------------------------------------\n"
  message(sprintf("%s %s, %s %s", st1, chr2, date(), st1))  
  
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
  
  # -------------------------------------------------------------------
  # read in info
  # -------------------------------------------------------------------
  
  datF  = sprintf("../phased_geno_EA/%s_MAF_0.02.txt", chr2)
  infF  = sprintf("../phased_geno_EA/%s_MAF_0.02_info.txt", chr2)
  ginfo = read.table(infF, comment.char = "", sep="\t", 
                     header = FALSE, as.is=TRUE)
  dim(ginfo)
  ginfo[1:2,]

  # -------------------------------------------------------------------
  # randomly choose 10 genes to plot in each case
  # -------------------------------------------------------------------
  
  aseEQ = eqtls[eqtls$Chr==chr1,]
  aseEQ = aseEQ[order(aseEQ$final_Pvalue),]
  aseEQ[1:10,]
  
  for(kk in 1:10){
    
    eID = aseEQ$GeneRowID[kk]
    gID = aseEQ$MarkerRowID[kk]
    
    ff2  = paste(datF, ".tmp1234", sep="")
    cmd1 = sprintf("cut -f %d %s > %s", gID+1, datF, ff2)
    system(cmd1)
    x1   = scan(ff2)
    cmd1 = sprintf("rm %s", ff2)
    system(cmd1)
    
    y0  = Y[,eID]
    y1  = Y1[,eID]
    y2  = Y2[,eID]

    tag   = "trecase_figure/x"
    sName = aseEQ$snp[kk]
    eName = aseEQ$gen[kk]
    
    fname = paste(tag, chr2, eName, sName, sep="_")
    fname = paste(fname, ".png", sep="")  
    xlab  = sName
    ylab  = eName
    
    mm1   = sprintf("pJ=%.1e,b(T/A)=%.1f/%.1f,pCis=%.0e", aseEQ$final_Pvalue[kk],
                    aseEQ$TReC_b[kk], aseEQ$ASE_b[kk], aseEQ$trans_Pvalue[kk])

    png(fname, width=8.4, height=3, units="in", res=200)
    
    layout(matrix(1:3, nrow=1), widths = c(3,3,2.4))
    par(mar=c(4,4,3,1), bty="n", cex=0.85, cex.main=0.95)

    plotTReC(y0, x1, X, xlab, ylab)
    
    if(is.na(aseEQ$ASE_Pvalue[kk])){
      plot(0:1, 0:1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
    }else{
      plotASE(y1, y2, x1, xlab, ylab, sub=mm1, min.nTotal=10)
    }
    
    dev.off()
  }
  
}

