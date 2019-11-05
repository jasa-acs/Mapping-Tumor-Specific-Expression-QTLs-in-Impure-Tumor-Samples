
setwd("/lustre/scr/w/e/weisun/TCGA/bam/")

load("~/research/data/human/Homo_sapiens.GRCh37.66.nTE.RData")
dim(nTE)
nTE[1:2,]
length(unique(nTE$geneId))

ffs  = list.files(pattern="_asCounts_hap1.txt")
nn   = length(ffs)
ffs[1:2]
nn

sams  = gsub("_asCounts_hap1.txt", "", ffs)
couts = matrix(0, nrow=nrow(nTE), ncol=length(sams))
colnames(couts) = sams
rownames(couts) = nTE$geneId

for(idx in 1:length(sams)){
  
  cat(idx, date(), "\n")
  
  f1   = ffs[idx]
  dat  = scan(f1, what=character(0))  
  dat  = matrix(dat, ncol=2, byrow=TRUE);
 
  colNames = c("count", "exons")
  cN = sprintf("%s and %s", colNames[1], colNames[2])
  
  if(ncol(dat) != 2){
    stop(countFile, " should have 2 columns: ", cN, "\n")
  }
  
  colnames(dat) = colNames
  dim(dat)
  dat[1:2,]
  
  # --------------------------------------------------------- 
  # obtain transcript cluster IDs and gene IDs
  # --------------------------------------------------------- 
  
  groupIDs = strsplit(dat[,"exons"], split=";|\\|", perl=TRUE)
  
  splitFun <- function(vx){
    unique(vx[grep("ENSG", vx)])
  }
  
  date()
  geneIDs = lapply(groupIDs, splitFun)
  date()
  geneIDs[1:2]
  
  # --------------------------------------------------------- 
  # if only one geneID is shared by all the geneIDs, 
  # assign the count to this geneID, otherwise ignore it.
  # --------------------------------------------------------- 
  
  ngIDs   = sapply(geneIDs, length)
  table(ngIDs)
  
  w2check = which(ngIDs > 1)
  
  chkGIDs <- function(g1){
    gkp    = ""
    ncombo = sum(!grepl(":", g1, fixed=TRUE))
    
    if(ncombo <= 1){
      g2s = strsplit(g1, split=":", fixed=TRUE)
      gus = unique(unlist(g2s))
      foundONE = FALSE
      
      for(gu1 in gus){
        if (all(grepl(gu1, g1))){
          if(foundONE){  
            gkp = ""
            break
          }else{
            foundONE = TRUE
            gkp = gu1
          }
        }
      }
    }
    
    gkp
  }

  gIDchk  = sapply(geneIDs[w2check], chkGIDs)
  length(gIDchk)
  gIDchk[1:4]
  
  geneIDs[w2check] = gIDchk
  n1      = length(geneIDs)
  geneIDs = unlist(geneIDs)
  if(n1 != length(geneIDs)){ stop("non-unique geneIDs\n") }
  
  gID2rm = w2check[which(gIDchk=="")]
  str1   = "combinations are skipped because they belong to different genes"
  message(length(gID2rm), " exon ", str1)
  
  dim(dat)
  if(length(gID2rm) > 0){
    dat     = dat[-gID2rm,]
    geneIDs = geneIDs[-gID2rm]
  }
  
  dim(dat)
  length(unique(geneIDs))
  
  # --------------------------------------------------------- 
  # record the counts
  # --------------------------------------------------------- 
    
  cts     = as.numeric(dat[,"count"])
  geneCts = tapply(cts, geneIDs, sum)
  
  mat1    = match(names(geneCts), nTE$geneId)
  wNotNA  = which(! is.na(mat1))
  
  pp1 = round(sum(geneCts[-wNotNA])/sum(geneCts),4)
  nn1 = length(geneCts) - length(wNotNA)
  message(100*pp1, "% of reads @ ", nn1, " gene combinations are skipped\n")

  couts[mat1[wNotNA],idx] = geneCts[wNotNA]
}

outF = "../data/gene_asCounts_EA_hap1.txt"

write.table(couts, file = outF, append = FALSE, quote = FALSE, sep = "\t",
row.names = TRUE, col.names = TRUE)
