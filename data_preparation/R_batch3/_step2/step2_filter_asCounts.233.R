i = 234

library(asSeq, lib="/nas02/home/w/e/weisun/R/Rlibs/")

# -------------------------------------------------------------------------
# read in the list of the SNP to be excluded
# -------------------------------------------------------------------------

setwd("/lustre/scr/w/e/weisun/TCGA/hetSNP_EA/")
files = list.files(path = ".", pattern="hetSNP_")

sams = gsub("hetSNP_", "", files)
sams = gsub(".txt", "", sams, fixed=TRUE)

#for(i in 1:length(files)){

  f1   = files[i]
  sam1 = sams[i]

  cat("\n", sam1, date(), "\n")
  
  input   = sprintf("../bam/%s_sorted_by_name_uniq_filtered.bam", sam1)
  outputTag  = sprintf("../bam/%s_asCounts_hetSNP_EA", sam1)
  snpList = f1
  
  if(! file.exists(f1)){
    stop("snpList file does not exist") 
  }
  
  extractAsReads(input, snpList, outputTag)
#}
