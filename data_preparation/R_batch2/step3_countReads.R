
library(isoform, lib.loc="/nas02/home/w/e/weisun/R/Rlibs/") 
bedFile = "/nas02/home/w/e/weisun/research/data/human/Homo_sapiens.GRCh37.66.nonoverlap.exon.bed"

setwd("/lustre/scr/w/e/weisun/TCGA/bam/")

cmd  = "ls *_sorted_by_name_uniq_filtered.bam"
ffs  = system(cmd, intern=TRUE)
length(ffs)
head(ffs)
sams = gsub("_sorted_by_name_uniq_filtered.bam", "", ffs)

sam1 = sams[i]
cat(i, sam1, date(), "\n")

bamFile = ffs[i]
outFile = sprintf("%s_counts.txt", sam1)

countReads(bamFile, bedFile, outFile)
