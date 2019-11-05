i = 516

library(isoform, lib.loc="/nas02/home/w/e/weisun/R/Rlibs/") 
bedFile = "/nas02/home/w/e/weisun/research/data/human/Homo_sapiens.GRCh37.66.nonoverlap.exon.bed"

setwd("/lustre/scr/w/e/weisun/TCGA/bam/")

cmd  = "ls *_asCounts_hetSNP_EA_hap1.bam"
ffs  = system(cmd, intern=TRUE)
length(ffs)
head(ffs)
sams = gsub("_asCounts_hetSNP_EA_hap1.bam", "", ffs)

sam1 = sams[i]
cat(i, sam1, date(), "\n")

bamFile = ffs[i]
outFile = sprintf("%s_asCounts_hap1.txt", sam1)

countReads(bamFile, bedFile, outFile)


bamFile = gsub("_hap1", "_hap2", ffs[i], fixed=TRUE)
outFile = sprintf("%s_asCounts_hap2.txt", sam1)

countReads(bamFile, bedFile, outFile)
