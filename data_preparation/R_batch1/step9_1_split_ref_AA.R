
setwd("/lustre/scr/w/e/weisun/TCGA/1000G")

#-------------------------------------------------------------------------------
# splitRef, which can be found here 
# http://www.sph.umich.edu/csg/yli/splitRef/
# the 1000G reference data were downloaded from 
# ftp://share.sph.umich.edu/1000genomes/fullProject/2012.02.14/
#   v2.20101123.EUR.autosomes.hap.tgz
#   v2.20101123.autosomes.snps.tgz
#-------------------------------------------------------------------------------

bsubF = "ALL_submit_jobs.sh"
cat("", file=bsubF)

for (chr in 1:22){
  cmd = "software/splitRef.V002/splitRef_V002/splitRef.pl"
  cmd = sprintf("%s -windowSize 100000 -overlapSize 10000", cmd)
  cmd = sprintf("%s -hap hap/v2.20101123.chr%d.hap.gz", cmd, chr)
  cmd = sprintf("%s -snps snps/v2.20101123.chr%d.snps", cmd, chr)
  cmd = sprintf("%s -o ../data_AA/chr%d  > ../data_AA/chr%d.log\n", cmd, chr, chr)
  ff1 = sprintf("ALL_splitRef_chr%d.sh", chr)
  cat(cmd, file=ff1)
  cat(sprintf("bsub -q week bash %s\n", ff1), file=bsubF, append=TRUE)
}
