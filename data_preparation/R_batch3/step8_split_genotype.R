
setwd("/lustre/scr/w/e/weisun/TCGA/phased_geno_EA/")

chrs  = integer(0)
idx1s = idx2s = integer(0)

step  = 10000

for(chr1 in 1:22){
  
  cat(chr1, date(), "\n")
  ff1 = sprintf("chr%d_MAF_0.02.txt", chr1)
  dat = scan(ff1, what=character(0), sep="\t", nlines=1)

  nn   = length(dat) - 1
  idx1 = seq(1, nn, by=step)
  idx2 = seq(step, nn, by=step)
  kk   = length(idx1)
  
  if(length(idx2) == kk - 1){
    idx2 = c(idx2, as.numeric(nn))
  }

  for(j in 1:kk){
    ffj = sprintf("chr%d_MAF_0.02_SNP_%d_%d.txt", chr1, idx1[j], idx2[j])
    cmd = sprintf("cut -f 1,%d-%d %s > %s", idx1[j]+1, idx2[j]+1, ff1, ffj)
    system(cmd)
  }
  
  chrs  = c(chrs, rep(chr1, kk))
  idx1s = c(idx1s, idx1)
  idx2s = c(idx2s, idx2)
  
}

length(chrs)
length(idx1s)
length(idx2s)

st1 = paste(chrs,  collapse=" ")
st2 = paste(as.integer(idx1s), collapse=" ")
st3 = paste(as.integer(idx2s), collapse=" ")

message(st1)
message(st2)
message(st3)
