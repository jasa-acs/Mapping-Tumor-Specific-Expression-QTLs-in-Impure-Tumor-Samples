chrId = 5


#for(chrId in 1:22){
  chr1 = sprintf("chr%d", chrId)

  message(chr1, " ", date())

  outF = sprintf("/lustre/scr/w/e/weisun/TCGA/phased_geno_EA/chr%d.txt", chrId)

  # ----------------------------------------------------------------------
  # read in the SNP annotation 
  # ----------------------------------------------------------------------

  setwd("/lustre/scr/w/e/weisun/TCGA/1000G/anno_out/")

  annoFile = sprintf("%s.annotation.txt", chr1)
  anno1    = read.table(annoFile, header=TRUE, sep="\t", 
              comment.char="", as.is=TRUE)

  dim(anno1)
  anno1[1:5,]

  setwd("/lustre/scr/w/e/weisun/TCGA/1000G/snps/")

  snpFile = sprintf("v2.20101123.%s.snps", chr1)
  snp1    = scan(snpFile, what=character())
  snp1[1:5]

  ## some of the markers in snpFile are indels
  ## try to exclude them for now

  wsnp = match(anno1$SNP, snp1)
  all(snp1[wsnp] == anno1$SNP)

  if(any(anno1$X.chr != chr1)){
    stop("chromsome name does not match\n")
  }

  # ----------------------------------------------------------------------
  # read in the list of the SNP to be excluded
  # ----------------------------------------------------------------------

  setwd("/lustre/scr/w/e/weisun/TCGA/MACH_output_EA/")

  ff1 = sprintf("%s.hap", chr1)
  cmd = sprintf("wc -l %s", ff1)
  wc1 = system(cmd, intern=TRUE)
  wc1 = unlist(strsplit(wc1, split="\\s+"))[1]
  if(wc1 != "1102") stop("number of lines is not correct")

  for(k in 1:551){
    message("")
    message("  ", k, " ", date())
    
    ff2  = sprintf("%s.%d.%d.tmp", ff1, 2*k-1, 2*k)
    cmd1 = sprintf("sed -n '%d,%d p' %s > %s", 2*k-1, 2*k, ff1, ff2)
    system(cmd1)

    dat1 = scan(ff2, what=character())
    
    cmd1 = sprintf("rm %s", ff2)
    system(cmd1)

    if(dat1[1] != dat1[4]){ stop("sample names do not match\n") }
    
    if(dat1[2] != "HAPLO1"){ 
      stop(sprintf("expect HAPLO1 but see %s\n", dat1[2])) 
    }
    
    if(dat1[5] != "HAPLO2"){ 
      stop(sprintf("expect HAPLO2 but see %s\n", dat1[5])) 
    }
    
    hap1 = unlist(strsplit(dat1[3], split=""))
    hap2 = unlist(strsplit(dat1[6], split=""))
    
    if(length(hap1) != length(snp1)){
      stop("length of haplotype 1 does not match annotation\n")
    }
    
    if(length(hap2) != length(snp1)){
      stop("length of haplotype 2 does not match annotation\n")
    }
    
    hap1 = toupper(hap1[wsnp])
    hap2 = toupper(hap2[wsnp])
    
    ww0  = which(hap1==anno1$reference_allele & hap2==anno1$reference_allele)
    ww4  = which(hap1==anno1$observed_allele  & hap2==anno1$observed_allele)
    ww1  = which(hap1==anno1$reference_allele & hap2==anno1$observed_allele)
    ww3  = which(hap1==anno1$observed_allele  & hap2==anno1$reference_allele)

    length(ww0)
    length(ww1)
    length(ww3)
    length(ww4)
    
    geno1  = rep(NA, length(hap1))
    geno1[ww0] = 0
    geno1[ww1] = 1
    geno1[ww3] = 3
    geno1[ww4] = 4

    write.table( t(table(geno1, useNA="ifany")), quote=FALSE, row.names=FALSE, sep="\t")
    
    sam1  = unlist(strsplit(dat1[1], split="->"))[2]
    geno1 = c(sam1, as.character(geno1))
    
    cat(geno1, file = outF, sep="\t", append = (k>1))
    cat("\n", file=outF, append = TRUE)

  }
#}
