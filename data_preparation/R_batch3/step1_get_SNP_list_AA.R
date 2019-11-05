
# ----------------------------------------------------------------------
# gerate a list of heterzygous SNPs for each sample, which will be 
# used by function "asCounts"
#
# the output file is a tab-delimiated file with four columns, 
# chromosome, position, allele 1 and allele 2, without header. 
#
# ----------------------------------------------------------------------

for(chrId in 1:22){
  chr1  = sprintf("chr%d", chrId)

  message(chr1, " ", date())

  # ----------------------------------------------------------------------
  # read in the SNP annotation 
  # ----------------------------------------------------------------------

  setwd("/Volumes/Moon/TCGA_BRCA/1000G/anno_out/")

  annoFile = sprintf("%s.annotation.txt", chr1)
  anno1    = read.table(annoFile, header=TRUE, sep="\t", 
              comment.char="", as.is=TRUE)

  dim(anno1)
  anno1[1:5,]

  setwd("/Volumes/Moon/TCGA_BRCA/1000G/snps/")

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

  setwd("/Volumes/Moon/TCGA_BRCA/MACH_output_AA/")

  ff1 = sprintf("%s.hap", chr1)
  cmd = sprintf("wc -l %s", ff1)
  wc1 = system(cmd, intern=TRUE)
  wc1 = unlist(strsplit(wc1, split="\\s+"))[2]
  if(wc1 != "100") stop("number of lines is not correct")

  for(k in 1:50){
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
    
    alleles = c("A", "C", "G", "T")
    
    wdiff = which(hap1 != hap2 & hap1 %in% alleles & hap2 %in% alleles)
    
    hap1A = hap1[wdiff]
    hap1B = hap2[wdiff]
    pos1  = anno1$coordinate[wdiff]
    
    hetDat = data.frame(chr=rep(chr1, length(pos1)), pos=pos1, hap1A, hap1B)
    hetDat = hetDat[order(pos1),]
    
    sam1   = unlist(strsplit(dat1[1], split="->"))[2]
    ff2 = sprintf("../hetSNP_AA/hetSNP_%s.txt", sam1)
    
    if(chr1 == "chr1"){
      app = FALSE
    }else{
      app = TRUE
    }
    
    write.table(hetDat, file = ff2, append = app, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)

  }
}
