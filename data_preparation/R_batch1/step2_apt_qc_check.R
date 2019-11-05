
plotIt <- function(dat1, grps){

  par(mar=c(5,4,1,1), bty="n")
  plot(dat1$qc.call.rate.all, dat1$contrast.qc, cex=0.5, pch=20, 
       xlab="QC call rate", ylab="contrast QC", bty="n", type="n")
  abline(v=0.80, col="grey")
  abline(v=0.90, col="grey")
  abline(h=0.40, col="grey")
  
  
  t1 = table(grps)
  n1 = length(t1)
  u1 = names(t1)
  
  cols = rainbow(n1)
  
  if(n1 < 6){
    pnts = c(1, 15:(14+n1))
  }else{
    pnts = 1:n1
  }
  
  for(i in 1:length(u1)){
    ui   = u1[i]
    grp1 = which(as.character(grps) == ui)
    points(dat1$qc.call.rate.all[grp1], dat1$contrast.qc[grp1], 
           cex=0.8, col=cols[i], pch=pnts[i])
  }
  
  u2 = paste(u1, "(", t1, ")", sep="")
  legend("topleft", u2, col=cols, pch=pnts, bty="n")
  
}

# --------------------------------------------------------------
# read in QC results
# --------------------------------------------------------------

setwd("~/research/TCGA/BRCA/R_batch1")

dat1 = read.table("apt-geno-qc.txt", sep="\t", header=TRUE, as.is=TRUE)
dim(dat1)
dat1[1:2,]

sex = dat1$em.cluster.chrX.het.contrast_gender
table(sex)

# --------------------------------------------------------------
# read in sample information
# --------------------------------------------------------------

sam1 = read.table("../info/brca_samples2use.txt", sep="\t", header=TRUE, as.is=TRUE)
dim(sam1)
sam1[1:2,]

if(!all(dat1$cel_files == sam1$DNAnorml_arrayFile)){
  stop("sample mismatch!\n")  
}

table(sam1$DNAnorml_plate)
table(sam1$DNAnorml_institution)
table(sam1$DNAnorml_type)
table(sam1$DNAnorml_portion)

summary(dat1$contrast.qc)

anova(lm(dat1$contrast.qc ~ DNAnorml_plate + DNAnorml_institution 
  + DNAnorml_type + DNAnorml_portion + sex, data=sam1))

anova(lm(dat1$qc.call.rate.all ~ DNAnorml_plate + DNAnorml_institution 
  + DNAnorml_type + DNAnorml_portion + sex, data=sam1))

chisq.test(sam1$DNAnorml_institution, sam1$DNAnorml_plate)

# --------------------------------------------------------------
# http://www.affymetrix.com/support/help/faqs/
#          genotyping_console/quality_control/faq_2.jsp
# --------------------------------------------------------------

# Contrast QC measures how well experiments resolve SNP signals 
# into three genotype clusters. It uses a subset of probes and 
# measures the differences in contrast distributions for homozygote 
# and heterozygote genotypes. In high-quality data sets, the 
# homozygote distributions are well-resolved from the heterozygote 
# distribution and the Contrast QC metric is higher than the 0.4 
# 
# --------------------------------------------------------------
# http://www.affymetrix.com/support/help/faqs/
#          genotyping_console/quality_control/faq_4.jsp
# --------------------------------------------------------------

# The QC Call Rate threshold is 86 percent or greater for the SNP 
# Array 5.0, 93 percent or greater for the 500K Array Set, and 95 
# percent or greater for the 100K Set. In good-quality data sets, 
# 90 percent of samples should pass the QC Call Rate threshold 
# and the average QC Call Rate should be in the mid-90 percent 
# range. Occasionally, poor samples will pass the QC Call Rate 
# metric, which is why Contrast QC is recommended for the SNP 
# Array 6.0. To minimize the risk of poor samples passing the 
# QC Call Rate metric, samples with outlier values for 
# heterozygosity or Birdseed call rate should be rejected.


png("../figures/QC_check_sex.png", width=6, height=6, units="in", res=400)
plotIt(dat1, sex)
dev.off()

plate = sam1$DNAnorml_plate
png("../figures/QC_check_plate.png", width=6, height=6, units="in", res=400)
plotIt(dat1, plate)
dev.off()

dtype = sam1$DNAnorml_type
png("../figures/QC_check_type.png", width=6, height=6, units="in", res=400)
plotIt(dat1, dtype)
dev.off()

# --------------------------------------------------------------
# output the files after QC filtering 
# --------------------------------------------------------------

setwd("~/research/TCGA/BRCA/info")

table(dat1$contrast.qc > 0.4, dat1$qc.call.rate.all > 0.8)
w2kp = which(dat1$contrast.qc > 0.4 & dat1$qc.call.rate.all > 0.8)
length(w2kp)
db   = sam1[w2kp,]

write.table(db, file = "brca_samples2use_after_qc.txt", append = FALSE, 
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

fname = "cel_files_normal_after_qc.txt"
cat("cel_files\n", file=fname)
cels = paste("/Volumes/Moon/TCGA_BRCA/data_DNA/", db$DNAnorml_arrayFile, sep="")
cat(cels, file=fname, sep="\n", append=TRUE)


fname = "cel_files_tumor_after_qc.txt"
cat("cel_files\n", file=fname)
cels = paste("/Volumes/Moon/TCGA_BRCA/data_DNA/", db$DNAtumor_arrayFile, sep="")
cat(cels, file=fname, sep="\n", append=TRUE)

# --------------------------------------------------------------------
# prepare wget command
# --------------------------------------------------------------------

rna2download = db$RNA_rnaseqID.Tumor

cmd = character(length(rna2download))
adr = "https://webshare.bioinf.unc.edu/seqware/v2_pipeline/"

for(i in 1:length(rna2download)){
  ri = rna2download[i]
  ci = "wget --user=weisun --password=JSM@sd12! -O "
  ci = sprintf("%s %ssorted_genome_alignments.bam", ci, ri)
  ci = sprintf("%s %s%s/sorted_genome_alignments.bam", ci, adr, ri)
  cmd[i] = ci
}

cat(cmd, file="wget.sh", sep="\n")

