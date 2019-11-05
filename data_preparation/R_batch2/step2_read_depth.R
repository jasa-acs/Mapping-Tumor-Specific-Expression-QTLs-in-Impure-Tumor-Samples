
# -----------------------------------------------------------------
# read in sample information
# -----------------------------------------------------------------

setwd("~/research/TCGA/BRCA/info/")

sam = read.table("brca_samples2use_after_qc_female_caucasian.txt", 
header=TRUE, sep="\t", as.is=TRUE)
dim(sam)
sam[1,]

# -----------------------------------------------------------------
# read in count data
# -----------------------------------------------------------------

setwd("~/research/TCGA/BRCA/bam/_count")

files = list.files(pattern="count")
sams  = gsub("count_", "", files, fixed=TRUE)
sams  = gsub(".txt",  "", sams,  fixed=TRUE)
cts   = matrix(NA, nrow=length(sams), ncol=2)

for(i in 1:length(sams)){

  sam1 = sams[i]
  ff1  = files[i]
  
  ct   = scan(ff1, quiet=TRUE)
  cts[i,] = ct
}

cts = cts/1e6

summary((cts[,1] - cts[,2])/cts[,1] )

# -----------------------------------------------------------------
# summarize
# -----------------------------------------------------------------

setwd("~/research/TCGA/BRCA/figures")

pdf("TReC_vs_QCed_TReC.pdf", width=3, height=3)
par(mar=c(5,4,1,1))
plot(cts[,1], cts[,2], xlab="Total # of Reads (million)", 
  ylab="Total # of Reads after QC", bty="n", cex=0.5, pch=20, col="darkgrey")
w2kp = which(sams %in% sam$DNAnorml_patientID)
length(w2kp)
points(cts[w2kp,1], cts[w2kp,2], col="darkgreen", cex=0.5)

b = median(cts[,2])/median(cts[,1])
abline(0, b, col="darkred", lwd=2)
abline(0, 0.70, col="orange", lwd=2)

table(cts[w2kp,2]/cts[w2kp,1] < 0.7)

lg = c(sprintf("y=%.2fx", b), "y=0.70x")
legend("topleft", lg, lty=c(1,1), col=c("darkred", "orange"), bty="n")
dev.off()

female_caucasian = sams %in% sam$DNAnorml_patientID

db = data.frame(sample=sams, cts, round(cts[,2]/cts[,1],4), female_caucasian)

names(db) = c("sample", "TReC", "TReC_after_QC", "ratio", "femaleCau")
dim(db)
db[1:2,]

setwd("~/research/TCGA/BRCA/info")

write.table(db, file = "samples_count.txt", append = FALSE, quote = FALSE, 
  sep = "\t", row.names = FALSE, col.names = TRUE)


