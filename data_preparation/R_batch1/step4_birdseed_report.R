
setwd("~/research/TCGA/BRCA/genotype_normal")

dat1 = read.table("birdseed-v2.report.txt", sep="\t", 
  header=TRUE, as.is=TRUE)
dim(dat1)
names(dat1)
dat1[1:2,]

table(dat1$computed_gender, dat1$cn.probe.chrXY.ratio_gender)
table(dat1$computed_gender, dat1$em.cluster.chrX.het.contrast_gender)

summary(dat1[,3:5])
summary(dat1[,4] + dat1[,5])

# --------------------------------------------------------------
# read in sample information
# --------------------------------------------------------------

sam1 = read.table("../info/brca_samples2use_after_qc.txt", 
  sep="\t", header=TRUE, as.is=TRUE)
dim(sam1)
sam1[1:2,]

if(any(sam1$DNAnorml_arrayFile != dat1$cel_files)){
  stop("some samples do not have affy6 data!\n")
}

dat2 = read.table("../R_batch1/apt-geno-qc.txt", sep="\t", 
  header=TRUE, as.is=TRUE)
dim(dat2)
dat2[1:2,]

mat2 = match(sam1$DNAnorml_arrayFile, dat2$cel_files)
sex2  = dat2$em.cluster.chrX.het.contrast_gender[mat2]
table(dat1$computed_gender, sex2)
table(dat1$em.cluster.chrX.het.contrast_gender, sex2)

# --------------------------------------------------------------
# plot1
# --------------------------------------------------------------

png("../figures/birdseed_report_all.png", width=5, height=5, units="in", res=200)

par(mfrow=c(1,1), mar=c(5,4,1,1), bty="n")

plot(dat1$call_rate, dat1$het_rate,  typ="n",
xlab="call rate", ylab="het rate", bty="n")

grp1 = which(dat1$computed_gender == "female")
grp2 = which(dat1$computed_gender == "male")
grp3 = which(dat1$computed_gender == "unknown")

points(dat1$call_rate[grp1], dat1$het_rate[grp1], 
cex=0.5, col="red", pch=19)

points(dat1$call_rate[grp2], dat1$het_rate[grp2], 
cex=0.5, col="darkgreen", pch=3)

points(dat1$call_rate[grp3], dat1$het_rate[grp3], 
cex=1.0, col="darkblue", pch=10)

ll = c("female", "male", "unknown")
legend("bottomleft", ll, bty="n", pch=c(19,3,10),  
col=c("red", "darkgreen", "darkblue"))


dev.off()

# --------------------------------------------------------------
# plot2
# --------------------------------------------------------------

w1   = which(dat1$computed_gender == "female")
dat1 = dat1[w1,]
sam1 = sam1[w1,]

png("../figures/birdseed_report_female.png", width=9, height=6, units="in", res=200)

par(mfrow=c(2,1), mar=c(5,4,1,1), bty="n", cex=1.08)

boxplot(dat1$call_rate ~ sam1$DNAnorml_plate, xlab="", 
        ylab="Call rate", cex=0.7, las=3)
mtext("Plate", side=1, line=4) 

boxplot(dat1$het_rate  ~ sam1$DNAnorml_plate, xlab="", 
        ylab="Het rate", cex=0.7, las=3)
mtext("Plate", side=1, line=4) 

dev.off()

sam1$gender = dat1$computed_gender

# --------------------------------------------------------------
# write out gender informaiton
# --------------------------------------------------------------

write.table(sam1, file = "../info/brca_samples2use_after_qc_female.txt", 
  append = FALSE, quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = TRUE)

