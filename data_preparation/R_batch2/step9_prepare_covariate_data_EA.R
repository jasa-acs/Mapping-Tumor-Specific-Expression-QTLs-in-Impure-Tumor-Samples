
##
## take 0 gene expression PCs, plates and institution, and 6 genotype PCs
## 

# --------------------------------------------------------
# read in sample information
# --------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/info/")

sam = read.table("brca_samples2use_after_qc_female.txt", 
header=TRUE, sep="\t", as.is=TRUE)
dim(sam)
sam[1,]

pam50 = read.table("BRCA.819.DataFreeze.20120912_pam50scores.txt", 
header=TRUE, sep="\t", as.is=TRUE)
dim(pam50)
pam50[1,]

all(sam$RNA_RNAseq.Tumor %in% pam50$bcr_patient_barcode)

# --------------------------------------------------------
# read in clinical information
# --------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/clinical_brca/")

clin = read.table("clinical_patient_brca.txt", quote="",
header=TRUE, sep="\t", as.is=TRUE, na.string="[Not Available]")
dim(clin)
names(clin)
all(sam$DNAnorml_patientID %in% clin$patient_id)

# --------------------------------------------------------
# check gene expression PCs
# --------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/data/")

PCs = read.table("PCs_log_TReC_EA.txt", sep="\t", header=TRUE)
dim(PCs)
PCs[1:2,1:5]

PCs$id = paste("expression", PCs$id, sep="")

# --------------------------------------------------------
# take a subset of sample information
# --------------------------------------------------------

sam1 = sam[match(names(PCs)[-1], sam$DNAnorml_patientID),]
dim(sam1)
all(names(PCs)[-1] == sam1$DNAnorml_patientID)

table(sam1$RNA_institution, useNA="ifany")
table(sam1$RNA_type, useNA="ifany")
table(sam1$RNA_plate, useNA="ifany")
table(sam1$gender, useNA="ifany")

pam50 = pam50[match(sam1$RNA_RNAseq.Tumor, pam50$bcr_patient_barcode),]
dim(pam50)
all(sam1$RNA_RNAseq.Tumor == pam50$bcr_patient_barcode)

table(pam50$tissue, useNA="ifany")
table(pam50$Call, useNA="ifany")
table(pam50$Confidence==1)
summary(pam50$Confidence[pam50$Confidence < 1])

clin  = clin[match(sam1$DNAnorml_patientID, clin$patient_id),]
dim(clin)
all(sam1$DNAnorml_patientID == clin$patient_id)


table(!is.na(clin$days_to_death), !is.na(clin$days_to_last_known_alive), useNA="ifany")
table(clin$tumor_tissue_site, useNA="ifany")
table(clin$race, useNA="ifany")
table(clin$gender, useNA="ifany")
table(sam1$gender, clin$gender, useNA="ifany")


table(clin$breast_carcinoma_estrogen_receptor_status, useNA="ifany")
table(clin$vital_status, useNA="ifany")
summary(clin$age_at_initial_pathologic_diagnosis)
table(sam1$RNA_institution, clin$tissue_source_site, useNA="ifany")


table(clin$ajcc_tumor_stage_code, useNA="ifany")
table(clin$ajcc_neoplasm_disease_stage, useNA="ifany")
table(clin$ajcc_neoplasm_disease_stage, clin$ajcc_tumor_stage_code, useNA="ifany")
table(clin$ajcc_cancer_metastasis_stage_code, clin$ajcc_neoplasm_disease_stage, useNA="ifany")
table(clin$tissue_prospective_collection_indicator, useNA="ifany")


# --------------------------------------------------------
# check the relation between PCs and covariates
# --------------------------------------------------------

pvals = matrix(NA, nrow=20, ncol=7)

for(i in 1:20){
  PCi = as.numeric(PCs[i,-1])
  ai  = anova(lm(PCi ~  as.factor(sam1$RNA_institution) 
                 + as.factor(sam1$RNA_plate)
                 + as.factor(clin$ajcc_neoplasm_disease_stage)
                 + as.factor(clin$ajcc_tumor_stage_code)
                 + as.factor(clin$breast_carcinoma_estrogen_receptor_status)
                 + clin$age_at_initial_pathologic_diagnosis 
                 + as.factor(pam50$Call)
                 ))
  pvals[i,] = ai$Pr[1:7]
}

colnames(pvals) = c("institution", "plate", "disease_stage", "stage_code", "ER", "age", "pam50")
signif(pvals,2)

table(sam1$RNA_plate, pam50$Call)
table(sam1$RNA_institution, pam50$Call)

chisq.test(sam1$RNA_plate, pam50$Call)
chisq.test(sam1$RNA_institution, pam50$Call)
table(sam1$RNA_institution, useNA="ifany")
table(sam1$RNA_plate, useNA="ifany")

table(sam1$RNA_institution, sam1$RNA_plate)

# --------------------------------------------------------
# re-do PCA
# --------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/data/")

pDat = read.table("log_TReC_EA.txt", sep = "\t", header = TRUE, as.is=TRUE)
dim(pDat)
pDat[1:2,1:5]
nDat = data.matrix(pDat[,-1])

mm = model.matrix( ~  as.factor(sam1$RNA_institution) 
+ as.factor(sam1$RNA_plate) + clin$age_at_initial_pathologic_diagnosis)
dim(mm)
s1 = svd(mm)
s1$d

H = solve(t(mm) %*% mm) 
H = mm %*% H %*% t(mm)
dim(H)

datR14Pr = nDat - nDat %*% H
summary(rowMeans(datR14Pr))

datR14Pr[is.na(datR14Pr)] = 0 
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]
PC3 =  prdatR1$vectors[,3]

# --------------------------------------------------------
# plot PCA
# --------------------------------------------------------

table(pam50$Call, clin$breast_carcinoma_estrogen_receptor_status)

setwd("/Volumes/Moon/TCGA_BRCA/figures/")

pdf("new_PCs_log_TReC_EA_rm_institution_plate.pdf", width=6, height=6)

par(mar=c(5,4,1,1), mfrow=c(2,2))
barplot(prdatR1$values[1:20], main="", 
xlab="Index", ylab="Eigen-value")

types = unique(pam50$Call)
cols  = rainbow(length(types))

legend("topright", bty="n", legend=types, col=cols, pch=1)


par(mar=c(5,4,1,1))
plot(PC1, PC2,  bty="n", type="n")
j = 0 
for(t1 in types){
  j  = j + 1
  w1 = which(pam50$Call == t1)
  points(PC1[w1], PC2[w1], col=cols[j])
}

plot(PC1, PC3,  bty="n")
j = 0 
for(t1 in types){
  j  = j + 1
  w1 = which(pam50$Call == t1)
  points(PC1[w1], PC3[w1], col=cols[j])
}

plot(PC2, PC3,  bty="n")
j = 0 
for(t1 in types){
  j  = j + 1
  w1 = which(pam50$Call == t1)
  points(PC2[w1], PC3[w1], col=cols[j])
}

dev.off()

PCs = t(prdatR1$vectors)
PCs = data.frame(id=paste("PC", 1:nrow(PCs), sep=""), PCs)
names(PCs) = names(pDat)

dim(PCs)
PCs[1:2,1:5]

write.table(PCs, file = "../data/new_PCs_log_TReC_EA_rm_institution_plate.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

pvals = matrix(NA, nrow=20, ncol=7)

for(i in 1:20){
  PCi = as.numeric(PCs[i,-1])
  ai  = anova(lm(PCi ~  as.factor(sam1$RNA_institution) 
                 + as.factor(sam1$RNA_plate)
                 + as.factor(clin$ajcc_neoplasm_disease_stage)
                 + as.factor(clin$ajcc_tumor_stage_code)
                 + as.factor(clin$breast_carcinoma_estrogen_receptor_status)
                 + clin$age_at_initial_pathologic_diagnosis 
                 + as.factor(pam50$Call)
                 ))
  pvals[i,] = ai$Pr[1:7]
}

colnames(pvals) = c("institution", "plate", "disease_stage", "stage_code", "ER", "age", "pam50")
signif(pvals,2)


# --------------------------------------------------------
# check genotype PCs
# --------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/data/")

PCsDNA = read.table("PCs_genotype_EA.txt", sep="\t", header=TRUE)
dim(PCsDNA)
PCsDNA[1:2,1:5]

PCsDNA$id = paste("genotype", PCsDNA$id, sep="")

sam1 = sam[match(names(PCsDNA)[-1], sam$DNAnorml_patientID),]
dim(sam1)
all(names(PCs)[-1] == sam1$DNAnorml_patientID)

table(sam1$DNAnorml_plate)
table(sam1$DNAnorml_institution)
table(sam1$DNAnorml_type)
table(sam1$DNAnorml_portion)

pvals = matrix(NA, nrow=20, ncol=7)

for(i in 1:20){
  PCi = as.numeric(PCsDNA[i,-1])
  ai  = anova(lm(PCi ~  as.factor(sam1$DNAnorml_institution) 
                 + as.factor(sam1$DNAnorml_plate)
                 + as.factor(clin$ajcc_neoplasm_disease_stage)
                 + as.factor(clin$ajcc_tumor_stage_code)
                 + as.factor(clin$breast_carcinoma_estrogen_receptor_status)
                 + clin$age_at_initial_pathologic_diagnosis 
                 + as.factor(pam50$Call)
                 ))
                 
  pvals[i,] = ai$Pr[1:7]
  
}

colnames(pvals) = c("institution", "plate", "disease_stage", "stage_code", "ER", "age", "pam50")
signif(pvals,2)

table(sam1$RNA_institution, sam1$DNAnorml_institution)
table(sam1$RNA_plate, sam1$DNAnorml_plate)


# --------------------------------------------------------
# write out results
# --------------------------------------------------------

mm1 = t(mm[,-1])
dim(mm1)

nms = gsub("as.factor(sam1$RNA_institution)", "", rownames(mm1), fixed=TRUE)
nms = gsub("as.factor(sam1$RNA_plate)", "", nms, fixed=TRUE)
nms = gsub("clin$age_at_initial_pathologic_diagnosis", "age", nms, fixed=TRUE)
nms
rownames(mm1) = 1:nrow(mm1)
colnames(mm1) = names(PCsDNA)[-1]
mm1 = data.frame(id=nms, mm1)
dim(mm1)
mm1[1:2,1:5]

xx = rbind(mm1, PCsDNA[1:6,])
yy = data.matrix(xx[,-1])
dim(yy)
ss = svd(yy)
ss$d

write.table(rbind(mm1, PCsDNA[1:6,]), file = "covariates_EA.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

