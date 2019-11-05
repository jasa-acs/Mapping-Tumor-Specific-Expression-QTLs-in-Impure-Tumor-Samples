
##
## take 0 gene expression PCs and 2 genotype PC
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

PCs = read.table("PCs_log_TReC_AA.txt", sep="\t", header=TRUE)
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
# check genotype PCs
# --------------------------------------------------------

PCsDNA = read.table("PCs_genotype_AA.txt", sep="\t", header=TRUE)
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

pvals = matrix(NA, nrow=20, ncol=4)

for(i in 1:20){
  PCi = as.numeric(PCsDNA[i,-1])
  ai  = anova(lm(PCi ~  as.factor(sam1$DNAnorml_plate) 
                 + as.factor(sam1$DNAnorml_institution)
                 + as.factor(sam1$DNAnorml_type)
                 + as.factor(sam1$DNAnorml_portion)))
  pvals[i,] = ai$Pr[1:4]
  
}

signif(pvals,2)

# --------------------------------------------------------
# write out results
# --------------------------------------------------------

write.table(PCsDNA[1:2,], file = "covariates_AA.txt", append = FALSE, 
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
