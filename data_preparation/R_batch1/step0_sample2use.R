
# --------------------------------------------------------------------
# check the number of samples to use, based on the QC file from Katie
# brca.txt
# --------------------------------------------------------------------

setwd("~/research/TCGA/BRCA/info")

# --------------------------------------------------------------------
# read in datafreeze information
# --------------------------------------------------------------------

freeze1 = read.table("BRCA.DataFreeze.20120907.txt", 
  sep="\t", header=TRUE, as.is=TRUE)
dim(freeze1)
freeze1[1:2,]
length(which(freeze1$RNAseq.Normal != ""))

freeze2 = read.table("BRCA.819.DataFreeze.txt", 
  sep="\t", header=FALSE, as.is=TRUE)
names(freeze2) = c("sample", "rnaseqID")
dim(freeze2)
freeze2[1:2,]

rnaseq = c(freeze1$RNAseq.Tumor, freeze1$RNAseq.Normal)
rnaseq = rnaseq[which(rnaseq != "")]
all(rnaseq %in% freeze2$sample)

r1 = rnaseq[which(!(rnaseq %in% freeze2$sample))]
r1
r1 %in% freeze1$RNAseq.Normal

mat1 = match(freeze1$RNAseq.Tumor, freeze2$sample)
table(is.na(mat1))
if(!all(freeze1$RNAseq.Tumor == freeze2$sample[mat1])){stop("wrong!\n")}
r1 = freeze2$rnaseqID[mat1]

mat2 = match(freeze1$RNAseq.Normal, freeze2$sample)
table(is.na(mat2))
r2 = freeze2$rnaseqID[mat2]

info = data.frame(freeze1, rnaseqID.Tumor=r1, rnaseqID.Normal=r2, 
stringsAsFactors=FALSE)
dim(info)
head(info)

# --------------------------------------------------------------------
# check RNA sample information
# --------------------------------------------------------------------

sams = strsplit(info$RNAseq.Tumor, split="-")

table(sapply(sams, length))

sams = matrix(unlist(sams), byrow=TRUE, ncol=7)
dim(sams)
sams[1:2,]

table(sams[,1])
table(sams[,7])

sams = sams[,-c(1,7)]
colnames(sams) = c("institution", "patientID", "type", "portion", "plate")
dim(sams)
sams[1:2,]

apply(sams[,c(1,3:5)], 2, table)
table(table(sams[,2]))

info = data.frame(info, sams, stringsAsFactors=FALSE)
dim(info)
head(info)

# --------------------------------------------------------------------
# read in DNA information
# --------------------------------------------------------------------

DNAinfo = read.table("broad.mit.edu_BRCA.Genome_Wide_SNP_6.sdrf.txt", 
header=TRUE, sep="\t", as.is=TRUE)
dim(DNAinfo)
DNAinfo[1:2,1:9]
names(DNAinfo)

DNAsams = strsplit(DNAinfo$Comment..TCGA.Barcode., split="-")

table(sapply(DNAsams, length))

DNAsams = matrix(unlist(DNAsams), byrow=TRUE, ncol=7)
dim(DNAsams)
DNAsams[1:2,]

table(DNAsams[,1])
table(DNAsams[,7])

DNAsams = DNAsams[,-c(1,7)]
colnames(DNAsams) = c("institution", "patientID", "type", "portion", "plate")

apply(DNAsams[,c(1,3:5)], 2, table)
table(table(DNAsams[,2]))

# --------------------------------------------------------------------
# check the samples with only one DNA sample
# --------------------------------------------------------------------

DNAsams = data.frame(DNAsams, stringsAsFactors=FALSE)
t1 = table(DNAsams$patientID)

D1 = DNAsams[DNAsams$patientID %in% names(t1)[t1==1],]

table(D1$type)

table(D1$patientID %in% info$patientID)

# --------------------------------------------------------------------
# check the samples with two DNA sample
# --------------------------------------------------------------------

fun1 <- function(v){paste(sort(v), collapse="_")}

D2 = DNAsams[DNAsams$patientID %in% names(t1)[t1==2],]
types = tapply(D2$type, D2$patientID, fun1)
table(types)

# --------------------------------------------------------------------
# check the samples with three or four DNA sample
# --------------------------------------------------------------------

D3 = DNAsams[DNAsams$patientID %in% names(t1)[t1>=3],]
types = tapply(D3$type, D3$patientID, fun1)
table(types)

## every sample has type 01A
all(grepl("01A", types))

## remove 06A or 01B
D3 = D3[which(D3$type != "06A" & D3$type != "01B"),]
table(D3$type)
types = tapply(D3$type, D3$patientID, fun1)
table(types)

samsD31 = names(types)[grep("10A", types)]
samsD32 = names(types)[-grep("10A", types)]
samsD31
samsD32

ww1 = which(D3$patientID %in% samsD31 & D3$type %in% c("01A", "10A"))
ww2 = which(D3$patientID %in% samsD32)

D3 = D3[c(ww1, ww2),]
table(table(as.character(D3$patientID)))
types = tapply(D3$type, D3$patientID, fun1)
table(types)

# --------------------------------------------------------------------
# DNA 2use
# --------------------------------------------------------------------

D2kp = rbind(D2, D3)
dim(D2kp)
head(D2kp)
table(table(as.character(D2kp$patientID)))

samName = paste("TCGA", apply(D2kp, 1, paste, collapse="-"), "01", sep="-")
mat1    = match(samName, DNAinfo$Comment..TCGA.Barcode.)
all(samName == DNAinfo$Comment..TCGA.Barcode.[mat1])

D2kp[["arrayFile"]] = paste(DNAinfo$Hybridization.Name[mat1], "CEL", sep=".")

# --------------------------------------------------------------------
# read in DNA sample list
# --------------------------------------------------------------------

dna = scan("brca_cel_list.txt", what=character(0))
length(dna)
dna[1:2]

all(D2kp$arrayFile %in% D2kp$arrayFile)

# --------------------------------------------------------------------
# intersections
# --------------------------------------------------------------------

dim(D2kp)
head(D2kp)
fun1 <- function(v){paste(sort(v), collapse="_")}
types = tapply(D2kp$type, D2kp$patientID, fun1)
table(types)

table(as.character(D2kp$type))

sam2use = intersect(info$patientID, D2kp$patientID)
length(sam2use)

wtumor = which(D2kp$type %in% c("01A", "01B"))
Dtumor = D2kp[wtumor,]
Dnorml = D2kp[-wtumor,]

mat1   = match(sam2use, info$patientID)
mat2   = match(sam2use, Dtumor$patientID)
mat3   = match(sam2use, Dnorml$patientID)

all(info$patientID[mat1] == sam2use)
all(Dtumor$patientID[mat2] == sam2use)
all(Dnorml$patientID[mat3] == sam2use)

info   = info[mat1,]
Dtumor = Dtumor[mat2,]
Dnorml = Dnorml[mat3,]

RNAsam = paste("TCGA", apply(info, 1, paste, collapse="-"), "07", sep="-")
all(info$sample == RNAsam)

names(Dtumor) = paste("DNAtumor", names(Dtumor), sep="_")
names(Dnorml) = paste("DNAnorml", names(Dnorml), sep="_")
names(info)   = paste("RNA", names(info), sep="_")

db = data.frame(Dtumor, Dnorml, info, stringsAsFactors=FALSE)
dim(db)
db[1:2,]

if(! all(db$DNAtumor_patientID == info$RNA_patientID)){
  stop("sample ID do not match\n")
}

if(! all(db$DNAnormal_patientID == info$RNA_patientID)){
  stop("sample ID do not match\n")
}

# --------------------------------------------------------------------
# write out results
# --------------------------------------------------------------------

write.table(db, file = "brca_samples2use.txt", append = FALSE, quote = FALSE, 
  sep = "\t", row.names = FALSE, col.names = TRUE)

# --------------------------------------------------------------------
# prepare CEL file name list
# --------------------------------------------------------------------

fname = "cel_files_normal.txt"
cat("cel_files\n", file=fname)
cels = paste("/Volumes/Moon/TCGA_BRCA/data_DNA/", db$DNAnorml_arrayFile, sep="")
cat(cels, file=fname, sep="\n", append=TRUE)


fname = "cel_files_tumor.txt"
cat("cel_files\n", file=fname)
cels = paste("/Volumes/Moon/TCGA_BRCA/data_DNA/", db$DNAtumor_arrayFile, sep="")
cat(cels, file=fname, sep="\n", append=TRUE)


