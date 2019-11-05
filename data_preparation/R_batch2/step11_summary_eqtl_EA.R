
setwd("/Volumes/Moon/TCGA_BRCA/data/")

# ------------------------------------------------------------------------
# read in eQTL results
# ------------------------------------------------------------------------

pvs = read.table("eQTL_log_TReC_EA.txt", 
  header=TRUE, sep="\t", as.is=TRUE)
dim(pvs)
pvs[1:2,]

summary(pvs$p.value)
summary(pvs$FDR)
summary(pvs$p.value[pvs$FDR < 0.05])
summary(pvs$FDR[pvs$p.value < 1e-9])

table(pvs$FDR < 0.2)
table(pvs$p.value < 1e-12)
table(pvs$p.value < 1e-9)
table(pvs$p.value < 1e-8)

# ------------------------------------------------------------
# read in SNP annotation 
# ------------------------------------------------------------

setwd("/Users/suninsky/research/data/Affy/Affy6/anno/")

if(file.exists("affy6_snp_loci.txt")){
  infoD = read.table("affy6_snp_loci.txt", sep="\t", header=TRUE, 
                    as.is=TRUE, na.string="---")
}else{
  info = read.csv("GenomeWideSNP_6.na32.annot.csv", as.is=TRUE, 
                  comment.char="#")
  dim(info)
  info[1:2,1:4]

  infoD = info[,1:4]
  write.table(infoD, file = "affy6_snp_loci.txt", append = FALSE, 
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

names(infoD)
table(infoD$Chromosome)
summary(infoD$Physical.Position)

dim(infoD)
infoD[1:2,]

# ------------------------------------------------------------
# read in gene annotation 
# ------------------------------------------------------------

setwd("/Users/suninsky/research/data/human/")

infoE = read.table("ens_gene_loci.txt", as.is=TRUE, sep="\t", header=TRUE)
dim(infoE)
infoE[1:2,]
table(infoE$chr)

# ------------------------------------------------------------------------
# read in location information and plot it
# ------------------------------------------------------------------------

setwd("/Volumes/Moon/TCGA_BRCA/")

ptplot   = pvs[which(pvs$gene %in% infoE$id),]
geneID   = match(ptplot$gene, infoE$id)
markerID = match(ptplot$SNP,  infoD$Probe.Set.ID)
scores   = ptplot$p.value
scuts    = c(1e-9, 1e-14, 1e-19, 1e-24)
cols     = c("green", "blue", "red", "black")
eChr     = gsub("chr", "", infoE$chr)
ePos     = 0.5*(infoE$start + infoE$end)
mChr     = gsub("chr", "", infoD$Chromosome)
mPos     = infoD$Physical.Position
chroms   = 1:22


source("codes_Wonil/eQTL_Plot.R")
xlab="eQTL Location"; ylab="Transcript Location";
plot.hotspots=TRUE; hotspots.cut=10; score.type="p-value"

png("figures/eQTL_logTReC_EA.png", width=7.5, height=9, res=200, units="in")
eqtl.plot(geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
mPos, chroms, xlab="eQTL Location", ylab="Transcript Location",
plot.hotspots=TRUE, hotspots.cut=10, score.type="p-value")
dev.off()

# ------------------------------------------------------------------------
# check hotspots
# ------------------------------------------------------------------------

tb = table(pvs$SNP[pvs$p.value<1e-9])
table(tb)

tb1 = tb[tb>10]
tb1

mat1 = match(names(tb1), infoD$Probe.Set.ID)

info1 = data.frame(infoD[mat1,], nEQTL=as.numeric(tb1))
table(info1$Chromosome)
info1$Chromosome = as.numeric(info1$Chromosome)
info1 = info1[order(info1$Chromosome, info1$Physical.Position),]
info1

genes = pvs$gene[pvs$SNP=="SNP_A-2223505" & pvs$p.value < 1e-9]
cat(genes, sep="\n")

genes = pvs$gene[pvs$SNP=="SNP_A-2013353" & pvs$p.value < 1e-9]
cat(genes, sep="\n")



