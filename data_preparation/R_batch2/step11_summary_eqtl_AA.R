
setwd("/Volumes/Moon/TCGA_BRCA/data/")

# ------------------------------------------------------------------------
# read in eQTL results
# ------------------------------------------------------------------------

pvs = read.table("eQTL_log_TReC_AA.txt", 
  header=TRUE, sep="\t", as.is=TRUE)
dim(pvs)
pvs[1:2,]

summary(pvs$p.value)
summary(pvs$FDR)
table(pvs$FDR < 0.2)
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
table(infoD$Strand)

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
scuts    = c(1e-8, 1e-9, 1e-10, 1e-12)
cols     = c("green", "blue", "red", "black")
eChr     = gsub("chr", "", infoE$chr)
ePos     = 0.5*(infoE$start + infoE$end)
mChr     = gsub("chr", "", infoD$Chromosome)
mPos     = infoD$Physical.Position;
chroms   = 1:22


source("codes_Wonil/eQTL_Plot.R")

png("figures/eQTL_logTReC_AA.png", width=7.5, height=9, res=200, units="in")
eqtl.plot(geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
mPos, chroms, xlab="eQTL Location", ylab="Transcript Location",
plot.hotspots=TRUE, hotspots.cut=10, score.type="p-value")
dev.off()

tb = table(pvs$SNP[pvs$p.value<1e-8])
table(tb)
