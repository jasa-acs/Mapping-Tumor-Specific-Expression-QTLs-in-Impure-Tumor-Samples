
library(asSeq, lib="/nas02/home/w/e/weisun/R/Rlibs/")

# -----------------------------------------------------------------
# read in sample information
# -----------------------------------------------------------------

setwd("/lustre/scr/w/e/weisun/TCGA/info/")

sam = read.table("brca_samples2use_after_qc.txt", 
header=TRUE, sep="\t", as.is=TRUE)
dim(sam)
sam[1,]

bams = paste(sam$RNA_rnaseqID.Tumor, "sorted_genome_alignments.bam", sep="")
bams[1:3]

# -----------------------------------------------------------------
# check the bam files
# -----------------------------------------------------------------

setwd("/lustre/scr/w/e/weisun/TCGA/bam/")
ffs = list.files(pattern="_alignments.bam")

length(ffs)
ffs[1:5]

bam2use = intersect(ffs, bams)
length(bam2use)

all(sam$DNAnorml_patientID == sam$RNA_patientID)
samples = sam$RNA_patientID[match(bam2use, bams)]

  
sami   = samples[i]
  
# ----------------------------------------------------------
# counting
# ----------------------------------------------------------
ctF  = sprintf("_count/count_%s.txt", sami)
cmd1 = sprintf("samtools view %s | wc -l >> %s\n", bam2use[i], ctF)
system(cmd1)

# ----------------------------------------------------------
# sorting
# ----------------------------------------------------------
cmd2 = sprintf("samtools sort -n %s %s_sorted_by_name", bam2use[i], sami)
system(cmd2)
bamF = sprintf("%s_sorted_by_name.bam", sami)

# ----------------------------------------------------------
# getUnique and filtering
# ----------------------------------------------------------
prepareBAM(bamF, sprintf("%s_sorted_by_name", sami), sortIt=FALSE)

system(sprintf("rm %s", bamF))

# ----------------------------------------------------------
# counting again
# ----------------------------------------------------------
cmd3   = sprintf("samtools view %s_sorted_by_name_uniq_filtered.bam | wc -l >> %s\n", sami, ctF)
system(cmd3)



