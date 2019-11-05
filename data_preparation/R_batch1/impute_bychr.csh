#!bin/csh
#####################################
## step2: imputation ##
## Input file: pedigree file and dat file for Mach 
## input parameter $1: which chromsome to run $2: which region in the previously defined chromsome to run 
## submitted as: bsub -q week -M 20 csh impute_bychr2.csh 1 1
#######################
# note that one does NOT need to split target once reference is splitted because
#   markers in target but not in reference will be automatically ignored 
foreach chr ($1)
  @ numSegs = `wc -l < ./splitRefVCF/chr$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz.splitlog.csv` - 1
  echo "number of segments for chr$chr is $numSegs"
  foreach i ($2)
    @ corestart = `sed 1,1d ./splitRefVCF/chr$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz.splitlog.csv | head -$i | tail -1 | cut -f2 -d ',' `
    @ coreend   = `sed 1,1d  ./splitRefVCF/chr$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz.splitlog.csv | head -$i | tail -1 | cut -f3 -d ',' `
   /nas02/home/n/z/nzhao/softwares/linux/mach-admix --forceImputation -d /netscr/nzhao/hnsc_untrimmed/step4_MaCH/pedi_marker/marker.chr$chr -p /netscr/nzhao/hnsc_untrimmed/step4_MaCH/pedi_marker/pedi.chr$chr -h ./splitRefVCF/chr$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.region$i.gz --vcfRef --outputstart $corestart --outputend $coreend --phase --rounds 50 --states 268 --autoflip  -o my.chr$chr.region$i.out > my.chr$chr.region$i.log
  end
end

