apt-probeset-genotype \
  -o ../genotype_normal \
  -c /Volumes/Moon/TCGA/lib/GenomeWideSNP_6.cdf \
  --set-gender-method cn-probe-chrXY-ratio \
  --chrX-probes /Volumes/Moon/TCGA/lib/GenomeWideSNP_6.chrXprobes \
  --chrY-probes /Volumes/Moon/TCGA/lib/GenomeWideSNP_6.chrYprobes \
  --special-snps /Volumes/Moon/TCGA/lib/GenomeWideSNP_6.specialSNPs \
  --read-models-birdseed /Volumes/Moon/TCGA/lib/GenomeWideSNP_6.birdseed-v2.models \
  -a birdseed-v2 \
  --cel-files ../info/cel_files_normal_after_qc.txt
