
bsub -q bigmem -n 1 -M 20 /netscr/nzhao/hnsc_untrimmed/step4_MaCH/MaCH_sub_splitRef/submit_mach_chr1_10.sh 
bsub -q bigmem -n 1 -M 30 /netscr/nzhao/hnsc_untrimmed/step4_MaCH/sh/submit_mach_chr1_2.sh 

cd /netscr/nzhao/hnsc_untrimmed/step4_MaCH 
 /nas02/home/n/z/nzhao/softwares/mach.1.0.18/executables/mach1  -d chr20.6.dat  -p pedi.chr20  -s /netscr/nzhao/TCGA_hnsc/lib/1000Genome/snps/chr20.snps  -h /netscr/nzhao/TCGA_hnsc/lib/1000Genome/hap/all/20101123.chr20.hap.gz  --r 50 --states 258 --phase --autoflip --greedy --compact --interimInterval 10 -o /netscr/nzhao/hnsc_untrimmed/step4_MaCH/out/MaCH_chr20_6.out
