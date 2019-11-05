
library(MatrixEQTL)

setwd("/Volumes/Moon/TCGA_BRCA/data/")

grp='AA'

SNP_file_name        = sprintf("../genotype_normal/dat%s.txt", grp)
covariates_file_name = sprintf("covariates_%s.txt", grp)
expression_file_name = sprintf("log_TReC_%s.txt", grp)
output_file_name     = sprintf("eQTL_log_TReC_%s.txt", grp)

useModel              = modelLINEAR;
pvOutputThreshold     = 1e-5;

#errorCovariance = numeric();

snps = SlicedData$new(); 
snps$fileDelimiter = '\t'; # the TAB character 
snps$fileOmitCharacters = 'NA'; # denote missing values; 
snps$fileSkipRows = 1; # one row of column labels 
snps$fileSkipColumns = 1; # one column of row labels 
snps$fileSliceSize = 2000; # read file in pieces of 2000 rows 
snps$LoadFile( SNP_file_name );

genes = SlicedData$new(); 
genes$fileDelimiter = '\t'; # the TAB character 
genes$fileOmitCharacters = 'NA'; # denote missing values; 
genes$fileSkipRows = 1; # one row of column labels 
genes$fileSkipColumns = 1; # one column of row labels 
genes$fileSliceSize = 2000; # read file in pieces of 2000 rows 
genes$LoadFile( expression_file_name );


cvrt = SlicedData$new(); 
cvrt$fileDelimiter = '\t'; # the TAB character 
cvrt$fileOmitCharacters = 'NA'; # denote missing values; 
cvrt$fileSkipRows = 1; # one row of column labels 
cvrt$fileSkipColumns = 1; # one column of row labels 
cvrt$fileSliceSize = 2000; # read file in pieces of 2000 rows 
cvrt$LoadFile( covariates_file_name );

gc()

me = Matrix_eQTL_main(
snps, 
genes, 
cvrt, 
output_file_name = output_file_name, 
pvOutputThreshold = pvOutputThreshold,
useModel = modelLINEAR, 
errorCovariance = numeric(), 
verbose = TRUE, 
pvalue.hist = T
)


