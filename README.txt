
------------------------------------------------------------------------------------
preparation to run the analysis
------------------------------------------------------------------------------------

-- Package dependencies
[1] pTReCASE_0.99.0     RcppEigen_0.3.3.4.0 Rcpp_0.12.18       
[4] VGAM_1.0-5          multcomp_1.4-10     TH.data_1.0-10     
[7] MASS_7.3-50         survival_2.42-6     mvtnorm_1.0-10     

		All the packages are available from CRAN, expect pTReCASE_0.99.0, which is available 
		at https://github.com/Sun-lab/pTReCASE
		
-- R version used to generate the results in the paper
		
	Since the initial data processing was conducted several years ago. The R versions used 
	in this paper range from R2.11.0 to R3.5.0. We have tested run earlier code in R3.5.0 
	and do not see any problem. The R functions that we used are very basic, such as read 
	or write files and thus we do not expect any R version dependence. 

-- Approximate run times
	It takes about 6 hours for one simulation configuration. The computation time can be 
	found in log files, e.g., results_Power_Calcs2/Linux_Sim_P_1_psi1_call.txt

-- Instal the pTReCASE package
  Since the R package contains C code, a C complier is required for installation. With 
  both R and appropriate c complier installed, this R package can be installed using the 
  following command (in Mac Terminal window or Windows command window) 

	R CMD install pTReCASE_0.99.0.tar.gz

------------------------------------------------------------------------------------
file structure
------------------------------------------------------------------------------------

This folder contains the codes for data preparation, simulation studies, real data 
analysis, as well as example data for real data analysis. 

Here is the folder structure in top two levels.

├── data_preparation
│   ├── R_batch1
│   ├── R_batch2
│   └── R_batch3
├── ex_real_data_analysis
│   ├── data
│   ├── pTReCASE_analysis.R
│   ├── pTReCASE_analysis.Rout
│   └── pTReCASE_ex_output.txt
└── simulation
    ├── Step1_simulations.R
    ├── Step1_simulations.Rout
    ├── Step1_simulations_job_submission.R
    ├── Step2_summarize_results.R
    ├── Step2_summarize_results.Rout
    ├── _asSeq
    ├── compiled_pwr_foldS_1.0_psi0.1_v2.RData
    ├── compiled_pwr_foldT_1.0_psi0.1_v2.RData
    ├── compiled_pwr_foldT_1.0_psi0.2_v2.RData
    ├── compiled_pwr_folds_1.0_psi0.2_v2.RData
    ├── power.R
    ├── results_Power_Calcs2
    └── results_TypeI_Error2

------------------------------------------------------------------------------------
data_preparation
------------------------------------------------------------------------------------

	data preparation pipelines to prepare gene expression and genotype data for real data 
	analysis. Most files have associated .Rout files that show the intermediate output.
	
	R_batch1: 
	
		step0-step4: codes used to generate genotype data, including calling genotype from
		affymetrix 6.0 array raw data (CEL files)

		step5-step6: PCA
		
		step7-step11: genotype phasing and imputation
		 
	R_batch2: 
	
		step1-step6: prepare total expression per gene per sample from bam files
		
		step7: genotype PCA for african ancestor (AA) or european ancestor (EA) separately. 
		the analysis done in this paper only use data from EA.
		
		step8-step11: prepare expression data and run an initial eQTL analysis without using 
		allele-specific expression (ASE).
		
	R_batch3:
	
		step1-step5: prepare allele-specific read counts. Some steps were run in parallele, 
		and the output were saved in folders named as _step2, _step3 etc.
		
		step6-step10: prepare data and run eQTL mapping using both total expression and ASE, 
		which is the old trecase method that do not account for tumor purity information. 
		

------------------------------------------------------------------------------------
ex_real_data_analysis
------------------------------------------------------------------------------------

	data:
	
		data used in real data analysis. The gene expression data are complete. The genotype 
		data cannot be shared due to data access restriction. Simulated genotype data of 1000 
		SNPs were provided. 
		
	pTReCASE_analysis.R/pTReCASE_analysis.Rout
		
		R code to run real data analysis, using the simulated genotype data instead of real 
		genotype data
		
	pTReCASE_ex_output.txt 
	
		example output

------------------------------------------------------------------------------------
simulation
------------------------------------------------------------------------------------

    ├── Step1_simulations.R
    ├── Step1_simulations.Rout
    ├── Step1_simulations_job_submission.R
    ├── Step2_summarize_results.R
    ├── Step2_summarize_results.Rout
    ├── _asSeq
    ├── compiled_pwr_foldS_1.0_psi0.1_v2.RData
    ├── compiled_pwr_foldT_1.0_psi0.1_v2.RData
    ├── compiled_pwr_foldT_1.0_psi0.2_v2.RData
    ├── compiled_pwr_folds_1.0_psi0.2_v2.RData
    ├── power.R
    ├── results_Power_Calcs2
    └── results_TypeI_Error2

	Step1_simulations.R/Step1_simulations_job_submission.R
	
		Step1_simulations.R has the code used to run simulation in one configuration. 
		Step1_simulations_job_submission.R runs it on clusters by submitting multiple jobs, 
		one for each simulation configuration. Step1_simulations_job_submission.R is 
		configured for a specific cluster and thus need to modified for other clusters. 
		
	Step2_summarize_results.R
		Summarize simulation results to generate the files that saved a power/typeI error 
		table: compiled_pwr_fold*
	
	compiled_pwr_fold*
		Simulation results generated by combining simulation results from individual 
		configurations. foldS_1.0 means it is simulation for power (foldS is fold change in 
		normal cells). foldT_1.0 means it is simulation for type I error (foldT is fold 
		change in tumor cells). 
	
	power.R
		Functions used to simulate data and run analysis using different methods. 
		
	_asSeq
		A folder including several R functions used for simulation. Some of the functions are 
		R version of C/C++ implementation in other packages. 
			trecaseR.R is an R version of the function trecase in R package asSeq 
			(https://github.com/Sun-lab/asSeq). trecaseR.TS.R is an R version of function 
			'pTReCASE_multapp' from R package pTReCASE. 
			
			aseR.R			eQTL mapping using ASE only. 
			aseR.TS.R   eQTL mapping using ASE only, and account for tumor/normal mixture
			trecR.R			eQTL mapping using total read count only.
			trecR.TS.R	eQTL mapping using total read count only, and account for tumor/normal mixture
			trecaseR.R	eQTL mapping using total read count and ASE.
			trecaseR.TS.R	eQTL mapping using total read count and ASE, and account for tumor/normal mixture


	results_Power_Calcs2
		Simulation results as well as log files for power analysis. 
		
	results_Power_Calcs2
		Simulation results as well as log files for type I error analysis. 


	NOTE:
	In order to reproduce the exact results, one must go into each .rout file in 
	/results_Power_Calcs2/ and /results_Power_Calcs2/ to extract the printed random seed.

	Example:
	Find the code which starts:

	> message("Used Seed: ",tmp.seed)
	Used Seed: 361409668

	Then use set.seed(tmp.seed)



