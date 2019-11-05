#-------------------------------------------#
# Running the Simulation Quickly            #
#-------------------------------------------#
setwd("/netscr/drwilson/2016-05-17 Paper 2/Programs/")

#-------------------------------------------#
# Initialize Sim: No normal eQTL            #
#-------------------------------------------#
folds = seq(1.0,1.7,by = 0.1)
for(i in 1:length(folds)){
  cmd_text = sprintf("bsub -q week -o ./Run_Output/Linux_Sim_P_%s_psi1_call.txt R CMD BATCH --no-save --no-restore '--args foldT.val=%s foldS.val=%s psi_val=0.1' Step1_simulations.R ./Run_Output/Linux_sim_P_%s_psi1.rout",
                     i,folds[i],1,i)
  system(cmd_text)
}

folds = seq(1.0,2.0,by = 0.2)
for(i in 1:length(folds)){
  cmd_text = sprintf("bsub -q week -o ./Run_Output/Linux_Sim_T1_%s_psi1_call.txt R CMD BATCH --no-save --no-restore '--args foldT.val=%s foldS.val=%s psi_val=0.1' Step1_simulations.R ./Run_Output/Linux_sim_T1_%s_psi1.rout",
                     i,1,folds[i],i)
  system(cmd_text)
}

folds = seq(1.0,1.7,by = 0.1)
for(i in 1:length(folds)){
  cmd_text = sprintf("bsub -q week -o ./Run_Output/Linux_Sim_P_%s_psi2_call.txt R CMD BATCH --no-save --no-restore '--args foldT.val=%s foldS.val=%s psi_val=0.2' Step1_simulations.R ./Run_Output/Linux_sim_P_%s_psi2.rout",
                     i,folds[i],1,i)
  system(cmd_text)
}

folds = seq(1.0,2.0,by = 0.2)
for(i in 1:length(folds)){
  cmd_text = sprintf("bsub -q week -o ./Run_Output/Linux_Sim_T1_%s_psi2_call.txt R CMD BATCH --no-save --no-restore '--args foldT.val=%s foldS.val=%s psi_val=0.2' Step1_simulations.R ./Run_Output/Linux_sim_T1_%s_psi2.rout",
                     i,1,folds[i],i)
  system(cmd_text)
}

