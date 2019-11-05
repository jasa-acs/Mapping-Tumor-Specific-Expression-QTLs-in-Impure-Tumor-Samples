
library(multcomp)
library(MASS)
library(VGAM)
library(pTReCASE)


source("_asSeq/aseR.R")
source("_asSeq/aseR.TS.R")
source("_asSeq/trecaseR.R")
source("_asSeq/trecaseR.TS.R")
source("_asSeq/trecR.R")
source("_asSeq/trecR.TS.R")

#----------------------------------------------#
# Read in Command Line Arguments               #
#----------------------------------------------#
args<-commandArgs(TRUE)

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  foldT.val = 1.6
  foldS.val = 1.0
  psi_val   = 0.1
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#-------------------------------------------#
# Load the Functions                        #
#-------------------------------------------#

source("power.R")

#----------------------------------------------#
# Quick and Dirty Simulation                   #
#----------------------------------------------#

tmp.seed = as.integer(runif(1)*2e9)
set.seed(tmp.seed)
message("Used Seed: ",tmp.seed)

n      = 500
mu     = 100
rhoA   = 1.5
maf    = 0.2

phi    = 0.2
psi    = 0.2
theta  = 0.1
# use 10 simulation to demonstrate the code
# actual simulation results were based on 400 replicates
n.simu = 10
# n.simu = 400
alpha  = 0.05

pis    = 1-runif(n, 0, .5)
foldT  = foldT.val
foldS  = foldS.val

pwr = rep(0,7)

cat(date(), "\n")
pwr = powerNB(n, mu, foldT, foldS, phi, psi_val, n.simu, alpha, 
              maf, pis, rhoA)

file_name = sprintf("pwr_foldT_%s_foldS_%s_psi%s.RData", foldT.val, 
                    foldS.val, psi_val)

save(pwr, file = file_name)

sessionInfo()

q(save="no")


