#----------------------------------------------------#
# Concatenate Results Files                          #
#----------------------------------------------------#
setwd("results_Power_Calcs2/")

# Obtain file list and Order appropriately
file_listA = list.files(pattern = "pwr_")
file_listB = list.files(pattern = "foldS_1_psi0.1.RData")
file_list = intersect(file_listA,file_listB)
file_list

# Initialize Storage
pwr_mat = matrix(0,nrow = length(file_list),ncol = 7)

for(i in 1:length(file_list)){
  load(file_list[i])
  pwr_mat[i,] = pwr
}

save(pwr_mat,file = "../compiled_pwr_foldS_1.0_psi0.1_v2.RData")

# Obtain file list and Order appropriately
file_listA = list.files(pattern = "pwr_")
file_listB = list.files(pattern = "foldS_1_psi0.2.RData")
file_list = intersect(file_listA,file_listB)
file_list

# Initialize Storage
pwr_mat = matrix(0,nrow = length(file_list),ncol = 7)

for(i in 1:length(file_list)){
  load(file_list[i])
  pwr_mat[i,] = pwr
}

save(pwr_mat,file = "../compiled_pwr_folds_1.0_psi0.2_v2.RData")

#----------------------------------------------------#
# Concatenate Results Files                          #
#----------------------------------------------------#
setwd("../results_TypeI_Error2/")

# Obtain file list and Order appropriately
file_listA = list.files(pattern = "pwr_foldT_")
file_listB = list.files(pattern = "psi0.1.RData")
file_list = intersect(file_listA,file_listB)
file_list

# Initialize Storage
pwr_mat = matrix(0,nrow = length(file_list),ncol = 7)

for(i in 1:length(file_list)){
  load(file_list[i])
  pwr_mat[i,] = pwr
}

save(pwr_mat,file = "../compiled_pwr_foldT_1.0_psi0.1_v2.RData")

file_listA = list.files(pattern = "pwr_foldT_")
file_listB = list.files(pattern = "psi0.2.RData")
file_list = intersect(file_listA,file_listB)
file_list

# Initialize Storage
pwr_mat = matrix(0,nrow = length(file_list),ncol = 7)

for(i in 1:length(file_list)){
  load(file_list[i])
  pwr_mat[i,] = pwr
}

save(pwr_mat,file = "../compiled_pwr_foldT_1.0_psi0.2_v2.RData")

sessionInfo()

q(save="no")
