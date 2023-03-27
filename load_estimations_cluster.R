filename_sims       = paste0('./BVS.results.simulations.txt')
filename_params     = paste0('./BVS.results.parameters.txt')
output_sims         = read.csv(filename_sims,sep="\t")
output_params       = read.csv(filename_params,sep="\t")
output_simulation   = output_sims %>% dplyr::filter(type == 'simulation')
output_data         = output_sims %>% dplyr::filter(type == 'data')
# MCMC RESULTS - LOAD TO SAMPLE FROM THEM LATER
interactionMat_MCMC = readMat(paste0('./interactionMat_MCMC_',type,'.mat'))
growthRates_MCMC    = readMat(paste0('./growthRates_MCMC_',type,'.mat'))
burnin              = 2500 #burnin used for MCMC 

#### DATA AND SIMULATIONS TO DECIDE ON IC
# CHECK IC OF ALL TRAJECTORIES
allICs = c()
for (tj in unique(output_data$trajectory_ID)){
  output_data_original_traj    = output_data %>% filter(trajectory_ID==tj)
  timePoints                   = unique(output_data_original_traj$time)
  output_data_original_traj_IC = output_data_original_traj %>% filter(time==min(timePoints))
  output_data_original_traj_IC = as.numeric(output_data_original_traj_IC$abundance)
  allICs = rbind(allICs,output_data_original_traj_IC)
}
allICs_colmeans_keep = colMeans(allICs)
allICs_colsds_keep   = colSds(allICs)
allICs_colmeans      = allICs_colmeans_keep
allICs_colsds        = allICs_colsds_keep

output_growth   = output_params %>% filter(parameter_type=='growth_rate') 
output_interact = output_params %>% filter(parameter_type=='interaction') 
family_order    = output_growth$target_taxon #keep same order for interaction matrix
focalSpecies    = family_order
idxFocal        = match(focalSpecies,family_order)
colnamesVec     = focalSpecies



