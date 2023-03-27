# FIND A SYSTEM
is2use = FALSE
tried  = 1
while(!is2use){
  indSampleUse = sample(seq(1,22500), 1, replace = FALSE, prob = NULL) #sample index after burnin
  print( paste0('System type : ',type,', log10(popSize): ',popSizeLog,', simulation index : ',simIdx,', trying MCMC set # ',tried))
  source('sample_posteriors_cluster.R') #sample interaction terms from the 
  # OUTPUTS : growth_rates, interaction_matrix, and initial_conditions
  if(is2use){
  an.error.occured = FALSE
  source('run2ss_cluster.R')
  # OUTPUTS : out_all_final, out_final, SS, isStable, isReasonable, isSta1Stre1, is2use
  }
  tried=tried+1
}


