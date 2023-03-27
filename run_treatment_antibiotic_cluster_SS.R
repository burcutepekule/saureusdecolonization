########### ANTIBIOTIC-TREATMENT
fuptime                     = fuptime_tr #say 24mo.
d                           = d_tr; # here to actually kill
preTreatmentDeathRates      = rep(d,20) # INITIAL TREATMENT - REDUCE ICs 99%
preTreatmentDeathRates[2:3] = 0 # Corynebacteriaceae don't die
inoculumSizes_abs           = c(inocSize_ss,inocSize_ss) #if no inoculum, set to zero
source('run_model_save_plots_antibiotic_per_cluster_SS.R')

