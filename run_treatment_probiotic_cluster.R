########### PROBIOTICS-TREATMENT
fuptime                     = fuptime_tr #say 24mo.
d                           = d_tr; # here to actually kill
preTreatmentDeathRates      = rep(d,20) # INITIAL TREATMENT - REDUCE ICs 99%
preTreatmentDeathRates[2:3] = 0 # Corynebacteriaceae don't die
inoculumSizes_abs           = c(inocSize,inocSize) #if no inoculum, set to zero
source('run_model_save_plots_probiotic_per_cluster.R')
if(!an.error.occured){
  sta_all_df                 = as.data.frame(out_probiotic[,c(1,16)])
  sta_all_df$time_days       = sta_all_df$time/dt
  sta_all                    = out_probiotic[1:dim(out_probiotic)[1],c(16)]
  sta_all_after_treatment    = out_probiotic[(1+max(c(0,inoculumTimes))/dt):dim(out_probiotic)[1],c(16)]
  sta_all_after_treatment_df = sta_all_df[(1+max(c(0,inoculumTimes))/dt):dim(out_probiotic)[1],]
  
  sta_ext  = 0
  sta_min  = min(sta_all)
  sta_init = sta_all[1]
  time2min = out_probiotic[which(sta_all==sta_min),1]
  time2min = min(time2min)
  if(any(sta_all==0)){
    sta_ext  = 1
    time2min = min(time2min) #first extinction time point
    time2rec = Inf
  }else{
    sta_ext  = 0
    time2rec = sta_all_after_treatment_df[which(abs(sta_all_after_treatment-recratio*sta_init)==min(abs(sta_all_after_treatment-recratio*sta_init))),1]
    time2rec = time2rec-sta_all_after_treatment_df[1,]$time
    time2rec = min(time2rec)
    # time2rec = time2rec[which(time2rec>time2min)]#get the ones AFTER hitting the minimum
  }
  sta_ext_probiotic  = sta_ext
  sta_min_probiotic  = sta_min # for plotting
  time2min_probiotic = time2min # for plotting
  time2rec_probiotic = time2rec+sta_all_after_treatment_df[1,]$time # for plotting
  
  reccurance_limit_probiotic = recratio*sta_init # for plotting
  out_all_probiotic_final    = out_probiotic[dim(out_probiotic)[1],2:dim(out_probiotic)[2]]
  shannon_probiotic_final    = diversity(out_all_probiotic_final, index = "shannon", MARGIN = 1, base = exp(1))
  R_probiotic_species        = rowSums(out_all_probiotic_final != 0) # HOW TO FIND THE RICHNESS (NUM OF SPECIES)
}
