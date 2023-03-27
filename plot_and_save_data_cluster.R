plot_until_time = 0.25*(fuptime_tr/30)

if(plotOnTrt){
  source('plot_antibiotic_probiotic_cluster.R')
}

filnameSave=paste0(filname_base_data,'.RData')

if(saveData){
  save(simIdx,indSampleUse, dt, inoculumSpecs, recratio, numOfProTtreatmentDays, preTreatmentStart, resFrac, resFracVec, 
       proTreatmentStart, numOfBreakDays, numOfPreTtreatmentDays, 
       pm, pd , s2pVec, p2sVec, s2pMat, eps_in, fuptime_ss, d_ss, inocSize_ss, fuptime_tr, 
       d_tr, popSize, popSizeLog, eps, eps_time, s2p, p2s, numOfTotalCourses,inocSize, recratio, type,
       growth_rates, interaction_matrix, initial_conditions, SS, 
       sta_min_antibiotic,time2min_antibiotic,time2rec_antibiotic,
       out_antibiotic,out_all_antibiotic,out_all_antibiotic_final,shannon_antibiotic_final,R_antibiotic_species,
       sta_min_probiotic,time2min_probiotic,time2rec_probiotic,
       out_probiotic,out_all_probiotic,out_all_probiotic_final,shannon_probiotic_final,R_probiotic_species,
       plot_until_time,sta_ext_antibiotic,sta_ext_probiotic,
       inoculumTimes,preTreatmentTimes,treatmentTimes,file = filnameSave)
  
  print(paste0('Data saved to ',filnameSave))
}

