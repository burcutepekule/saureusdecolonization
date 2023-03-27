plot_until_time = 0.25*(fuptime_tr/30)

if(plotOnTrt){
  source('plot_antibiotic_probiotic_cluster.R')
}

filnameSave=paste0(filname_base_data,'.RData')

if(saveData){
  save(simIdx,indSampleUse, dt, inoculumSpecs, numOfProTtreatmentDays, preTreatmentStart, perFrac, perFracVec, 
       proTreatmentStart, numOfBreakDays, numOfPreTtreatmentDays, 
       pm, pd , s2pVec, p2sVec, s2pMat, eps_in, fuptime_ss, inocSize_ss, fuptime_tr, 
       d_tr, popSize, popSizeLog, eps, s2p, p2s, numOfTotalCourses,inocSize, type,
       growth_rates, interaction_matrix, initial_conditions, SS, 
       out_all_antibiotic_steady,out_all_probiotic_steady,file = filnameSave)
  
  print(paste0('Data saved to ',filnameSave))
}

