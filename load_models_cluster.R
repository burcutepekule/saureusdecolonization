################################################ ODE MODELS ################################################ 
nasal.model.challenge <- function (t, y, params) { #same as probiotic
  
  eps_in    = params$eps_in
  dt_in     = params$dt_in
  muVec     = params$mu
  alpha     = params$alpha
  challengeTime   = params$challengeTime
  inoculumSizeVec = params$inoculumSizes
  inoculumSpecVec = params$inoculumSpecs
  s2p             = params$s2p #rate of transition from sensitive to persister
  p2s             = params$p2s #rate of transition persister to sensitive
  pm              = params$pm #modifier : decay in growth rate
  pd              = params$pd #modifier : decay in death rate
  # s2p << p2s (difference is around the order of 5)
  # MODIFY FOR PERSISTENT VARIANTS
  numSpecies      = 0.5*length(muVec) #keep the original number of actual species then divide to persister/sensitive
  numCompartments = 2*numSpecies #1:numSpecies -> sensitive, numSpecies+1:2*numSpecies -> persister
  
  #########################################################################################
  ##### FOR NOW, PERSISTERS APPLY ONLY FOR SA
  
  muVec[numSpecies+c(15)]      = pm*muVec[numSpecies+c(15)]# persisters grow slower! (only for Sta1)
  muVec[numSpecies+c(1:14,16:20)] = 0
  s2pVec = rep(0,numCompartments)
  s2pVec[c(15)] = s2p #SENSITIVE SA TO PERSISTERS
  p2sVec = rep(0,numCompartments)
  p2sVec[numSpecies+c(15)] = p2s #PERSISTERS SA TO SENSITIVES
  #########################################################################################
  
  dydtVec = rep(0,numCompartments) #add an additional compartment for each resistant variant
  
  for(i in 1:numCompartments){
    
    sumInteractions = 0
    for(j in 1:numCompartments){
      sumInteractions = sumInteractions+y[i]*y[j]*alpha[i,j]
    }
    
    if(i<=numSpecies){ 
      dydtVec[i] = +muVec[i]*y[i] - s2pVec[i]*y[i] + p2sVec[i+numSpecies]*y[i+numSpecies]+ sumInteractions
    }else{#in persister population
      # persisters don't die but also do not grow ***when there is treatment***
      dydtVec[i] = +muVec[i]*y[i] + s2pVec[i-numSpecies]*y[i-numSpecies] - p2sVec[i]*y[i]+ sumInteractions
    }
  }
  dydt <- c(dydtVec)
  list(dydt)
}



nasal.model.probiotic <- function (t, y, params) {
  
  eps_in    = params$eps_in
  dt_in     = params$dt_in
  muVec     = params$mu
  alpha     = params$alpha
  inoculumTimeVec = params$inoculumTimes
  inoculumSizeVec = params$inoculumSizes
  inoculumSpecVec = params$inoculumSpecs
  s2p             = params$s2p #rate of transition from sensitive to persister
  p2s             = params$p2s #rate of transition persister to sensitive
  pm              = params$pm #modifier : decay in growth rate
  pd              = params$pd #modifier : decay in death rate
  # s2p << p2s (difference is around the order of 5)
  # MODIFY FOR PERSISTENT VARIANTS
  numSpecies      = 0.5*length(muVec) #keep the original number of actual species then divide to persister/sensitive
  numCompartments = 2*numSpecies #1:numSpecies -> sensitive, numSpecies+1:2*numSpecies -> persister
  
  #########################################################################################
  ##### FOR NOW, PERSISTERS APPLY ONLY FOR SA
  
  muVec[numSpecies+c(15)]  = pm*muVec[numSpecies+c(15)]# persisters grow slower! (only for Sta1)
  muVec[numSpecies+c(1:14,16:20)] = 0
  s2pVec = rep(0,numCompartments)
  s2pVec[c(15)] = s2p #SENSITIVE SA TO PERSISTERS
  p2sVec = rep(0,numCompartments)
  p2sVec[numSpecies+c(15)] = p2s #PERSISTERS SA TO SENSITIVES
  #########################################################################################

  dydtVec = rep(0,numCompartments) #add an additional compartment for each resistant variant

  for(i in 1:numCompartments){
    
    sumInteractions = 0
    for(j in 1:numCompartments){
      sumInteractions = sumInteractions+y[i]*y[j]*alpha[i,j]
    }
    
    if(i<=numSpecies){ 
      dydtVec[i] = +muVec[i]*y[i] - s2pVec[i]*y[i] + p2sVec[i+numSpecies]*y[i+numSpecies]+ sumInteractions
    }else{#in persister population
      # persisters don't die but also do not grow ***when there is treatment***
      dydtVec[i] = +muVec[i]*y[i] + s2pVec[i-numSpecies]*y[i-numSpecies] - p2sVec[i]*y[i]+ sumInteractions
    }
  }
  dydt <- c(dydtVec)
  list(dydt)
}


nasal.model.antibiotic <- function (t, y, params) {
  
  eps_in    = params$eps_in
  dt_in     = params$dt_in
  muVec     = params$mu
  alpha     = params$alpha
  treatmentTimeVec     = params$treatmentTimes
  preTreatmentDeathVec = params$preTreatmentDeathRates
  s2p             = params$s2p #rate of transition from sensitive to persister
  p2s             = params$p2s #rate of transition persister to sensitive
  pm              = params$pm #modifier : decay in growth rate
  pd              = params$pd #modifier : decay in death rate
  # s2p << p2s (difference is around the order of 5)
  # MODIFY FOR PERSISTENT VARIANTS
  numSpecies      = 0.5*length(muVec) #keep the original number of actual species then divide to persister/sensitive
  numCompartments = 2*numSpecies #1:numSpecies -> sensitive, numSpecies+1:2*numSpecies -> persister
  #########################################################################################
  ##### FOR NOW, PERSISTERS ONLY FOR SA

  muVec[numSpecies+c(15)]      = pm*muVec[numSpecies+c(15)]# persisters grow slower! (only for Sta1)
  muVec[numSpecies+c(1:14,16:20)] = 0
  preTreatmentDeathVec = c(preTreatmentDeathVec,preTreatmentDeathVec)
  preTreatmentDeathVec[numSpecies+c(15)] = pd*preTreatmentDeathVec[numSpecies+c(15)]   # persisters die slower! (only for SA)
  # ONLY FOR SA FOR NOW
  s2pVec = rep(0,numCompartments)
  s2pVec[c(15)] = s2p #SENSITIVE SA TO PERSISTERS
  p2sVec = rep(0,numCompartments)
  p2sVec[numSpecies+c(15)] = p2s #PERSISTERS SA TO SENSITIVES

  #########################################################################################
  
  dydtVec   = rep(0,numCompartments) #add an additional compartment for each resistant variant

  # for trt variables
  dydtVec         = c(dydtVec,0,0)
  antiTreatmentOn = as.numeric((y[numCompartments+1]+y[numCompartments+2])>0) #either pre or post trt

  for(i in 1:numCompartments){
    
    sumInteractions = 0
    for(j in 1:numCompartments){
      sumInteractions = sumInteractions+y[i]*y[j]*alpha[i,j]
    }
    
    if(i<=numSpecies){ #in sensitive population
      dydtVec[i] = +muVec[i]*y[i] - s2pVec[i]*y[i] + p2sVec[i+numSpecies]*y[i+numSpecies]+ sumInteractions- antiTreatmentOn*preTreatmentDeathVec[i]*y[i]
    }else{#in persister population
      # persisters don't die but also do not grow ***when there is treatment***
      dydtVec[i] = +muVec[i]*y[i] + s2pVec[i-numSpecies]*y[i-numSpecies] - p2sVec[i]*y[i]+ sumInteractions - antiTreatmentOn*preTreatmentDeathVec[i]*y[i]
    }
  }
  dydt <- c(dydtVec)
  list(dydt)
}


