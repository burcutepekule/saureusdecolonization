# Setup
rm(list=ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(zoo)
library(Rcpp)
library(lubridate)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(prodlim)
library(matrixStats)
library(dplyr)
library(plyr)
library(readxl)
library(writexl)
# library(hydroGOF)
library(ggpubr)
library(ggplot2)
library(useful)
library(gridExtra)
# library(grid)
library(lattice)
library(gtable)
library(vegan)
library(deSolve)
library(tcpl)
library(Matrix)
library(R.matlab)
library(tryCatchLog)
library(rootSolve)
# library(logr)
type       = 'FORCED' #'ORIGINAL','FORCED'

# for (simIdx in 0:999){
for (simIdx in 0:9){
  
  ##################### ALL PARAMETERS 
  dt                     = 1/32; # 1 day (data is monthly (~32 days), so the interaction and growth parameters)
  inoculumSpecs          = c(1,2) # Commensal incoulum indexes for "Carnobacteriaceae" and "C1"
  numOfProTtreatmentDays = 5; #number of probiotic treatment days for one treatment course
  preTreatmentStart      = 15; #the day that the pre-treatment starts
  perFrac                = 0 #fraction of persistence
  perFracVec             = c(rep(0,14),perFrac,perFrac,perFrac,perFrac,rep(0,2)) #% persistent indexes : Sta1-4
  proTreatmentStart      = 0; #the day that the probiotic treatment starts (after preTreatment)
  numOfBreakDays         = 0; #number of break days between treatments  
  numOfPreTtreatmentDays = 5; #number of pre-treatment days for mupirocin treatment
  
  ########## PERSISTENCE PARAMETERS
  # ARTICLE : "Identifying determinants of persistent MRSA bacteremia using mathematical modeling" - TABLE 1 FOR SA
  # Death and growth
  pm  = 0.05; #growth rate modifier : less growth (Assumed that growth rate of persisters was 1/20 of that of normal SA.)
  pd  = 0; #death rate modifier : less death (they do not die at all in the paper above)
  # Switching rate
  s2pVec = 32*24*c(1e-5,1e-4,1e-3) # hourly to monthly conversion
  p2sVec = 32*24*c(1e-2,1e-1) # hourly to monthly conversion
  s2pMat = expand.grid(s2pVec,p2sVec)
  colnames(s2pMat)=c('s2p','p2s')
  s2pMat = rbind(c(0,0),s2pMat) #add 0-0 as control
  
  ########## PROBIOTIC TREATMENT PARAMETERS
  # inocSize is relative to the population size -> IF TOTAL LOAD IS 10^6, AND INOC IS 10^9, THEN inocSize = 10^3
  # if we assume that the biomass per sample (per swab) in our dataset is constant, then relative abundance = normalized abundance
  inocSizeVec   = 10^c(-1,0,1,2,3) #number of inoculum sizes to loop over
  inocSizeVec   = sort(inocSizeVec) 
  trtVec        = seq(0,8,1) #number of treatment courses to loop over
  
  # ## Steady State 
  fuptime_ss             = 12*30*10 #run long to get the steady state years
  d_tr                   = 15.5 #death rate by mupirocin treatment
  round_time             = 6
  plotOn_ss              = 0; #plot the steady state results
  solverMethod           = "impAdams"
  saveData               = TRUE #ONLY SAVE SUMMARY TO BE FASTER
  
  source('load_estimations_cluster.R') #Load the parameter estimation results
  source('load_models_cluster.R') #Load the ODE models
  options(warn=-1) #TO TURN BACK ON AGAIN, USE options(warn=0)
  
  ### FIND THE SYSTEMS THAT WORK FIRST
  # popSizeLogVec = seq(8,3,-1) #pop-size vector to loop over
  popSizeLog   = 6 #to make it faster
  popSize      = 10^popSizeLog
  eps          = 1/popSize #ODE sensitivity
  eps_in       = 100*popSize #ODE sensitivity
  print(simIdx)
  try.again     = TRUE
  while(try.again){
    
    # SAMPLE AN MCMC SYSTEM
    source('sample_system_cluster.R') # ONLY POP SIZE MATTERS, s2p=0, p2s=0 when finding a stable reasonable set
    dim(SS_all)
    is2useKeep = c()
    
    if(!an.error.occured){
      print(paste0('Sampled a new stable, reasonable, and Staph-Strep co-exsistent set. indSampleUse : ',indSampleUse))
      
      try.again = FALSE
      print(paste0('HERE : try.again? ',try.again))
    }
  }
  
  print(paste0('This simIdx : ',simIdx ,' uses indSampleUse ',indSampleUse))
  # SAMPLE AN MCMC SYSTEM
  filname_base_metadata = paste0('./SIM_META_SP_H_',type,'_',simIdx)
  filnameSave           = paste0(filname_base_metadata,'.RData')
  save(simIdx,indSampleUse,growth_rates,interaction_matrix,initial_conditions,SS_all,file = filnameSave)
}
