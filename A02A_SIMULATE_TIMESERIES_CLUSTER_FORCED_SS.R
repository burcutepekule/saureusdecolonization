# Setup
rm(list=ls())
dir = getwd()
setwd(dir)
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
# library(logr)
library(rootSolve)

# args <- commandArgs(TRUE)
# print(args)

# simIdx     = as.integer(args[1]) #0-999
# ind        = as.integer(args[2]) #1-7 #this needs to vary for different switching rates between susceptibles / persisters
simIdx     = 0
ind        = 2 #1-7 #pick one switching rate for local run
type       = 'FORCED' #'ORIGINAL','FORCED'

### REMEMBER THAT YOU NEED TO LOOP OVER ALL VARIABLES TO HAVE THE SUMMARY STATISTICS
### THIS WAS DONE IN A NON-LOCAL CLUSTER, THAT'S WHY HERE ONE FILE IS DEMONSTRATED

##################### ALL PARAMETERS HERE
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
fuptime_tr             = 12*30*10 #run long to get the steady state years
inocSize_ss            = 0 #inocize for steady state calc. -> 0
d_tr                   = 15.5 #death rate by mupirocin treatment
round_time             = 6
plotOn_ss              = 0; #plot the steady state results
solverMethod           = "impAdams"
saveData               = TRUE #ONLY SAVE SUMMARY TO BE FASTER

source('load_estimations_cluster.R')
source('load_models_cluster.R')

popSizeLog   = 6 
popSize      = 10^popSizeLog
eps          = 1/popSize
eps_in       = 100*popSize

filname_base_metadata = paste0('./SIM_META_SP_',type,'_',simIdx)
load(paste0(filname_base_metadata,'.RData'))
indSampleUseKeep = indSampleUse

# to make it faster

s2p             = s2pMat[ind,1]
p2s             = s2pMat[ind,2]
spcode          = -1*(10*log10(s2p/(32*24))+log10(p2s/(32*24)))
SS              = SS_all[ind,]

filname_base_SS = paste0('./SIM_SS_SP_H_',spcode,'_POPSIZE_',popSizeLog,'_',type,'_',simIdx,'_SS')
totalNumOfLoops = length(inocSizeVec)*length(trtVec)


out_all_keep_probiotic     = as.data.frame(matrix(NA, nrow = length(inocSizeVec)*length(trtVec), ncol = 6 + 40 + 1))
out_all_keep_antibiotic    = as.data.frame(matrix(NA, nrow = length(trtVec), ncol = 6 + 40 + 1))
counter_row_ant = 1
loopIdx         = 1;
print('lets go')

for (numOfTotalCourses in trtVec){
  
  an.error.occured = FALSE
  # RUN ANTIBIOTIC-TREATMENT
  inocSize=NA
  print('entered antibiotic')
  source('run_treatment_antibiotic_cluster_SS.R')
  out_all_antibiotic_steady = out_all_antibiotic
  print('antibiotic done')
  
  if(!an.error.occured){
    
    tempRow_out_all_keep_antibiotic = c(simIdx,indSampleUse,inocSize,numOfTotalCourses,s2p,p2s,
                                        out_all_antibiotic_steady,popSize)
    
    inocSizeCounter                 = 0
    counter_row_pro                 = 1
    out_all_keep_probiotic_temp     = as.data.frame(matrix(NA, nrow = length(inocSizeVec), ncol = 6 + 40 + 1))
    
    for(inocSize in inocSizeVec){
      
      print('entered probiotic')
      
      ########### PROBIOTIC-TREATMENT
      source('run_treatment_probiotic_cluster_SS.R')
      out_all_probiotic_steady = out_all_probiotic
      print('probiotic done')
      
      if(!an.error.occured){

        
        out_all_keep_probiotic_temp[counter_row_pro,]    = c(simIdx,indSampleUse,inocSize,numOfTotalCourses,s2p,p2s,
                                                             out_all_probiotic_steady,popSize)
        
        counter_row_pro = counter_row_pro + 1
        inocSizeCounter = inocSizeCounter + 1
        
        # PLOT AND SAVE #   THIS WAS FOR THE CASE WHERE I SAVED EACH SIM
        filname_base_data = paste0('./SIM_COURSES_',numOfTotalCourses,'_INOCSIZE_',log10(inocSize),
                                   '_SP_',spcode,'_POPSIZE_',popSizeLog,'_',type,'_H_',simIdx,'_SS')
        plotOnTrt = FALSE
        source('plot_and_save_data_cluster_SS.R')
        
        loopIdx = loopIdx + 1;
        
      }
      
    }
    
    # IF PROBIOTIC TRT WORKED FOR ALL INOC SIZES, SAVE ANTIBIOTIC
    if(inocSizeCounter==length(inocSizeVec)){
      out_all_keep_antibiotic[counter_row_ant,]    = tempRow_out_all_keep_antibiotic
      
      counter_row_ant = counter_row_ant + 1
      
      out_all_keep_probiotic    = rbind(out_all_keep_probiotic,out_all_keep_probiotic_temp)
    }
  }
}

filname_base_data_summary = paste0('./SIM_SUMMARY_SP_H_',spcode,'_POPSIZE_',popSizeLog,'_',type,'_',simIdx,'_SS.RData')
save(out_all_keep_antibiotic,out_all_keep_probiotic,file=filname_base_data_summary)
print('SS Summary saved.')

