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
args <- commandArgs(TRUE)
print(args)

# NEED TO LOOP OVER THESE IN YOUR OWN NON-LOCAL CLUSTER. LOCAL COMPUTATIONS WOULD COST A LOT OF TIME
simIdx     = as.integer(args[1]) #0-999
ind        = as.integer(args[2]) #1-7
type       = 'FORCED' #'ORIGINAL','FORCED'

##################### ALL PARAMETERS HERE
dt                     = 1/32; # 1 day (data is monthly, therefore the interaction and growth parameters -> 1/30 leads to numerical issues!)
inoculumSpecs          = c(1,2) # "Carnobacteriaceae" and "C1"
numOfProTtreatmentDays = 5; #days
preTreatmentStart      = 15; # days
resFrac                = 0
resFracVec             = c(rep(0,14),resFrac,resFrac,resFrac,resFrac,rep(0,2)) #% persistent Sta1-4
proTreatmentStart      = 0; #days (after preTreatment)
numOfBreakDays         = 0; #days
numOfPreTtreatmentDays = 5; #days

########## PERSISTENCE PARAMETERS
pm     = 0.05; #less growth (Assumed that growth rate of persister was 1/20 of that of normal SA.)
pd     = 0; #less death (they do not die at all in the paper above)
s2pVec = 32*24*c(1e-5,1e-4,1e-3) # hourly to monthly (BUT conversion rate shouldn't be more than 1?)
p2sVec = 32*24*c(1e-2,1e-1) # hourly to monthly
s2pMat = expand.grid(s2pVec,p2sVec)
colnames(s2pMat)=c('s2p','p2s')
s2pMat = rbind(c(0,0),s2pMat) #add 0-0 as control
# s2pVec = 0 # hourly to monthly (BUT conversion rate shouldn't be more than 1?)
# p2sVec = 0 # hourly to monthly

########## PROBIOTIC TREATMENT PARAMETERS
# IF TOTAL LOAD IS 10^6, AND INOC IS 10^9, THEN inocSize = 10^3
# if we assume that the biomass per sample (per swab) in our dataset is constant, then relative abundance = normalized abundance
# inocSizeVec   = 10^c(-4,-3,-2,-1,0,1)
inocSizeVec   = 10^c(-1,0,1,2,3)
inocSizeVec   = sort(inocSizeVec)
trtVec        = seq(0,8,1)

## SS

fuptime_ss             = 12*30*100 #100 ?! years
d_ss                   = 0
inocSize_ss            = 0
fuptime_tr             = 30*12*3 #2Y is good (MAKE IT 3Y TO BE ON THE SAFE SIDE)
d_tr                   = 15.5
eps_time               = 1e-12
round_time             = 6
recratio               = 0.75 #recurrance ratio
plotOn_ss              = 0;
solverMethod           = "impAdams_d"
saveData               = TRUE #ONLY SAVE SUMMARY TO BE FASTER

source('load_estimations_cluster.R')
source('load_models_cluster.R')

popSizeLog   = 6 
popSize      = 10^popSizeLog
eps          = 1/popSize
eps_in       = 100*popSize

filname_base_metadata = paste0('./SIM_META_SP_H_',type,'_',simIdx)
load(paste0(filname_base_metadata,'.RData'))
print('loaded')
indSampleUseKeep = indSampleUse

# to make it faster

s2p             = s2pMat[ind,1]
p2s             = s2pMat[ind,2]
spcode          = -1*(10*log10(s2p/(32*24))+log10(p2s/(32*24)))
SS              = SS_all[ind,]

filname_base_SS = paste0('./SIM_SS_SP_H_',spcode,'_POPSIZE_',popSizeLog,'_',type,'_',simIdx,'_2Y')
totalNumOfLoops = length(inocSizeVec)*length(trtVec)

out_merged_keep_probiotic  = as.data.frame(matrix(NA, nrow = length(inocSizeVec)*length(trtVec), ncol = 12 + 21 +1))
out_merged_keep_antibiotic = as.data.frame(matrix(NA, nrow = length(trtVec), ncol = 12 + 21 + 1))
out_all_keep_probiotic     = as.data.frame(matrix(NA, nrow = length(inocSizeVec)*length(trtVec), ncol = 12 + 41 + 1))
out_all_keep_antibiotic    = as.data.frame(matrix(NA, nrow = length(trtVec), ncol = 12 + 41 + 1))
counter_row_ant = 1
loopIdx         = 1;

for (numOfTotalCourses in trtVec){
  
  an.error.occured = FALSE
  # RUN ANTIBIOTIC-TREATMENT
  inocSize=NA
  source('run_treatment_antibiotic_cluster.R')
  
  if(!an.error.occured){
    
    tempRow_out_merged_keep_antibiotic = c(simIdx,indSampleUse,inocSize,numOfTotalCourses,s2p,p2s,sta_ext_antibiotic,
                                           sta_min_antibiotic,time2min_antibiotic,time2rec_antibiotic,shannon_antibiotic_final,
                                           R_antibiotic_species,out_antibiotic[dim(out_antibiotic)[1],],popSize)
    
    tempRow_out_all_keep_antibiotic = c(simIdx,indSampleUse,inocSize,numOfTotalCourses,s2p,p2s,sta_ext_antibiotic,
                                        sta_min_antibiotic,time2min_antibiotic,time2rec_antibiotic,shannon_antibiotic_final,
                                        R_antibiotic_species,out_all_antibiotic[dim(out_all_antibiotic)[1],],popSize)
    
    
    inocSizeCounter                 = 0
    counter_row_pro                 = 1
    out_merged_keep_probiotic_temp  = as.data.frame(matrix(NA, nrow = length(inocSizeVec), ncol = 12 + 21 +1))
    out_all_keep_probiotic_temp     = as.data.frame(matrix(NA, nrow = length(inocSizeVec), ncol = 12 + 41 + 1))
    
    for(inocSize in inocSizeVec){
      
      ########### PROBIOTIC-TREATMENT
      source('run_treatment_probiotic_cluster.R')
      
      if(!an.error.occured){
        
        
        out_merged_keep_probiotic_temp[counter_row_pro,] = c(simIdx,indSampleUse,inocSize,numOfTotalCourses,s2p,p2s,sta_ext_probiotic,
                                                             sta_min_probiotic,time2min_probiotic,time2rec_probiotic,shannon_probiotic_final,
                                                             R_probiotic_species,out_probiotic[dim(out_probiotic)[1],],popSize) 
        
        out_all_keep_probiotic_temp[counter_row_pro,]    = c(simIdx,indSampleUse,inocSize,numOfTotalCourses,s2p,p2s,sta_ext_probiotic,
                                                             sta_min_probiotic,time2min_probiotic,time2rec_probiotic,shannon_probiotic_final,
                                                             R_probiotic_species,out_all_probiotic[dim(out_all_probiotic)[1],],popSize)
        
        counter_row_pro = counter_row_pro + 1
        inocSizeCounter = inocSizeCounter + 1
        
        # PLOT AND SAVE #   THIS WAS FOR THE CASE WHERE I SAVED EACH SIM
        filname_base_data = paste0('./SIM_COURSES_',numOfTotalCourses,'_INOCSIZE_',log10(inocSize),
                                   '_SP_',spcode,'_POPSIZE_',popSizeLog,'_',type,'_H_',simIdx,'_2Y')
        plotOnTrt = FALSE
        source('plot_and_save_data_cluster.R')
        
        loopIdx = loopIdx + 1;
        
      }
      
    }
    
    # IF PROBIOTIC TRT WORKED FOR ALL INOC SIZES, SAVE ANTIBIOTIC
    if(inocSizeCounter==length(inocSizeVec)){
      out_merged_keep_antibiotic[counter_row_ant,] =  tempRow_out_merged_keep_antibiotic
      out_all_keep_antibiotic[counter_row_ant,]    = tempRow_out_all_keep_antibiotic
      
      counter_row_ant = counter_row_ant + 1
      
      out_all_keep_probiotic    = rbind(out_all_keep_probiotic,out_all_keep_probiotic_temp)
      out_merged_keep_probiotic = rbind(out_merged_keep_probiotic,out_merged_keep_probiotic_temp)
    }
  }
}

filname_base_data_summary = paste0('./SIM_SUMMARY_SP_H_',spcode,'_POPSIZE_',popSizeLog,'_',type,'_',simIdx,'_2Y.RData')
save(out_merged_keep_antibiotic,out_all_keep_antibiotic,out_merged_keep_probiotic,out_all_keep_probiotic,file=filname_base_data_summary)
print('Summary saved.')



