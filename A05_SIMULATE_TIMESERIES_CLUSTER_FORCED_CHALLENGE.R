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
library(rootSolve)

# library(logr)
args <- commandArgs(TRUE)
print(args)
simIdxCheck= as.integer(args[1]) #0-9999
ind        = as.integer(args[2]) #1-4
ntcCheck   = as.integer(args[3]) #0-8
####################################################

speciesNames_sp=c("Carnobacteriaceae_s","C1_s","C2_s","M1_s","M2_s","M3_s","M4_s","M5_s",
                  "P1_s","P2_s","P3_s","P4_s","P5_s","P6_s",
                  "Sta1_s","Sta2_s","Sta3_s","Sta4_s",
                  "Stre1_s","Stre2_s",
                  "Carnobacteriaceae_p",
                  "C1_p","C2_p","M1_p","M2_p","M3_p","M4_p","M5_p",
                  "P1_p","P2_p","P3_p","P4_p","P5_p","P6_p",
                  "Sta1_p","Sta2_p","Sta3_p","Sta4_p",
                  "Stre1_p","Stre2_p")

type         ='FORCED'
popSizePick  = 1e6
inocSizeVec  = 10^c(-1,0,1,2,3)
inocSizeStaphVec = c(0,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3)

s2pVec = 32*24*c(1e-5,1e-4,1e-3) # hourly to monthly (BUT conversion rate shouldn't be more than 1?)
p2sVec = 32*24*c(1e-2,1e-1) # hourly to monthly
s2pMat = expand.grid(s2pVec,p2sVec)
colnames(s2pMat)=c('s2p','p2s')
s2pMat = rbind(c(0,0),s2pMat) #add 0-0 as control                 
s2p             = s2pMat[ind,1]
p2s             = s2pMat[ind,2]
spcode          = -1*(10*log10(s2p/(32*24))+log10(p2s/(32*24)))

# ADD CHALLENGE
inoculumSpecs          = c(15) # This time staph
solverMethod           = "impAdams_d"
saveData               = TRUE
dt                     = 1/32; # 1 day (data is monthly, therefore the interaction and growth parameters -> 1/30 leads to numerical issues!)
eps_time               = 1e-12
round_time             = 6
popSize                = 1e6

source('load_models_cluster.R')

keepAll_antibiotic=c()
keepAll_probiotic=c()

for (incComm in inocSizeVec){
  filename_data = paste0('./SIM_COURSES_',ntcCheck,'_INOCSIZE_',log10(incComm),'_SP_',spcode,'_POPSIZE_6_FORCED_H_',simIdxCheck,'_2Y.RData')
  
  if(file.exists(filename_data)){
    load(filename_data)
    out_all_probiotic_2Y  = out_all_probiotic
    out_all_antibiotic_2Y = out_all_antibiotic
    
    for (inocSizeStaph in inocSizeStaphVec){
      for(fup in c(1,5,15,30,60,90,120)){
        challengeStart         = inoculumTimes[length(inoculumTimes)]+fup*dt
        challengeTime          = challengeStart
        fuptime_ch             = challengeTime+12*30*2 #2 years fup #give a small time and then check steady state
        
        if(dim(out_all_probiotic_2Y)[1]>0 & dim(out_all_antibiotic_2Y)[1]>0){
          
          print(c(dim(out_all_probiotic_2Y)[1],challengeStart))
          out_all_probiotic_FUP  = out_all_probiotic_2Y %>% filter(time > challengeStart)
          out_all_antibiotic_FUP = out_all_antibiotic_2Y %>% filter(time > challengeStart)
          
          yStartProbiotic       = out_all_probiotic_FUP[1,speciesNames_sp]
          ystart                = as.numeric(yStartProbiotic)
          an.error.occured      = FALSE
          source('run_treatment_challenge_cluster_SS.R')
          tempRow_probiotic   = c(simIdxCheck,ntcCheck,popSize,spcode,s2p,p2s,fup,incComm,inocSizeStaph,out_all_challenge)
          if(!an.error.occured){
            keepAll_probiotic   = rbind(keepAll_probiotic,tempRow_probiotic)
          }
          
          yStartAntibiotic     = out_all_antibiotic_FUP[1,speciesNames_sp]
          ystart               = as.numeric(yStartAntibiotic)
          an.error.occured     = FALSE
          source('run_treatment_challenge_cluster_SS.R')
          tempRow_antibiotic   = c(simIdxCheck,ntcCheck,popSize,spcode,s2p,p2s,fup,incComm,inocSizeStaph,out_all_challenge)
          if(!an.error.occured){
            keepAll_antibiotic   = rbind(keepAll_antibiotic,tempRow_antibiotic)
          }
        }
      }
    }
  }
}


filname_base_data_summary = paste0('./SIM_SUMMARY_SP_H_',spcode,'_',ntcCheck,'_',simIdxCheck,'_SS_FUP.RData')
save(keepAll_antibiotic,keepAll_probiotic,file=filname_base_data_summary)
print('SS Summary saved.')






