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

# args <- commandArgs(TRUE)
# print(args)
# popSizeLog = as.integer(args[1]) #3-8
# ind        = as.integer(args[2]) #1-7
# simIdx     = as.integer(args[3]) #0-999

type   = 'FORCED' #'ORIGINAL','FORCED'
s2pVec = 32*24*c(1e-5,1e-4,1e-3) # hourly to monthly (BUT conversion rate shouldn't be more than 1?)
p2sVec = 32*24*c(1e-2,1e-1) # hourly to monthly
s2pMat = expand.grid(s2pVec,p2sVec)
colnames(s2pMat)=c('s2p','p2s')
s2pMat = rbind(c(0,0),s2pMat) #add 0-0 as control

out_merged_keep_probiotic_all = c()
out_merged_keep_antibiotic_all= c()
out_all_keep_probiotic_all    = c()
out_all_keep_antibiotic_all   = c()
popSizeLog                    = 6

for (simIdx in seq(0,999)){
  filename_metadata = paste0('./SIM_META_SP_H_',type,'_',simIdx,'.RData')
  
  for (ind in 1:dim(s2pMat)[1]){
    
    s2p             = s2pMat[ind,1]
    p2s             = s2pMat[ind,2]
    spcode          = -1*(10*log10(s2p/(32*24))+log10(p2s/(32*24)))
    filename_data   = paste0('./SIM_SUMMARY_SP_H_',spcode,'_POPSIZE_',popSizeLog,'_',type,'_',simIdx,'_SS.RData')
    
    if(file.exists(filename_data)){
      # print(paste0(filename_data,' EXISTS'))
      print(simIdx)
      
      load(filename_data)

      out_all_keep_probiotic_all    = rbind(out_all_keep_probiotic_all, out_all_keep_probiotic)
      out_all_keep_antibiotic_all   = rbind(out_all_keep_antibiotic_all, out_all_keep_antibiotic)
    }
  }
}

# add treatment column
out_all_keep_probiotic_all[,dim(out_all_keep_probiotic_all)[2]+1] ='probiotic'
out_all_keep_antibiotic_all[,dim(out_all_keep_antibiotic_all)[2]+1]='antibiotic'

speciesNames=c('Carnobacteriaceae','C1','C2','M1','M2','M3','M4','M5',
               'P1','P2','P3','P4','P5','P6',
               'Sta1','Sta2','Sta3','Sta4',
               'Stre1','Stre2')

speciesNames_sp=c("Carnobacteriaceae_s","C1_s","C2_s","M1_s","M2_s","M3_s","M4_s","M5_s",
                  "P1_s","P2_s","P3_s","P4_s","P5_s","P6_s",
                  "Sta1_s","Sta2_s","Sta3_s","Sta4_s",
                  "Stre1_s","Stre2_s",
                  "Carnobacteriaceae_p","C1_p","C2_p","M1_p","M2_p","M3_p","M4_p","M5_p",
                  "P1_p","P2_p","P3_p","P4_p","P5_p","P6_p",
                  "Sta1_p","Sta2_p","Sta3_p","Sta4_p",
                  "Stre1_p","Stre2_p")

colnames(out_all_keep_probiotic_all) = c('simIdx','indSampleUse','inocSize','numOfTotalCourses','s2p','p2s',speciesNames_sp,'popSize','trt')
colnames(out_all_keep_antibiotic_all) = c('simIdx','indSampleUse','inocSize','numOfTotalCourses','s2p','p2s',speciesNames_sp,'popSize','trt')


out_all_keep_probiotic_all = na.omit(out_all_keep_probiotic_all)
out_all_keep_antibiotic_all$inocSize=-1
out_all_keep_antibiotic_all = na.omit(out_all_keep_antibiotic_all)
out_all_keep_antibiotic_all$inocSize=NA

save(out_all_keep_probiotic_all,file = paste0('./out_all_keep_probiotic_all_H_',type,'_SS.RData'))
save(out_all_keep_antibiotic_all,file= paste0('./out_all_keep_antibiotic_all_H_',type,'_SS.RData'))

length(unique(out_all_keep_probiotic_all$simIdx))
length(unique(out_all_keep_antibiotic_all$simIdx))
