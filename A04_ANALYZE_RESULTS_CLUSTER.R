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

myfunc_1 <- function(vec){
  if(all(is.infinite(vec))){
    return(Inf)
  }else{
    return(mean(vec[!is.infinite(vec)]))
  }
}

myfunc_2 <- function(vec){
  if(length(vec[is.infinite(vec)])>length(vec[!is.infinite(vec)])){
    return(Inf)
  }else{
    return(mean(vec[!is.infinite(vec)]))
  }
}

myfunc_3 <- function(vec){
  return(median(vec))
}

myfunc_4 <- function(vec){
  return(sd(vec[!is.infinite(vec)]))
}

mymin<- function(vec){
  return(min(vec[!is.infinite(vec)]))
}

mymax<- function(vec){
  return(max(vec[!is.infinite(vec)]))
}
# pick
myfunc=myfunc_3
myfunc2=myfunc_4

popSizePick =1e6
type='FORCED'
suffix='H' #this is the version index for the code

load(paste0('./out_merged_keep_probiotic_all_H_',type,'.RData'))
out_merged_keep_probiotic_all_keep  = out_merged_keep_probiotic_all

load(paste0('./out_merged_keep_antibiotic_all_H_',type,'.RData'))
out_merged_keep_antibiotic_all_keep  = out_merged_keep_antibiotic_all

### FOR s2p==0 & p2s==0

out_antibiotic_agg_all_Inf = out_merged_keep_antibiotic_all_keep %>% filter(s2p==0 & p2s==0 & popSize==popSizePick)
out_antibiotic_agg_all_Inf = out_antibiotic_agg_all_Inf[c('simIdx','indSampleUse','numOfTotalCourses','sta_ext','time2rec','R_species','trt','popSize','s2p','p2s')]

out_probiotic_agg_all_Inf = out_merged_keep_probiotic_all_keep %>% filter(s2p==0 & p2s==0 & popSize==popSizePick)
out_probiotic_agg_all_Inf = out_probiotic_agg_all_Inf[c('simIdx','indSampleUse','numOfTotalCourses','inocSize','sta_ext','time2rec','R_species','trt','popSize','s2p','p2s')]

freqOfInds_pro    = table(out_probiotic_agg_all_Inf$indSampleUse)/(length(unique(out_probiotic_agg_all_Inf$numOfTotalCourses))*length(unique(out_probiotic_agg_all_Inf$inocSize)))
inds2exclude_pro  = as.numeric(names(which(!(freqOfInds_pro %% 1 == 0))))
freqOfInds_anti   = table(out_antibiotic_agg_all_Inf$indSampleUse)/(length(unique(out_antibiotic_agg_all_Inf$numOfTotalCourses))*length(unique(out_antibiotic_agg_all_Inf$inocSize)))
inds2exclude_anti = as.numeric(names(which(!(freqOfInds_anti %% 1 == 0))))
inds2exclude_both = c(inds2exclude_pro,inds2exclude_anti)

######## PLOTTING - NEEDS MODIFICATION
out_antibiotic_agg_all_Inf$inocSize = NA
out_antibiotic_agg_all_Inf = out_antibiotic_agg_all_Inf[c(colnames(out_probiotic_agg_all_Inf))]
out_agg_all_Inf = rbind(out_probiotic_agg_all_Inf,out_antibiotic_agg_all_Inf)
out_agg_all_Inf = out_agg_all_Inf %>% filter(!indSampleUse %in% inds2exclude_both)

# discard the ones that died in the antibiotic pretreatment?
out_agg_all_Inf_check   = out_agg_all_Inf %>% filter(numOfTotalCourses==0 & time2rec == Inf)
inds_discard_probiotic  = out_agg_all_Inf_check %>% filter(trt=='probiotic')
inds_discard_antibiotic = out_agg_all_Inf_check %>% filter(trt=='antibiotic')
indSampleNotUseBoth = unique(c(unique(inds_discard_probiotic$indSampleUse),unique(inds_discard_antibiotic$indSampleUse)))
out_probiotic_agg_all_Inf_use   = out_probiotic_agg_all_Inf %>% filter(! (indSampleUse %in% indSampleNotUseBoth))
out_antibiotic_agg_all_Inf_use  = out_antibiotic_agg_all_Inf %>% filter(! (indSampleUse %in% indSampleNotUseBoth))

useSamples1 = sort(unique(out_probiotic_agg_all_Inf_use$indSampleUse))
useSamples2 = sort(unique(out_antibiotic_agg_all_Inf_use$indSampleUse))

useSims1 = sort(unique(out_probiotic_agg_all_Inf_use$simIdx))
useSims2 = sort(unique(out_antibiotic_agg_all_Inf_use$simIdx))

# take the avg first
out_probiotic_agg_all_Inf_use_keep  = out_probiotic_agg_all_Inf_use
out_antibiotic_agg_all_Inf_use_keep = out_antibiotic_agg_all_Inf_use

out_probiotic_agg_all_Inf_use_keep_0 = out_probiotic_agg_all_Inf_use_keep %>% filter(numOfTotalCourses==0)
out_antibiotic_agg_all_Inf_use_keep_0 = out_antibiotic_agg_all_Inf_use_keep %>% filter(numOfTotalCourses==0)

out_probiotic_agg_all_Inf_use_keep_0_d = distinct(out_probiotic_agg_all_Inf_use_keep_0[c( "simIdx","indSampleUse","numOfTotalCourses","sta_ext","time2rec","R_species")])
out_antibiotic_agg_all_Inf_use_keep_0_d= distinct(out_antibiotic_agg_all_Inf_use_keep_0[c( "simIdx","indSampleUse","numOfTotalCourses","sta_ext","time2rec","R_species")])

colnames(out_probiotic_agg_all_Inf_use_keep_0_d)[4:6]=paste0(colnames(out_probiotic_agg_all_Inf_use_keep_0_d)[4:6],'_pro')
colnames(out_antibiotic_agg_all_Inf_use_keep_0_d)[4:6]=paste0(colnames(out_antibiotic_agg_all_Inf_use_keep_0_d)[4:6],'_anti')

out_both_agg_all_Inf_use_keep_0_d = merge(out_probiotic_agg_all_Inf_use_keep_0_d,out_antibiotic_agg_all_Inf_use_keep_0_d,by=c("simIdx","indSampleUse","numOfTotalCourses"))

N = 27
source('plot_heatmap_time2recDiff.R')
N = 21
source('plot_heatmap_speciesRichnessDiff.R')


load(paste0('./out_merged_keep_probiotic_all_H_',type,'.RData'))
out_merged_keep_probiotic_all_keep  = out_merged_keep_probiotic_all

load(paste0('./out_merged_keep_antibiotic_all_H_',type,'.RData'))
out_merged_keep_antibiotic_all_keep = out_merged_keep_antibiotic_all
unique(table(out_merged_keep_antibiotic_all_keep$simIdx))

out_probiotic_agg_all_use   = out_merged_keep_probiotic_all_keep %>% filter((indSampleUse %in% useSamples1))
out_probiotic_agg_all_use   = na.omit(out_probiotic_agg_all_use)

out_antibiotic_agg_all_use  = out_merged_keep_antibiotic_all_keep %>% filter((indSampleUse %in% useSamples1))
out_antibiotic_agg_all_use$inocSize=0
out_antibiotic_agg_all_use=na.omit(out_antibiotic_agg_all_use)
out_antibiotic_agg_all_use$inocSize=NA

source('plot_barcombo_time2recDiff.R')

out_antibiotic_agg_all_use = out_antibiotic_agg_all_use %>% filter(s2p==0 & p2s==0 & popSize==popSizePick)
out_probiotic_agg_all_use  = out_probiotic_agg_all_use %>% filter(s2p==0 & p2s==0 & popSize==popSizePick)

N = 21
source('plot_heatmap_Stre1Diff.R')

source('plot_ext_prob_heatmap.R')




