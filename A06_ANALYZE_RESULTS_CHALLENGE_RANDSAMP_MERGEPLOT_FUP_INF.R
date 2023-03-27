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


filename="./ANALYZE_CHALLENGE_H_INF.RData"

if(!file.exists(filename)){
  print('File doesnt exists, will generate')

  keepData_pro_keep  = readRDS(paste0('./keepData_CH_samples_sim_ALL_FUP_PRO_ALL_FILTERED_H_INF.rds')) # for all combos
  keepData_anti_keep = readRDS(paste0('./keepData_CH_samples_sim_ALL_FUP_ANTI_ALL_FILTERED_H_INF.rds'))
  
  c(length(unique(keepData_pro_keep$simIdx)),length(unique(keepData_anti_keep$simIdx)))
  # FROM HEATMAP PLOTTING SCRIPT G8095_ANALYZE_RESULTS_CLUSTER
  useSims1    = readRDS('./useSims_H.rds')
  useSamples1 = readRDS('./useSamples_H.rds')
  keepData_pro_keep=keepData_pro_keep%>%filter(simIdx %in% useSims1)
  keepData_anti_keep=keepData_anti_keep%>%filter(simIdx %in% useSims1)
  c(length(unique(keepData_pro_keep$simIdx)),length(unique(keepData_anti_keep$simIdx))) #954
  
  keepData_pro_keep_use  = keepData_pro_keep %>% filter(inocStaph>0)
  keepData_anti_keep_use = keepData_anti_keep %>% filter(inocStaph>0)
  
  keepData_pro_keep_use  = keepData_pro_keep_use[c("simIdx","ntc","spcode","fup","incComm","inocStaph","Treatment","ext_pro_sta","ext_pro_stre")]
  keepData_anti_keep_use = keepData_anti_keep_use[c("simIdx","ntc","spcode","fup","incComm","inocStaph","Treatment","ext_anti_sta","ext_anti_stre")]
  colnames(keepData_pro_keep_use)[8:9]=c('ext_sta','ext_stre')
  colnames(keepData_anti_keep_use)[8:9]=c('ext_sta','ext_stre')
  
  keepData_both = rbind(keepData_pro_keep_use,keepData_anti_keep_use)
  
  keepData_both_plot_agg_all = aggregate(ext_sta~ntc+fup+incComm+inocStaph+Treatment+ext_stre,keepData_both,FUN=mean)
  keepData_both_plot_agg_trt = aggregate(ext_sta~Treatment,keepData_both,FUN=mean)
  keepData_both_plot_agg_fup = aggregate(ext_sta~Treatment+fup,keepData_both,FUN=mean)
  keepData_both_plot_agg_ntc = aggregate(ext_sta~Treatment+ntc,keepData_both,FUN=mean)
  keepData_both_plot_agg_incc = aggregate(ext_sta~Treatment+incComm,keepData_both,FUN=mean)
  keepData_both_plot_agg_incs = aggregate(ext_sta~Treatment+inocStaph,keepData_both,FUN=mean)
  keepData_both_plot_agg_stre = aggregate(ext_sta~Treatment+ext_stre,keepData_both,FUN=mean)
  
  rm(keepData_anti_keep_use)
  rm(keepData_pro_keep_use)
  rm(keepData_anti_keep)
  rm(keepData_pro_keep)
  save.image(file = "./ANALYZE_CHALLENGE_H_INF.RData")
}else{
  load(filename)
}
###############################################################################################
length(unique(keepData_both$simIdx)) #954

###############################################################################################
library(ggtext)
library(glue)
# y_sub = expression(paste("Mean ", italic("S. aureus")))
# my_y_title = glue("Mean <i>S. aureus</i> </b><br> extinction probability (%)")
my_y_title <- "Mean S. aureus \nextinction probability (%)"
###############################################################################################

plot = ggplot(keepData_both_plot_agg_trt, aes(x =  Treatment, y = 100*ext_sta, fill = Treatment)) + geom_bar(stat = "identity",width=0.4,alpha=0.8,position = "dodge") + theme_light()+scale_fill_brewer(palette="Set1")+
  labs(y = my_y_title)+theme(axis.text=element_text(size=4))+
  labs(x = "Treatment")

filname=paste0('./ext_prob_ch_TRT_H_INF.png')
ggsave(filname,plot,width = 4, height = 2, dpi = 300, units = "in", device='png')

my_x_title <- expression(paste(italic("S. pneumoniae"), " extinction"))

plot = ggplot(keepData_both_plot_agg_stre, aes(x =  ext_stre, y = 100*ext_sta, fill = Treatment)) + geom_bar(stat = "identity",width=0.4,alpha=0.8,position = "dodge") + 
  theme_light()+scale_fill_brewer(palette="Set1")+theme(axis.text=element_text(size=4))+
  labs(y = my_y_title)+ 
  labs(x = my_x_title)

filname=paste0('./ext_prob_ch_STRE_H_INF.png')
ggsave(filname,plot,width = 4, height = 2, dpi = 300, units = "in", device='png')

