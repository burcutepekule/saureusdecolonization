out_merged_keep_probiotic_all_keep  = out_probiotic_agg_all_use
out_merged_keep_antibiotic_all_keep = out_antibiotic_agg_all_use


colorRange_1 <- colorRampPalette(c("#FFB396","#FF4646"))
colorRange_2 <- colorRampPalette(c("#FF4646","#3B14A7"))
colorsAll_1  <- colorRange_1(3)
colorsAll_2  <- colorRange_2(6)                         
cols = c('#93ABD3',rev(c(colorsAll_1,colorsAll_2[2:length(colorsAll_2)])))

############################################################################################################################################################
out_forced_merged_keep_pro = out_merged_keep_probiotic_all_keep %>% filter(popSize==popSizePick)

out_forced_merged_keep_both_plot_all_bar_f_agg_mean_pro  = aggregate(time2rec ~ numOfTotalCourses + inocSize + trt, out_forced_merged_keep_pro, FUN= function(z) myfunc(z))
out_forced_merged_keep_both_plot_all_bar_f_agg_std_pro  = aggregate(time2rec ~ numOfTotalCourses + inocSize + trt, out_forced_merged_keep_pro, FUN= function(z) myfunc2(z))
out_forced_merged_keep_both_plot_all_bar_f_agg_min_pro  = aggregate(time2rec ~ numOfTotalCourses + inocSize + trt, out_forced_merged_keep_pro, FUN= function(z) mymin(z))
out_forced_merged_keep_both_plot_all_bar_f_agg_max_pro  = aggregate(time2rec ~ numOfTotalCourses + inocSize + trt, out_forced_merged_keep_pro, FUN= function(z) mymax(z))
colnames(out_forced_merged_keep_both_plot_all_bar_f_agg_mean_pro)[4]='time2rec_median'
colnames(out_forced_merged_keep_both_plot_all_bar_f_agg_std_pro)[4]='time2rec_std'
colnames(out_forced_merged_keep_both_plot_all_bar_f_agg_min_pro)[4]='time2rec_min'
colnames(out_forced_merged_keep_both_plot_all_bar_f_agg_max_pro)[4]='time2rec_max'

out_forced_merged_keep_both_plot_all_bar_f_agg_pro_1 = merge(out_forced_merged_keep_both_plot_all_bar_f_agg_mean_pro,out_forced_merged_keep_both_plot_all_bar_f_agg_std_pro,
                                                             by=c('trt','inocSize','numOfTotalCourses'))
out_forced_merged_keep_both_plot_all_bar_f_agg_pro_2 = merge(out_forced_merged_keep_both_plot_all_bar_f_agg_min_pro,out_forced_merged_keep_both_plot_all_bar_f_agg_max_pro,
                                                             by=c('trt','inocSize','numOfTotalCourses'))

out_forced_merged_keep_both_plot_all_bar_f_agg_pro   = merge(out_forced_merged_keep_both_plot_all_bar_f_agg_pro_1,out_forced_merged_keep_both_plot_all_bar_f_agg_pro_2,
                                                             by=c('trt','inocSize','numOfTotalCourses'))

# out_forced_merged_keep_both_plot_all_bar_f_agg_pro[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_pro$time2rec_median)),]$time2rec_std=0
# out_forced_merged_keep_both_plot_all_bar_f_agg_pro[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_pro$time2rec_median)),]$time2rec_median=-0.1

out_forced_merged_keep_both_plot_all_bar_f_agg_pro[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_pro$time2rec_median)),]$time2rec_std=0
out_forced_merged_keep_both_plot_all_bar_f_agg_pro[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_pro$time2rec_median)),]$time2rec_median=Inf

out_forced_merged_keep_anti = out_merged_keep_antibiotic_all_keep %>% filter(popSize==popSizePick)

combo_status <- c(
  '31' = "s2p-p2s \n 1e-3 - 1e-1",
  '32' = "s2p-p2s \n 1e-3 - 1e-2",
  '41' = "s2p-p2s \n 1e-4 - 1e-1",
  '42' = "s2p-p2s \n 1e-4 - 1e-2",
  '51' = "s2p-p2s \n 1e-5 - 1e-1",
  '52' = "s2p-p2s \n 1e-5 - 1e-2",
  '0' = "s2p-p2s \n 0 - 0"
)
out_forced_merged_keep_anti      = out_forced_merged_keep_anti %>% rowwise() %>% mutate(combo=-1*(10*log10(s2p/(32*24))+log10(p2s/(32*24))))
out_forced_merged_keep_anti[which(is.infinite(out_forced_merged_keep_anti$combo)),]$combo=0
out_forced_merged_keep_anti = out_forced_merged_keep_anti[, -which(names(out_forced_merged_keep_anti) %in% c('inocSize'))]

out_forced_merged_keep_both_plot_all_bar_f_agg_mean_anti  = aggregate(time2rec ~ numOfTotalCourses + combo , out_forced_merged_keep_anti, FUN= function(z) myfunc(z))
out_forced_merged_keep_both_plot_all_bar_f_agg_std_anti  = aggregate(time2rec ~ numOfTotalCourses + combo , out_forced_merged_keep_anti, FUN= function(z) myfunc2(z))
out_forced_merged_keep_both_plot_all_bar_f_agg_min_anti  = aggregate(time2rec ~ numOfTotalCourses + combo , out_forced_merged_keep_anti, FUN= function(z) mymin(z))
out_forced_merged_keep_both_plot_all_bar_f_agg_max_anti  = aggregate(time2rec ~ numOfTotalCourses + combo , out_forced_merged_keep_anti, FUN= function(z) mymax(z))

colnames(out_forced_merged_keep_both_plot_all_bar_f_agg_mean_anti)[3] = 'time2rec_median'
colnames(out_forced_merged_keep_both_plot_all_bar_f_agg_std_anti)[3]='time2rec_std'
out_forced_merged_keep_both_plot_all_bar_f_agg_anti = merge(out_forced_merged_keep_both_plot_all_bar_f_agg_mean_anti,out_forced_merged_keep_both_plot_all_bar_f_agg_std_anti,by=c('combo','numOfTotalCourses'))
# out_forced_merged_keep_both_plot_all_bar_f_agg_anti[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_anti$time2rec_median)),]$time2rec_std=0
# out_forced_merged_keep_both_plot_all_bar_f_agg_anti[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_anti$time2rec_median)),]$time2rec_median=-0.1

out_forced_merged_keep_both_plot_all_bar_f_agg_anti[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_anti$time2rec_median)),]$time2rec_std=0
out_forced_merged_keep_both_plot_all_bar_f_agg_anti[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_anti$time2rec_median)),]$time2rec_median=Inf

out_forced_merged_keep_both_plot_all_bar_f_agg_anti$combo =factor(out_forced_merged_keep_both_plot_all_bar_f_agg_anti$combo,levels = c('31','32','41','42','51','52','0'))
out_forced_merged_keep_both_plot_all_bar_f_agg_anti = distinct(merge(out_forced_merged_keep_both_plot_all_bar_f_agg_anti,out_forced_merged_keep_anti[c('s2p','p2s','combo')],by='combo'))

a=fct_reorder(out_forced_merged_keep_both_plot_all_bar_f_agg_anti$combo, out_forced_merged_keep_both_plot_all_bar_f_agg_anti$s2p)
levs=levels(a)
srlabels=c()
srlabels[1]='(0,0)'
for(k in 2:length(levs)){
  substr(levs[k],1,1)
  srlabels[k]=paste0('(1e-',substr(levs[k],1,1),',1e-',substr(levs[k],2,2),')')
}


############################################################################################################################################################
################################################################# PLOT #####################################################################################
############################################################################################################################################################

vecTempAnti  = out_forced_merged_keep_both_plot_all_bar_f_agg_anti$time2rec_median+out_forced_merged_keep_both_plot_all_bar_f_agg_anti$time2rec_std
vecTempPro = out_forced_merged_keep_both_plot_all_bar_f_agg_pro$time2rec_median+out_forced_merged_keep_both_plot_all_bar_f_agg_pro$time2rec_std

ymaxAnti = max(vecTempAnti[is.finite(vecTempAnti)])
ymaxPro  = max(vecTempPro[is.finite(vecTempPro)])
ymaxBoth = max(ymaxAnti,ymaxPro)+1

cols_back=rep('grey90',length(cols))

out_forced_merged_keep_both_plot_all_bar_f_agg_anti_fwd=out_forced_merged_keep_both_plot_all_bar_f_agg_anti
out_forced_merged_keep_both_plot_all_bar_f_agg_anti_fwd[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_anti_fwd$time2rec_median)),]$time2rec_std=0
out_forced_merged_keep_both_plot_all_bar_f_agg_anti_fwd[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_anti_fwd$time2rec_median)),]$time2rec_median=0

data_back = out_forced_merged_keep_both_plot_all_bar_f_agg_anti %>% mutate(combo = fct_reorder(combo, s2p))
data_fwd  = out_forced_merged_keep_both_plot_all_bar_f_agg_anti_fwd %>% mutate(combo = fct_reorder(combo, s2p)) 

barplot_temp_back = ggplot(data_back,aes(x = as.factor(numOfTotalCourses), y = time2rec_median, fill=factor(combo))) +
  geom_bar(stat = "identity", position="dodge",width=0.8) + ylim(NA, ymaxBoth) +
  geom_errorbar(aes(x=as.factor(numOfTotalCourses), ymin=time2rec_median-time2rec_std, ymax=time2rec_median+time2rec_std), width=0.5,
                colour="black", alpha=1, size=0.5,  position= position_dodge(width = 0.8))+
  # geom_hline(data=out_forced_merged_keep_both_plot_all,aes(yintercept=stre1_ss), linetype="dashed", color = "black", size=0.5)+
  scale_fill_manual(name="Switching rates \n(s2p,p2s)",values = cols_back,   labels = srlabels)+
  xlab("Number of treatment courses (1 course = 5 days)")+ylab("Time to S. aureus \nrecurrence (months)")+
  # facet_wrap(~ inocSize, labeller = labeller(incoulum_size = inocSize), nrow = 1) +
  # facet_grid(. ~combo, labeller = labeller(combo = combo_status))+
  theme_bw()+theme( panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA))

barplot_temp_fwd = ggplot(data=data_fwd,aes(x = as.factor(numOfTotalCourses),y = time2rec_median, fill=factor(combo))) +
  geom_bar(stat = "identity", position="dodge",width=0.8) + ylim(NA, ymaxBoth) +
  geom_errorbar(aes(x=as.factor(numOfTotalCourses), ymin=time2rec_median-time2rec_std, ymax=time2rec_median+time2rec_std), width=0.5,
                colour="black", alpha=1, size=0.5,  position= position_dodge(width = 0.8))+
  # geom_hline(data=out_forced_merged_keep_both_plot_all,aes(yintercept=stre1_ss), linetype="dashed", color = "black", size=0.5)+
  scale_fill_manual(name="Switching rates \n(s2p,p2s)",values = cols,   labels = srlabels)+
  xlab("Number of treatment courses (1 course = 5 days)")+ylab("Time to S. aureus \nrecurrence (months)")+
  # facet_wrap(~ inocSize, labeller = labeller(incoulum_size = inocSize), nrow = 1) +
  # facet_grid(. ~combo, labeller = labeller(combo = combo_status))+
  theme_bw()+theme( panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid = element_blank())


filname = paste0('./SS_time2rec_inocSize_antibiotic_grid_',
                 'popSize_',log(popSizePick,10),'_',type,'_fwd_',suffix,'.png')
ggsave(filname,barplot_temp_fwd,width = 8.5, height = 2.5, dpi = 300, units = "in", device='png')

filname = paste0('./SS_time2rec_inocSize_antibiotic_grid_',
                 'popSize_',log(popSizePick,10),'_',type,'_back_',suffix,'.png')
ggsave(filname,barplot_temp_back,width = 8.5, height = 2.5, dpi = 300, units = "in", device='png')


out_forced_merged_keep_both_plot_all_bar_f_agg_pro_fwd=out_forced_merged_keep_both_plot_all_bar_f_agg_pro
out_forced_merged_keep_both_plot_all_bar_f_agg_pro_fwd[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_pro_fwd$time2rec_median)),]$time2rec_std=0
out_forced_merged_keep_both_plot_all_bar_f_agg_pro_fwd[which(is.infinite(out_forced_merged_keep_both_plot_all_bar_f_agg_pro_fwd$time2rec_median)),]$time2rec_median=0

data_back = out_forced_merged_keep_both_plot_all_bar_f_agg_pro
data_fwd  = out_forced_merged_keep_both_plot_all_bar_f_agg_pro_fwd

barplot_temp_back = ggplot(data_back,aes(x = as.factor(numOfTotalCourses), y = time2rec_median, fill=factor(inocSize))) +
  geom_bar(stat = "identity", position="dodge",width=0.8) + ylim(NA, ymaxBoth) +
  geom_errorbar(aes(x=as.factor(numOfTotalCourses), ymin=time2rec_median-time2rec_std, ymax=time2rec_median+time2rec_std), width=0.5,
                colour="black", alpha=1, size=0.5,  position= position_dodge(width = 0.8))+
  # geom_hline(data=out_forced_merged_keep_both_plot_all,aes(yintercept=stre1_ss), linetype="dashed", color = "black", size=0.5)+
  scale_fill_manual(name="Normalized \ninoculum size",values = cols_back)+
  xlab("Number of treatment courses (1 course = 5 days)")+ylab("Time to S. aureus \nrecurrence (months)")+
  # facet_wrap(~ inocSize, labeller = labeller(incoulum_size = inocSize), nrow = 1) +
  # facet_grid(. ~combo, labeller = labeller(combo = combo_status))+
  theme_bw()+theme( panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA))

barplot_temp_fwd = ggplot(data=data_fwd,aes(x = as.factor(numOfTotalCourses),y = time2rec_median, fill=factor(inocSize))) +
  geom_bar(stat = "identity", position="dodge",width=0.8) + ylim(NA, ymaxBoth) +
  geom_errorbar(aes(x=as.factor(numOfTotalCourses), ymin=time2rec_median-time2rec_std, ymax=time2rec_median+time2rec_std), width=0.5,
                colour="black", alpha=1, size=0.5,  position= position_dodge(width = 0.8))+
  # geom_hline(data=out_forced_merged_keep_both_plot_all,aes(yintercept=stre1_ss), linetype="dashed", color = "black", size=0.5)+
  scale_fill_manual(name="Normalized \ninoculum size",values = cols)+
  xlab("Number of treatment courses (1 course = 5 days)")+ylab("Time to S. aureus \nrecurrence (months)")+
  # facet_wrap(~ inocSize, labeller = labeller(incoulum_size = inocSize), nrow = 1) +
  # facet_grid(. ~combo, labeller = labeller(combo = combo_status))+
  theme_bw()+theme( panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid = element_blank())


filname = paste0('./SS_time2rec_inocSize_probiotic_grid_',
                 'popSize_',log(popSizePick,10),'_',type,'_fwd_',suffix,'.png')
ggsave(filname,barplot_temp_fwd,width = 8.5, height = 2.5, dpi = 300, units = "in", device='png')

filname = paste0('./SS_time2rec_inocSize_probiotic_grid_',
                 'popSize_',log(popSizePick,10),'_',type,'_back_',suffix,'.png')
ggsave(filname,barplot_temp_back,width = 8.5, height = 2.5, dpi = 300, units = "in", device='png')


# 
# barplot_temp = out_forced_merged_keep_both_plot_all_bar_f_agg_pro %>%
#   ggplot(aes(x = as.factor(numOfTotalCourses), y = time2rec_median, fill=factor(inocSize))) +
#   geom_bar(stat = "identity", position="dodge",width=0.8)  + ylim(NA, ymaxBoth) +
#   geom_errorbar(aes(x=as.factor(numOfTotalCourses), ymin=time2rec_median-time2rec_std, ymax=time2rec_median+time2rec_std), width=0.5,
#                 # geom_errorbar(aes(x=as.factor(numOfTotalCourses), ymin=time2rec_min, ymax=time2rec_max), width=0.5,
#                 colour="black", alpha=1, size=0.5,  position= position_dodge(width = 0.8))+
#   # geom_hline(data=out_forced_merged_keep_both_plot_all,aes(yintercept=stre1_ss), linetype="dashed", color = "black", size=0.5)+
#   scale_fill_manual(name="Normalized \ninoculum size",values = cols)+
#   xlab("Number of treatment courses (1 course = 5 days)")+ylab("Time to S. aureus \nrecurrence (months)")+
#   # facet_wrap(~ inocSize, labeller = labeller(incoulum_size = inocSize), nrow = 1) +
#   # facet_grid(. ~combo, labeller = labeller(combo = combo_status))+
#   theme_bw()
# 
# 
# filname = paste0('./SS_time2rec_inocSize_probiotic_grid_',
#                  'popSize_',log(popSizePick,10),'_',type,'_',suffix,'.png')
# ggsave(filname,barplot_temp,width = 8.5, height = 2.5, dpi = 300, units = "in", device='png')
