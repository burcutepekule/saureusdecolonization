
out_probiotic_agg_all_Inf_use  = aggregate(time2rec ~ numOfTotalCourses + inocSize + trt, out_probiotic_agg_all_Inf_use_keep, FUN= function(z) myfunc(z))
out_antibiotic_agg_all_Inf_use = aggregate(time2rec ~ numOfTotalCourses  + trt, out_antibiotic_agg_all_Inf_use_keep, FUN= function(z) myfunc(z))
out_antibiotic_agg_all_Inf_use$inocSize = NA
out_antibiotic_agg_all_Inf_use = out_antibiotic_agg_all_Inf_use[c(colnames(out_antibiotic_agg_all_Inf_use))]
out_agg_all_Inf                = rbind(out_probiotic_agg_all_Inf_use,out_antibiotic_agg_all_Inf_use)

time2rec_antibiotic_df = out_antibiotic_agg_all_Inf_use[,c('numOfTotalCourses','time2rec')]
colnames(time2rec_antibiotic_df)[2] = 'time2rec_antibiotic'
time2rec_merged = merge(out_agg_all_Inf,time2rec_antibiotic_df,by='numOfTotalCourses')
time2rec_merged =  time2rec_merged %>% rowwise() %>% mutate(time2rec_diff = time2rec - time2rec_antibiotic)
time2rec_merged_prob = time2rec_merged %>% filter(trt=='probiotic')
inoc_sizes = sort(unique(time2rec_merged_prob$inocSize))
dataPlot   = time2rec_merged_prob[c('numOfTotalCourses','inocSize','time2rec_diff','time2rec')]
dataPlot   = dataPlot %>% filter(inocSize<=inoc_sizes[length(inoc_sizes)])
dataPlot_keep=dataPlot

min_val = min((dataPlot$time2rec_diff[!(is.infinite(dataPlot$time2rec_diff) | is.na(dataPlot$time2rec_diff))]))
max_val = max((dataPlot$time2rec_diff[!(is.infinite(dataPlot$time2rec_diff) | is.na(dataPlot$time2rec_diff))]))
max_val = max(abs(min_val),max_val)

colorRangeRed  <- colorRampPalette(c("#ac0e28","#FFFFFF"))   # red to white
colorRangeBlue <- colorRampPalette(c("#FFFFFF","#013766"))   # white to blue

a           = 2
brk         = sort(do.breaks(c(-(max_val+1), (max_val+1)), N))
paletteSize = length(brk)
colorsRed  <- colorRangeRed(a*(paletteSize-1)/2 - 1)                          # Extract 100 color codes
colorsBlue <- colorRangeBlue(a*(paletteSize-1)/2 - 1)                         # Extract 100 color codes
colorsAll  = c(colorsRed[1:length(colorsRed)-1],"#FFFFFF",colorsBlue[2:length(colorsBlue)])
cols_l     = colorsAll
cols       = cols_l[(N-N/a):(N+N/a)]

idx     = which(brk<min_val)
idx     = idx[length(idx)]
brk     = brk[(idx):length(brk)]
cols    = c(cols[idx:length(cols)]) #make the negs res

brk        = c(-Inf,brk,Inf)
cols       = c(cols_l[1],cols,cols_l[length(cols_l)])

dataPlot_use = dataPlot
# dataPlot_use[which(is.infinite(dataPlot_use$time2rec_diff)),]$time2rec_diff=1e10
dataPlot_use[which(dataPlot_use$time2rec_diff==Inf),]$time2rec_diff=1e10
dataPlot_use[which(dataPlot_use$time2rec_diff==-Inf),]$time2rec_diff=-1e10


hm_temp=levelplot((time2rec_diff) ~ numOfTotalCourses*log10(inocSize), data=dataPlot_use  ,
                  xlab=list(label="Number of Treatment Courses", cex=0.8), ylab=list(label='log10(Relative Inoculum Size)',cex=0.8),
                  main=list(label="Median Time Difference for Recurrence (months) \n (probiotic - antibiotic, 75% of initial abundance)",cex=0.8),
                  col.regions=cols,at = brk, colorkey = list(col = cols, at = brk), par.settings=list(panel.background=list(col="gray")),
                  panel = function(...){
                    panel.levelplot(...)
                    panel.abline(v=(0.5:7.5),col="black")              
                    # panel.abline(h=(-3.5:0.5),col="black")    #G2
                    panel.abline(h=(-3.5:3.5),col="black")    #G3
                    panel.abline(v=(3.5),col="black")              
                    
                    #sp.polygons(imap)         
                  })

png(paste0('./SS_time2rec_diff_heatmap_combo_0_',
           'popSize_',log(popSizePick,10),'_',type,'_',suffix,'.png'), width = 4, height = 4, units = "in", res =1200)
print(hm_temp)
dev.off()
