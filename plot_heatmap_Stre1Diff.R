out_probiotic_agg_all_Inf_use  = aggregate(Stre1 ~ numOfTotalCourses + inocSize + trt, out_probiotic_agg_all_use, FUN= function(z) myfunc(z))
out_antibiotic_agg_all_Inf_use = aggregate(Stre1 ~ numOfTotalCourses  + trt, out_antibiotic_agg_all_use, FUN= function(z) myfunc(z))
out_antibiotic_agg_all_Inf_use$inocSize = NA
out_antibiotic_agg_all_Inf_use = out_antibiotic_agg_all_Inf_use[c(colnames(out_antibiotic_agg_all_Inf_use))]
out_agg_all_Inf                = rbind(out_probiotic_agg_all_Inf_use,out_antibiotic_agg_all_Inf_use)

Stre1_antibiotic_df = out_antibiotic_agg_all_Inf_use[,c('numOfTotalCourses','Stre1')]
colnames(Stre1_antibiotic_df)[2] = 'Stre1_antibiotic'
Stre1_merged = merge(out_agg_all_Inf,Stre1_antibiotic_df,by='numOfTotalCourses')
Stre1_merged =  Stre1_merged %>% rowwise() %>% mutate(Stre1_diff = Stre1 - Stre1_antibiotic)
Stre1_merged_prob = Stre1_merged %>% filter(trt=='probiotic')
inoc_sizes = sort(unique(Stre1_merged_prob$inocSize))
dataPlot   = Stre1_merged_prob[c('numOfTotalCourses','inocSize','Stre1_diff','Stre1')]
dataPlot   = dataPlot %>% filter(inocSize<=inoc_sizes[length(inoc_sizes)])
dataPlot_keep=dataPlot

min_val = min((dataPlot$Stre1_diff[!(is.infinite(dataPlot$Stre1_diff) | is.na(dataPlot$Stre1_diff))]))
max_val = max((dataPlot$Stre1_diff[!(is.infinite(dataPlot$Stre1_diff) | is.na(dataPlot$Stre1_diff))]))
max_val = max(abs(min_val),max_val)

colorRangeRed  <- colorRampPalette(c("#ac0e28","#FFFFFF"))   # red to white
colorRangeBlue <- colorRampPalette(c("#FFFFFF","#013766"))   # white to blue

brk         = sort(do.breaks(c(-(max_val+1), (max_val+1)), N))
paletteSize = length(brk)
colorsRed  <- colorRangeRed((paletteSize-1)/2)                          # Extract 100 color codes
colorsBlue <- colorRangeBlue((paletteSize-1)/2)                         # Extract 100 color codes
colorsAll  = c(colorsRed[1:length(colorsRed)-1],"#FFFFFF",colorsBlue[2:length(colorsBlue)])
cols       = colorsAll

idx     = which(brk<min_val)
idx     = idx[length(idx)]
brk     = brk[(idx):length(brk)]
cols    = c(cols[idx:length(cols)]) #make the negs res


dataPlot_use = dataPlot
# # dataPlot_use[which(is.infinite(dataPlot_use$Stre1_diff)),]$Stre1_diff=1e10
# dataPlot_use[which(dataPlot_use$Stre1_diff==Inf),]$Stre1_diff=1e10
# dataPlot_use[which(dataPlot_use$Stre1_diff==-Inf),]$Stre1_diff=-1e10


hm_temp=levelplot((Stre1_diff) ~ numOfTotalCourses*log10(inocSize), data=dataPlot_use  ,
                  xlab=list(label="Number of Treatment Courses", cex=0.8), ylab=list(label='log10(Relative Inoculum Size)',cex=0.8),
                  main=list(label="Median Final Strep. Difference \n (probiotic - antibiotic)",cex=0.8),
                  col.regions=cols,at = brk, colorkey = list(col = cols, at = brk), par.settings=list(panel.background=list(col="gray")),
                  panel = function(...){
                    panel.levelplot(...)
                    panel.abline(v=(0.5:7.5),col="black")              
                    # panel.abline(h=(-3.5:0.5),col="black")    # G2
                    panel.abline(h=(-3.5:3.5),col="black")  # G3  
                    panel.abline(v=(3.5),col="black")              
                    
                    #sp.polygons(imap)         
                  })

png(paste0('./SS_Stre1_diff_heatmap_combo_0_',
           'popSize_',log(popSizePick,10),'_',type,'_',suffix,'.png'), width = 4, height = 4, units = "in", res =1200)
print(hm_temp)
dev.off()