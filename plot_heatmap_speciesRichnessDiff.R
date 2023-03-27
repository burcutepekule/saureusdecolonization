out_probiotic_agg_all_Inf_use  = aggregate(R_species ~ numOfTotalCourses + inocSize + trt, out_probiotic_agg_all_Inf_use_keep, FUN= function(z) myfunc(z))
out_antibiotic_agg_all_Inf_use = aggregate(R_species ~ numOfTotalCourses  + trt, out_antibiotic_agg_all_Inf_use_keep, FUN= function(z) myfunc(z))
out_antibiotic_agg_all_Inf_use$inocSize = NA
out_antibiotic_agg_all_Inf_use = out_antibiotic_agg_all_Inf_use[c(colnames(out_antibiotic_agg_all_Inf_use))]
out_agg_all_Inf                = rbind(out_probiotic_agg_all_Inf_use,out_antibiotic_agg_all_Inf_use)

R_species_antibiotic_df = out_antibiotic_agg_all_Inf_use[,c('numOfTotalCourses','R_species')]
colnames(R_species_antibiotic_df)[2] = 'R_species_antibiotic'
R_species_merged = merge(out_agg_all_Inf,R_species_antibiotic_df,by='numOfTotalCourses')
R_species_merged =  R_species_merged %>% rowwise() %>% mutate(R_species_diff = R_species - R_species_antibiotic)
R_species_merged_prob = R_species_merged %>% filter(trt=='probiotic')
inoc_sizes = sort(unique(R_species_merged_prob$inocSize))
dataPlot   = R_species_merged_prob[c('numOfTotalCourses','inocSize','R_species_diff','R_species')]
dataPlot   = dataPlot %>% filter(inocSize<=inoc_sizes[length(inoc_sizes)])
dataPlot_keep=dataPlot

min_val = min((dataPlot$R_species_diff[!(is.infinite(dataPlot$R_species_diff) | is.na(dataPlot$R_species_diff))]))
max_val = max((dataPlot$R_species_diff[!(is.infinite(dataPlot$R_species_diff) | is.na(dataPlot$R_species_diff))]))
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
# # dataPlot_use[which(is.infinite(dataPlot_use$R_species_diff)),]$R_species_diff=1e10
# dataPlot_use[which(dataPlot_use$R_species_diff==Inf),]$R_species_diff=1e10
# dataPlot_use[which(dataPlot_use$R_species_diff==-Inf),]$R_species_diff=-1e10


hm_temp=levelplot((R_species_diff) ~ numOfTotalCourses*log10(inocSize), data=dataPlot_use  ,
                  xlab=list(label="Number of Treatment Courses", cex=0.8), ylab=list(label='log10(Relative Inoculum Size)',cex=0.8),
                  main=list(label="Median Microbial Diversity Difference \n (probiotic - antibiotic)",cex=0.8),
                  col.regions=cols,at = brk, colorkey = list(col = cols, at = brk), par.settings=list(panel.background=list(col="gray")),
                  panel = function(...){
                    panel.levelplot(...)
                    panel.abline(v=(0.5:7.5),col="black")              
                    # panel.abline(h=(-3.5:0.5),col="black")    # G2
                    panel.abline(h=(-3.5:3.5),col="black")  # G3  
                    panel.abline(v=(3.5),col="black")              
                    
                    #sp.polygons(imap)         
                  })

png(paste0('./SS_R_species_diff_heatmap_combo_0_',
           'popSize_',log(popSizePick,10),'_',type,'_',suffix,'.png'), width = 4, height = 4, units = "in", res =1200)
print(hm_temp)
dev.off()
