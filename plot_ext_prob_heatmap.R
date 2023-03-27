
colorRange_1 <- colorRampPalette(c("#FFB396","#FF4646"))
colorRange_2 <- colorRampPalette(c("#FF4646","#3B14A7"))
colorsAll_1  <- colorRange_1(3)
colorsAll_2  <- colorRange_2(6)                         
cols = c('#93ABD3',rev(c(colorsAll_1,colorsAll_2[2:length(colorsAll_2)])))

col_pro   = cols[3]
col_anti  = cols[7]
col_fill  = cols[1]

out_pro  = out_probiotic_agg_all_use
out_anti = out_antibiotic_agg_all_use

completeSimIdxs_pro  = as.numeric(names(which(table(out_pro$simIdx)==length(unique(out_pro$inocSize))*length(unique(out_pro$numOfTotalCourses)))))
completeSimIdxs_anti = as.numeric(names(which(table(out_anti$simIdx)==length(unique(out_anti$numOfTotalCourses)))))
completeSimsBoth = intersect(completeSimIdxs_pro,completeSimIdxs_anti)
out_pro  = out_pro %>% filter(simIdx %in% completeSimsBoth)
out_anti = out_anti %>% filter(simIdx %in% completeSimsBoth)

unique(table(out_anti$simIdx))
unique(table(out_pro$simIdx))


out_pro_pre_ext  = out_pro %>% filter(numOfTotalCourses==0 & (Sta1==0))
out_anti_pre_ext = out_anti %>% filter(numOfTotalCourses==0 & (Sta1==0))
out_all_challenge_exclude_indSamples = sort(unique(c(unique(out_pro_pre_ext$indSampleUse),unique(out_anti_pre_ext$indSampleUse))))


out_anti      = out_anti %>% filter(!(indSampleUse %in% out_all_challenge_exclude_indSamples))
out_pro       = out_pro %>% filter(!(indSampleUse %in% out_all_challenge_exclude_indSamples))
# out_anti      = out_anti %>% filter(numOfTotalCourses>0)
# out_pro       = out_pro %>% filter(numOfTotalCourses>0)
out_anti      = out_anti %>% dplyr::rowwise() %>% dplyr::mutate(ext_anti  = as.numeric((Sta1)<1e-6))
out_pro       = out_pro %>% dplyr::rowwise() %>% dplyr::mutate(ext_pro  = as.numeric((Sta1)<1e-6))
out_anti_simple = out_anti[c("simIdx","indSampleUse","inocSize","numOfTotalCourses","s2p","p2s","ext_anti")]
out_pro_simple  = out_pro[c("simIdx","indSampleUse","inocSize","numOfTotalCourses","s2p","p2s","ext_pro")]
out_both_simple = merge(out_anti_simple,out_pro_simple,by=c("simIdx","indSampleUse","numOfTotalCourses","s2p","p2s"))
out_both_simple = subset(out_both_simple, select=-c(inocSize.x))
names(out_both_simple)[names(out_both_simple) == 'inocSize.y'] <- 'inocSize'


unique(table(out_both_simple$simIdx))
c(mean(out_both_simple$ext_anti),mean(out_both_simple$ext_pro))
length(unique(out_both_simple$simIdx))

out_agg_pro   = aggregate(ext_pro ~ numOfTotalCourses + inocSize, out_both_simple, FUN=mean)
out_agg_anti  = aggregate(ext_anti ~ numOfTotalCourses , out_both_simple, FUN=mean)
merged        = merge(out_agg_pro,out_agg_anti)

inocVec          = rev(sort(unique(out_agg_pro$inocSize)))
trtVec           = sort(unique(out_agg_anti$numOfTotalCourses))
merged_hypo      = expand.grid(trtVec,inocVec)
colnames(merged_hypo) = colnames(merged)[1:2]
merged_hypo$ext_pro  = NA
merged_hypo$ext_anti = NA
merged_hypo_pro = merged_hypo %>% left_join(merged, by= c("numOfTotalCourses","inocSize")) 
merged_hypo_pro = merged_hypo_pro[c("numOfTotalCourses","inocSize","ext_pro.y","ext_anti.y")]

merged = merged_hypo_pro[order((merged_hypo_pro[,2]), merged_hypo_pro[,1] ),]
merged$inocSize = log10(merged$inocSize)
# need to update
inocVec          = rev(sort(unique(merged$inocSize)))
trtVec           = sort(unique(merged$numOfTotalCourses))
############# # OPTION 1 #############################################

mat_Stre1_use_1       = matrix(merged$ext_pro, nrow = length(inocVec), byrow = TRUE)
keepDim_use_1         = dim(mat_Stre1_use_1)
mat_Stre1_use_1       = unlist(mat_Stre1_use_1)
dim(mat_Stre1_use_1)  = keepDim_use_1
mat_Stre1_use_1       = t(mat_Stre1_use_1)


mat_Stre1_use_2       = matrix(merged$ext_anti, nrow = length(inocVec), byrow = TRUE)
keepDim_use_2         = dim(mat_Stre1_use_2)
mat_Stre1_use_2       = unlist(mat_Stre1_use_2)
dim(mat_Stre1_use_2)  = keepDim_use_2
mat_Stre1_use_2       = t(mat_Stre1_use_2)

### "ncol1" -> (upper-left triangles)  "ncol2" -> (bottom-right triangles)
numRows = dim(mat_Stre1_use_1)[1]
numCols = dim(mat_Stre1_use_1)[2]

### "ncol1" -> (upper-left triangles)  "ncol2" -> (bottom-right triangles)
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

# colorRangeBlue     = colorRampPalette(c("#FFFFFF","#ff6666"))   # white to pink
col_fill='#3090C7'
# colorRangeBlue     = colorRampPalette(c("#FFFFFF",col_fill))   # white to pink
colorRangeBlue     = colorRampPalette(c("#FFE6E8",col_fill))   # white to pink

colorRangeBluePal  = colorRangeBlue(8)  
colorRangeBluePal  = colorRangeBluePal[1:5]

upperlim=max(max(mat_Stre1_use_1[!is.na(mat_Stre1_use_1)]),max(mat_Stre1_use_2[!is.na(mat_Stre1_use_2)]))

ncol1=matrix(map2color(as.vector(mat_Stre1_use_1[,(1:numCols)]),colorRangeBluePal,limits=c(0,upperlim)),nrow=numRows,ncol=numCols)
ncol2=matrix(map2color(as.vector(mat_Stre1_use_2[,(1:numCols)]),colorRangeBluePal,limits=c(0,upperlim)),nrow=numRows,ncol=numCols)

y=seq(1,numCols)
x=seq(1,numRows)


png(paste0('./Split_chance_ext_',
           'popSize_',log(popSizePick,10),'_SP_0_',type,'_SS_',suffix,'.png'), width = 6.6, height = 6.6, units = "in", res =600)

# setHook("grid.newpage", function() pushViewport(viewport(x=0.97,y=0.7,width=0.9, height=0.8, name="vp", just=c("right","top"))), action="prepend")
# par(mfrow=c(1,1))

par(oma=c(1,1,1,0)) # all sides have 3 lines of space
par(mar=c(1,1,1,0))

plot(0,0,xlim=c(-0.1,numRows+1),ylim=c(0,numCols+1), type="n",bty="n" ,xaxt="n",yaxt="n",xlab="",ylab="")
offset=0.5
for(i in 1:numRows){ 
  for(j in 1:numCols){
    if(!is.na(ncol1[i,j])){
      polygon(x[i]-offset+c(0,1,1),y[j]-offset+c(0,0,1),col=ncol1[i,j])
    }else{
      polygon(x[i]-offset+c(0,1,1),y[j]-offset+c(0,0,1),col="#E5E4E2")
    }
    
    if(!is.na(ncol1[i,j])){
      polygon(x[i]-offset+c(0,1,0),y[j]-offset+c(0,1,1),col=ncol2[i,j])
    }else{
      polygon(x[i]-offset+c(0,1,0),y[j]-offset+c(0,1,1),col="#E5E4E2")
    }
  }
}
for(i in 1:numRows){ 
  for(j in 1:numCols){
    if(!is.na(mat_Stre1_use_2[i,j])){
      text(x[i]-0.1,y[j]+0.5*offset,labels=as.character(round(mat_Stre1_use_2[i,j],3)),col=col_anti,font=2,cex=0.9)
    }else{
      text(x[i]-0.1,y[j]+0.5*offset,labels='x',col=col_anti,font =2)
    }
    
    if(!is.na(mat_Stre1_use_2[i,j])){
      text(x[i]+0.1,y[j]-0.5*offset,labels=as.character(round(mat_Stre1_use_1[i,j],3)),col=col_pro,font=2,cex=0.9)
    }else{
      text(x[i]+0.1,y[j]-0.5*offset,labels='x',col=col_pro,font= 2)
    }
  }
}

rowLabels=as.character(inocVec)
colLabels=as.character(trtVec)
axis(side=2, at= 1:numCols,labels=rev(rowLabels), font=0.15, outer=F, pos=c(0.5,0.5,0.5,0.5))
axis(side=1, at= 1:numRows,labels=colLabels,font=0.15, outer=F, pos=c(0.5,0,0,0))
# abline(h=seq(from=offset,to=offset,by=1),lty=1,col="black",lwd=1)
# abline(v=seq(from=offset,to=offset,by=1),lty=1,col="black",lwd=1)
# setHook("grid.newpage", NULL, "replace")

mtext("Number of Treatment Courses", side=1, line=-1, cex=1, col="black", outer=TRUE)
mtext("log10(Relative Inoculum Size)", side=2, line=-1, cex=1, col="black", outer=TRUE)  
mtext("Probability of decolonization after antibiotic treatment.", side=3, line=-3, cex=1,font=2, col=col_anti, outer=TRUE)  
mtext("Probability of decolonization after probiotic treatment.", side=3, line=-2, cex=1, font=2, col=col_pro, outer=TRUE)  

dev.off()
