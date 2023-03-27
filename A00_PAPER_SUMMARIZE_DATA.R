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
library(dplyr)
library(plyr)

dir=setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

infant_data = (read.csv("./ABUNDANCE_DATA.csv") %>% tbl_df())
meta_data   = (read.csv("./META_DATA.csv") %>% tbl_df())
use_data    = (read.csv("./USE_DATA.csv") %>% tbl_df())

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

colNames_all  = colnames(infant_data);
colNames      = colNames_all[2:(length(colNames_all)-4)]; #last 4 cols removed
# FILTER OUT _2 SAMPLES TO PREVENT DUPLICATES
idx2remove    = which(str_detect(colNames, "_", negate = FALSE));
colNames      = colNames[ -c(idx2remove)]
infant_data   = infant_data[,-c(idx2remove+1,dim(infant_data)[2]-3,dim(infant_data)[2]-2,
                                dim(infant_data)[2]-1,dim(infant_data)[2])]#last 4 cols removed
##### SANITY CHECK ##########
which(str_detect(colnames(infant_data), "_", negate = FALSE)) #0
##### SANITY CHECK ##########

rowNames      = infant_data[,1];
families      = as.character(unlist(rowNames));

colNamesSub      = sub(".*X", "",colNames);
infantIndexes    = as.numeric((sub("\\..*","",colNamesSub)));
infantIndexesUnq = as.numeric(unique(sub("\\..*","",colNamesSub)));
weekIndexes      = as.numeric(sub(".*\\.","",colNamesSub));

# TRANSFORM INTO A BETTER TABLE - TRANSPOSE
abundance_data     = tibble();
colnames_abundance = c('infant_idx','week_idx',families);
for(i in seq(1,length(colNames))){
  infant_idx     = infantIndexes[i];
  week_idx       = weekIndexes[i];
  all_abundances = as.numeric(unlist(infant_data[,i+1]));
  print(c(infant_idx,week_idx,colnames(infant_data[,i+1])))
  abundance_data = rbind(abundance_data,c(infant_idx,week_idx,all_abundances))
  colnames(abundance_data) = colnames_abundance;
}
saveRDS(abundance_data, file = "abundance_data.rds")
abundance_data  = readRDS('abundance_data.rds')

# abundance_456   = abundance_data %>% filter(infant_idx==456)
# abundance_456_s = abundance_456[,c(1,2,175)]
# abundance_456_s = abundance_456_s[order(abundance_456_s$week_idx),]

# USE-DATA FIRST COLUMN -> SEPARATE TO INFANT IDX AND WEEK IDX
useID               = as.character(use_data$ID);
useID_infantIndexes = as.numeric((sub("\\-.*","",useID)));
useID_weekIndexes   = as.numeric((sub(".*\\-","",useID)));
infant_idx    = useID_infantIndexes;
week_idx      = useID_weekIndexes;
use_data_all  = tibble(infant_idx=as.numeric(),week_idx=as.numeric());
use_data_all  = add_row(use_data_all,infant_idx,week_idx)

abundance_rows = abundance_data[c('infant_idx','week_idx')];
use_rows       = use_data_all[c('infant_idx','week_idx')];
use_rows       = as.data.frame(use_rows)
matched_rows   = prodlim::row.match(abundance_rows,use_rows,nomatch = NA);
abundance_use  = abundance_data[which(!is.na(matched_rows)==TRUE),];
saveRDS(abundance_use, file = "abundance_use.rds")

# abundance_456   = abundance_use %>% filter(infant_idx==456)
# abundance_456_s = abundance_456[,c(1,2,175)]
# abundance_456_s = abundance_456_s[order(abundance_456_s$week_idx),]

# META-DATA FIRST COLUMN -> SEPARATE TO INFANT IDX AND WEEK IDX
metaID                = as.character(meta_data$ID)
metaID_infantIndexes_1= as.numeric((sub("\\-.*","",metaID)));
metaID_infantIndexes_2= as.numeric((sub("\\_.*","",metaID)));
metaID_infantIndexes_1[is.na(metaID_infantIndexes_1)]=0;
metaID_infantIndexes_2[is.na(metaID_infantIndexes_2)]=0;
metaID_infantIndexes = metaID_infantIndexes_1 + metaID_infantIndexes_2;

metaID_weekIndexes_1= as.numeric((sub(".*\\-","",metaID)));
metaID_weekIndexes_2= as.numeric((sub(".*\\_","",metaID)));
metaID_weekIndexes_1[is.na(metaID_weekIndexes_1)]=0;
metaID_weekIndexes_2[is.na(metaID_weekIndexes_2)]=0;
metaID_weekIndexes = metaID_weekIndexes_1 + metaID_weekIndexes_2;

infant_idx    = metaID_infantIndexes;
week_idx      = metaID_weekIndexes;
meta_data_all = cbind(infant_idx,week_idx,meta_data)
meta_data_all = meta_data_all[c('infant_idx','week_idx','age','season')]

meta_rows           = meta_data_all[c('infant_idx','week_idx')];
abundance_use_rows  = abundance_use[c('infant_idx','week_idx')];
matched_rows        = row.match(abundance_use_rows,meta_rows,nomatch = NA);

# DIFFERENCE   : abundance_use_rows[which(is.na(matched_rows)==TRUE),]
# DIFFERENCE 1 : 456-4 is included in the data to be used, but no match for season and age
# DIFFERENCE 2 : 425-12 is included in the data to be used, but no match for season and age

abundance_use_meta       = abundance_use[which(!is.na(matched_rows)==TRUE),];
# abundance_456   = abundance_use_meta %>% filter(infant_idx==456)
# abundance_456_s = abundance_456[,c(1,2,175)]
# abundance_456_s = abundance_456_s[order(abundance_456_s$week_idx),]

abundance_use_age_season = join(abundance_use_meta,meta_data_all,by=c('infant_idx','week_idx'));
# abundance_456   = abundance_use_age_season %>% filter(infant_idx==456)
# abundance_456_s = abundance_456[,c(1,2,175,dim(abundance_456)[2]-1,dim(abundance_456)[2])]
# abundance_456_s = abundance_456_s[order(abundance_456_s$week_idx),]
# REORDER COLUMNS
abundance_use_age_season=abundance_use_age_season[,c(1,2,dim(abundance_use_age_season)[2]-1,dim(abundance_use_age_season)[2],seq(3,dim(abundance_use_age_season)[2]-2))]
saveRDS(abundance_use_age_season, file = "abundance_use_age_season.rds")

candidates          = c("Corynebacteriaceae", "Moraxellaceae","Staphylococcaceae", "Streptococcaceae", "Carnobacteriaceae");
pick                = 3;
pick_col_idx        = grep(candidates[pick], colnames(abundance_use_age_season));
tbl_pick            = abundance_use_age_season[c('infant_idx','week_idx','age','season',candidates[pick])];
tbl_pick            = tbl_pick %>% dplyr::rename('abundance'=candidates[pick])
tbl_pick_season_age = tbl_pick[c('season','age','abundance')];

# tbl_pick_season_age$season = as.factor(tbl_pick_season_age$season);
# tbl_pick_season_age$age    = as.factor(tbl_pick_season_age$age);

sumstat_season = summarySE(tbl_pick_season_age, measurevar="abundance", groupvars="season")
sumstat_age    = summarySE(tbl_pick_season_age, measurevar="abundance", groupvars="age")

ggplot(sumstat_age, aes(x=age, y=abundance)) + ylim(0,50) +
  geom_errorbar(aes(ymin=abundance-se, ymax=abundance+se), width=.1) +
  geom_line() +
  geom_point() + ggtitle(candidates[pick]) +
  xlab("Age") + ylab("Relative Abundance")

ggplot(sumstat_season, aes(x=season, y=abundance)) + ylim(0,50) +
  geom_errorbar(aes(ymin=abundance-se, ymax=abundance+se), width=.1) +
  geom_line() +
  geom_point() + ggtitle(candidates[pick]) +
  xlab("Season") + ylab("Relative Abundance")

############
# # COMPARE WITH WHAT MARKUS SENT
# meta_data_sort   = meta_data[c('abundance','age','season')];
# sumstat_season_2 = summarySE(meta_data_sort, measurevar="abundance", groupvars="season")
# sumstat_age_2    = summarySE(meta_data_sort, measurevar="abundance", groupvars="age")
# 
# ggplot(sumstat_age_2, aes(x=age, y=abundance)) + ylim(0,50) +
#   geom_errorbar(aes(ymin=abundance-se, ymax=abundance+se), width=.1) +
#   geom_line() +
#   geom_point() + ggtitle(candidates[pick]) +
#   xlab("Age") + ylab("Relative Abundance")


##### PLOT ALL GIVEN SEASON
candidates          = c("Corynebacteriaceae", "Staphylococcaceae", "Streptococcaceae", "Carnobacteriaceae");
pick                = 1:4;
sumstat_season      = c();
for (p in pick){
  tbl_pick            = abundance_use_age_season[c('infant_idx','week_idx','age','season',candidates[p])];
  tbl_pick            = tbl_pick %>% dplyr::rename('abundance'=candidates[p])
  tbl_pick_season_age = tbl_pick[c('season','age','abundance')];
  tbl_pick_season_age$abundance=tbl_pick_season_age$abundance/100
  sumstat_season_temp = summarySE(tbl_pick_season_age, measurevar='abundance',groupvars="season")
  sumstat_season_temp$taxon = candidates[p]
  sumstat_season = rbind(sumstat_season,sumstat_season_temp)
}

colorRange_1 <- colorRampPalette(c("#FFB396","#FF4646"))
colorRange_2 <- colorRampPalette(c("#FF4646","#3B14A7"))

colorsAll_1  <- colorRange_1(2)
colorsAll_2  <- colorRange_2(3)                         

cols = c('#93ABD3',rev(c(colorsAll_1,colorsAll_2[2:length(colorsAll_2)])))

cols = c("#0A97B0","#4CD3C2","#2541B2","#F54291") #sta1, stre1, carno, c1

p=ggplot(sumstat_season, aes(x=season, y=abundance, fill=taxon)) + ylim(0,0.3) +
  geom_errorbar(aes(ymin=abundance-se, ymax=abundance+se,colour = taxon), width=.1) +
  geom_line(aes(colour = taxon)) + scale_color_manual(values=cols)+
  geom_point(aes(colour = taxon)) + ggtitle('') +
  xlab("Season") + ylab("Normalized Abundance")

filname = paste0('./norm_abundance.png')
ggsave(filname,p,width = 4.5, height = 2, dpi = 300, units = "in", device='png')
    

