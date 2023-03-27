# sample from the posterior for growth rates
numSpecies   = length(allICs_colmeans_keep)
mu_in        = c()

for(spc in 1:numSpecies){
  vecTemp = unlist(growthRates_MCMC$growthRates[spc]);
  # mu_in[spc]=sample(vecTemp[(burnin+1):length(vecTemp)], 1, replace = FALSE, prob = NULL)
  vecTemp   = vecTemp[(burnin+1):length(vecTemp)]
  mu_in[spc]= vecTemp[indSampleUse] 
}
growth_rates = c(mu_in,mu_in)

# sample from the posterior for interaction parameters
alpha_in=matrix(nrow = numSpecies, ncol = numSpecies)
for(spc in seq(1,numSpecies*numSpecies)){
  # indexes goes by rows
  # index 20 in the list -> row=20, column=1
  # index 21 in the list -> row=1, column=2
  # index 41 in the list -> row=1, column=3
  
  colNumber = ceiling(spc/numSpecies)
  rowNumber = spc - (colNumber-1)*numSpecies;
  
  vecTemp = unlist(interactionMat_MCMC$interactionMat[spc]);
  
  # alpha_in[rowNumber,colNumber]=sample(vecTemp[(burnin+1):length(vecTemp)], 1, replace = FALSE, prob = NULL)
  
  vecTemp   = vecTemp[(burnin+1):length(vecTemp)]
  alpha_in[rowNumber,colNumber] = vecTemp[indSampleUse] # pick given the simulation index (starts from 0)
}
interaction_matrix = rbind(cbind(alpha_in,alpha_in),cbind(alpha_in,alpha_in))
orderedIC = abs(rnorm(length(allICs_colmeans), mean = allICs_colmeans, sd = 2*allICs_colsds))
orderedIC = orderedIC/sum(orderedIC)
######### ADDING 19TH OF JUNE
orderedIC[1]=0 # No D. Pigrum
######### ADDING 19TH OF JUNE
initial_conditions = c(orderedIC,rep(0,length(orderedIC))) # add persisters - 0


######### ADDING 16TH OF MAY
speciesNames=c("Carnobacteriaceae","C1","C2","M1","M2","M3","M4","M5",
               "P1","P2","P3","P4","P5","P6",
               "Sta1","Sta2","Sta3","Sta4",
               "Stre1","Stre2")

# Carnobacteriaceae -> D. Pigrum

speciesMat   = as.data.frame(expand.grid(speciesNames,speciesNames))

speciesMat$interaction=as.vector(interaction_matrix[1:20,1:20])
colnames(speciesMat)=c('target_taxon','source_taxon','interaction')

speciesMat_focal_1 = speciesMat %>% filter(source_taxon=='Carnobacteriaceae' & target_taxon=='Sta1') #NEG
speciesMat_focal_2 = speciesMat %>% filter(source_taxon=='Carnobacteriaceae' & target_taxon=='Stre1') #NEG
speciesMat_focal_3 = speciesMat %>% filter(source_taxon=='C1' & target_taxon=='Sta1') #NEG
speciesMat_focal_4 = speciesMat %>% filter(source_taxon=='C1' & target_taxon=='Stre1') #NEG
speciesMat_focal_5 = speciesMat %>% filter(source_taxon=='Sta1' & target_taxon=='Stre1') #NEG
speciesMat_focal_6 = speciesMat %>% filter(source_taxon=='Stre1' & target_taxon=='Sta1') #NEG
speciesMat_focal_7 = speciesMat %>% filter(source_taxon=='C1' & target_taxon=='Carnobacteriaceae') #POS


i_1 = speciesMat_focal_1$interaction
i_2 = speciesMat_focal_2$interaction
i_3 = speciesMat_focal_3$interaction
i_4 = speciesMat_focal_4$interaction
i_5 = speciesMat_focal_5$interaction
i_6 = speciesMat_focal_6$interaction
i_7 = speciesMat_focal_7$interaction

cond_1 = i_1<0
cond_2 = i_2<0
cond_3 = i_3<0
cond_4 = i_4<0
cond_5 = i_5<0
cond_6 = i_6<0
cond_7 = i_7>0
cond_8 = abs(i_1)>abs(i_2)
all_conds = c(cond_1,cond_2,cond_3,cond_4,cond_5,cond_6,cond_7,cond_8)
if(all(all_conds)){
  is2use=TRUE
}else{
  is2use=FALSE
}


