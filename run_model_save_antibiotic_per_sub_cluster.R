# YOU CAN ALWAYS GET THE FAMILY ORDER FROM 
# filename_params     = paste0('/net/cephfs/data/btepek/NASAL_R/ESTIMATIONS_MCMC/',type,'/MDSINE-source-code-v1.3/output_nose_oligo_th_25/BVS.results.parameters.txt')
# output_params       = read.csv(filename_params,sep="\t")
# output_growth       = output_params %>% filter(parameter_type=='growth_rate')
# family_order        = output_growth$target_taxon #keep same order for interaction matrix
family_order   = c("Carnobacteriaceae","C1","C2","M1","M2","M3","M4","M5",
                   "P1","P2","P3","P4","P5","P6",
                   "Sta1","Sta2","Sta3","Sta4",
                   "Stre1","Stre2")

# I just saved it as 'family_order.RDS' to make it quick
# family_order   = readRDS('/net/cephfs/home/btepek/NASAL_R/family_order.rds')

focalSpecies   = family_order
focalSpecies_p = paste0(family_order,'_p')
focalSpecies_s = paste0(family_order,'_s')
idxFocal       = match(focalSpecies,family_order)

mu_in          = growth_rates
alpha_in       = interaction_matrix
ystart         = as.numeric(SS)
dt_in          = dt

params   = list(mu=mu_in,alpha=alpha_in,preTreatmentTimes=preTreatmentTimes,preTreatmentDeathRates=preTreatmentDeathRates,
                treatmentTimes=treatmentTimes,s2p=s2p,p2s=p2s,pd=pd,pm=pm,eps_in=eps_in,dt_in=dt_in)
for(i in 1:length(family_order)){ #first fix that
  if (ystart[i] < eps) ystart[i] = 0
}
## ADD THE TREATMENT VARIABLES
ystart    = c(ystart,0,0) #pre trt, anti trt
# print(ystart)

times     = seq(from=0,to=round(max(preTreatmentTimes))+fuptimeIn,by=dt) # returns a sequence
times_all = unique(round(sort(c(times,treatmentTimes)),round_time))
times_pre = seq(from=0,to=max(preTreatmentTimes),by=dt) # returns a sequence
times_for = setdiff(times_all,times_pre)
times_after_trt = setdiff(times_for,treatmentTimes)
time_trt_end    = times_after_trt[1]

# # OLD : BUT WORKED
# rootfun_antibiotic <- function (t, y, pars) {
#   y         = round(y,log10(eps_in))
#   root1a    = 1-(((y[1:(2*length(family_order))])<eps) & (y[1:(2*length(family_order))]>0))
#   root1b    = 1-((y[1:(2*length(family_order))]<0))
#   # root1    = 1-(((y[1:(2*length(family_order))])<eps))
#   root2    = abs(t-preTreatmentTimes)
#   root3    = abs(t-treatmentTimes)
#   root4    = abs(t-time_trt_end)
#   rootBoth = c(root1a,root1b,root2,root3,root4)
#   # yAll   = y[1:(2*length(family_order))]
#   # indMin = which(yAll==min(yAll[yAll>0]))
#   # print(c(t,yAll[indMin]))
#   # print(c(root1a[indMin]))
#   # print(c(root1b[indMin]))
#   return(rootBoth)
# }

# NEW
rootfun_antibiotic <- function (t, y, pars) {

  # OPTION 1
  yBacteria  = y[1:(2*length(family_order))]
  root1cInc  = any(yBacteria<0)
  yBacteria  = round(yBacteria,log10(eps_in))
  yNew       = log10(yBacteria[!is.infinite(log10(yBacteria)) & !is.na(log10(yBacteria))])
  root1aInc  = any(yNew<log10(eps))
  root1bInc  = any(is.na(log10(yBacteria)))

  if(root1aInc | root1bInc | root1cInc){
    root1 = 0
  }else{
    root1 = 1
  }

  root2    = abs(t-preTreatmentTimes)
  root3    = abs(t-treatmentTimes)
  root4    = abs(t-time_trt_end)
  rootBoth = c(root1,root2,root3,root4)

  return(rootBoth)
}

eventfun_antibiotic_simple <- function(t, y, parms){ #first and second commensals 
  # root2    = abs(t-preTreatmentTimes) # these lead to differences!!
  # root3    = abs(t-treatmentTimes)
  # root4    = abs(t-time_trt_end) # these lead to differences!!
  
  root2    = abs(t<=preTreatmentTimes[length(preTreatmentTimes)] & t>= preTreatmentTimes[1])
  root3    = abs(t<=treatmentTimes[length(treatmentTimes)] & t>= treatmentTimes[1])
  
  if(is_null(treatmentTimes)){
    root4    = abs(t>=preTreatmentTimes[length(preTreatmentTimes)])
  }else{
    root4    = abs(t>=treatmentTimes[length(treatmentTimes)])
  }
  

  for(i in 1:(2*length(family_order))){
    if (y[i] < eps) y[i] = 0
  }
  
  if(any(root2)){
    y[(2*length(family_order))+1] = 1
  }else{
    y[(2*length(family_order))+1] = 0
  }
  
  if(any(root3)){
    y[(2*length(family_order))+2] = 1
  }else{
    y[(2*length(family_order))+2] = 0
  }
  
  if(any(root4)){
    y[(2*length(family_order))+1] = 0
    y[(2*length(family_order))+2] = 0
  }
  return(y)
}

out = tryCatch({ode(
  func    = nasal.model.antibiotic,
  y       = ystart,
  times   = times_all,
  parms   = params,
  method  = solverMethod,
  rootfun = rootfun_antibiotic,
  events  = list(func = eventfun_antibiotic_simple, root=TRUE),
  atol = 1e-6,
  rtol = 1e-6
)}, error = function(e) {an.error.occured <<- TRUE})

if(!an.error.occured){
  out                   = as.data.frame(out)
  outkeep               = out
  out_all_antibiotic    = out[,1:(dim(out)[2]-2)] #discard binary trt variables
  out_merged_antibiotic = cbind(out[,1],out[,2:(length(family_order)+1)] + out[,length(family_order)+(2:(length(family_order)+1))])
  out_antibiotic        = out_merged_antibiotic
  colnames(out_antibiotic)[1:dim(out_antibiotic)[2]] = c('time',focalSpecies)
  colnames(out_all_antibiotic)[1:dim(out_all_antibiotic)[2]] = c('time',c(focalSpecies_s,focalSpecies_p))
  print('NO errors for the antibiotic treatment, simulations continuing.')
}else{
  print('An error occured for the antibiotic treatment, terminating this set.')
}

