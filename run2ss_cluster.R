########### NO-TREATMENT - FIND SS
fuptime                     = fuptime_ss #say 24mo.
numOfTotalCourses           = 0
d                           = d_tr; # here to actually kill
preTreatmentDeathRates      = rep(d,20) # INITIAL TREATMENT - REDUCE ICs 99%
preTreatmentDeathRates[2:3] = 0 # Corynebacteriaceae don't die
inoculumSizes_abs           = c(0,0) #if no inoculum, set to zero


#convert days to months
preTreatmentStartIn    = preTreatmentStart*dt 
proTreatmentStartIn    = proTreatmentStart*dt 
addOfPreTtreatmentDays = numOfPreTtreatmentDays*dt;
addOfProTtreatmentDays = numOfProTtreatmentDays*dt;
addOfBreakDays         = numOfBreakDays*dt;
fuptimeIn              = fuptime*dt

#####################################################################################################
numCourses         = 0
preTreatmentTimes  = c()
treatmentTimes     = c()

family_order   = c("Carnobacteriaceae","C1","C2","M1","M2","M3","M4","M5",
                   "P1","P2","P3","P4","P5","P6",
                   "Sta1","Sta2","Sta3","Sta4",
                   "Stre1","Stre2")

focalSpecies   = family_order
focalSpecies_p = paste0(family_order,'_p')
focalSpecies_s = paste0(family_order,'_s')
idxFocal       = match(focalSpecies,family_order)
dt_in          = dt
# # HERE SET s2p=0 and p2s=0 because p2s>s2p always and therefore no persistence will occur without treatment pressure
# params         = list(mu=growth_rates,alpha=interaction_matrix,preTreatmentTimes=preTreatmentTimes,preTreatmentDeathRates=preTreatmentDeathRates,
#                       treatmentTimes=treatmentTimes,s2p=0,p2s=0,pd=pd,pm=pm,eps_in=eps_in,dt_in=dt_in)

ystart         = initial_conditions

for(i in 1:(2*length(family_order))){ #first fix that
  if (ystart[i] < eps) ystart[i] = 0
}
## ADD THE TREATMENT VARIABLES
ystart    = c(ystart,0,0) #pre trt, anti trt
times     = seq(from=0,to=fuptimeIn,by=dt) # returns a sequence
times_all = unique(round(sort(c(times)),round_time))
time_trt_end = c()

SS_all = c()
for (ind in 1:7){
  
  print(ind)
  
  s2p             = s2pMat[ind,1]
  p2s             = s2pMat[ind,2]
  spcode          = -1*(10*log10(s2p/(32*24))+log10(p2s/(32*24)))
  
  params = list(mu=growth_rates,alpha=interaction_matrix,preTreatmentTimes=preTreatmentTimes,preTreatmentDeathRates=preTreatmentDeathRates,
                treatmentTimes=treatmentTimes,s2p=s2p,p2s=p2s,pd=pd,pm=pm,eps_in=eps_in,dt_in=dt_in)
  
  out=tryCatch({runsteady(
    func    = nasal.model.antibiotic,
    y       = ystart,
    times   = c(0,Inf),
    parms   = params,
    # method  = solverMethod,
    # rootfun = rootfun_antibiotic_notrt,
    # events  = list(func = eventfun_antibiotic_simple_notrt, root=TRUE),
    # atol = 1e-6,
    # rtol = 1e-6,
    stol=1e-8,
    mf=10
  )}, error = function(e) {an.error.occured <<- TRUE})
  
  if(!an.error.occured){
    
    out            = as.data.frame(t(as.data.frame(out)))
    out            = out[,1:(dim(out)[2]-2)] #discard binary trt variables
    colnames(out)  = c(focalSpecies_s,focalSpecies_p)
    out[out<(1/popSize)]=0
    out_all_final  = out
    
    print('NO errors for the steady state integration, simulations continuing.')
  }else{
    print('An error occured for the steady state integration, terminating this set.')
  }
  
  
  #####################################################################################################
  if(!an.error.occured){
    print('made it here')
    # SS=out_all_final[dim(out_all_final)[1],2:dim(out_all_final)[2]]
    SS=out_all_final
    # Check linear stability of the equilibrium
    D         = diag(SS);
    S         = D*interaction_matrix;
    e         = eigen(S);
    isStable  = !(max(Re(e$values))>0);
    
    # Check if all abundances<1 and abundances>=0
    isReasonable  = !(any(SS<0) | any(SS>1));
    
    # Check if Sta1 and Stre1 coexist in the equilibrium
    isSta1Stre1 = ((SS$Sta1_s>eps) & (SS$Stre1_s>eps))
    
    # Combine
    is2use = (isStable & isReasonable & isSta1Stre1)
    
    if(is2use){
      if(plotOn_ss==1 & is2use==1){
        source('plot_convergence2SS.R')
      }
      print('... and the set is usable.')
      SS_all = rbind(SS_all,SS)
    }else{
      print(c(isStable,isReasonable,isSta1Stre1))
      print('...but not useable, terminating this set.')
      break
    }
  }else{
    is2use = FALSE
    break
  }
  
}