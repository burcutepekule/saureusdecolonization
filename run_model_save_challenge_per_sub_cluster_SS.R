inoculumSizes_eff = inoculumSizes_abs # this is already relative 
inoculumSizes     = inoculumSizes_eff
fuptimeIn         = fuptime

# YOU CAN ALWAYS GET THE FAMILY ORDER FROM 
# filename_params     = paste0('./BVS.results.parameters.txt')
# output_params       = read.csv(filename_params,sep="\t")
# output_growth       = output_params %>% filter(parameter_type=='growth_rate')
# family_order        = output_growth$target_taxon #keep same order for interaction matrix

# OR HARDCODE 
family_order   = c("Carnobacteriaceae","C1","C2","M1","M2","M3","M4","M5",
                   "P1","P2","P3","P4","P5","P6",
                   "Sta1","Sta2","Sta3","Sta4",
                   "Stre1","Stre2")

focalSpecies   = family_order
focalSpecies_p = paste0(family_order,'_p')
focalSpecies_s = paste0(family_order,'_s')
idxFocal       = match(focalSpecies,family_order)

mu_in          = growth_rates
alpha_in       = interaction_matrix
dt_in          = dt

params   = list(mu=mu_in,alpha=alpha_in,
                challengeTime=challengeTime,inoculumSizes=inoculumSizes,
                inoculumSpecs=inoculumSpecs,
                s2p=s2p,p2s=p2s,pd=pd,pm=pm, eps_in=eps_in, dt_in=dt_in)

for(i in 1:length(family_order)){ #first fix that
  if (ystart[i] < eps) ystart[i] = 0
}

times     = seq(from=0,to=fuptimeIn,by=dt) # returns a sequence
times_all = unique(round(sort(c(times,challengeTime)),round_time))
pseudoLog10 <- function(x) { asinh(x/2)/log(10) }


rootfun_challenge <- function (t, y, pars) {
  
  root1cInc  = any(y<0)
  y          = round(y,log10(eps_in))
  yNew       = log10(y[!is.infinite(log10(y))])
  root1aInc  = any(yNew<log10(eps))
  root1bInc  = any(is.na(log10(y)))
  
  if(root1aInc | root1bInc | root1cInc){
  # if(root1aInc){
    root1 = 0
  }else{
    root1 = 1
  }
  
  root2    = abs(t-challengeTime)
  
  rootBoth = c(root1,root2)
  
  # print(c(t,yNew))

  return(rootBoth)
}

eventfun_challenge_simple <- function(t, y, parms){ #first and second commensals 
  root2    = abs(t-challengeTime)
  
  if(any(root2<eps_time)){
    y[15] = y[15] + inoculumSizes
  }
  
  for(i in 1:length(y)){
    if (y[i] < eps) y[i] = 0
  }
  
  return(y)
}


out_challenge_pre=tryCatch({ode(
  func    = nasal.model.challenge,
  y       = ystart,
  times   = times_all,
  parms   = params,
  method  = solverMethod,
  rootfun = rootfun_challenge,
  events  = list(func = eventfun_challenge_simple, root=TRUE),
  atol = 1e-6,
  rtol = 1e-6
)}, error = function(e) {an.error.occured <<- TRUE})


if(!an.error.occured){
  out_challenge_pre = as.data.frame(out_challenge_pre)
  ystart_cut = unlist(out_challenge_pre[dim(out_challenge_pre)[1],2:dim(out_challenge_pre)[2]])
}

if(!an.error.occured){
  out=tryCatch({runsteady(
    func    = nasal.model.challenge,
    y       = ystart_cut,
    times   = c(0,Inf),
    parms   = params,
    stol = 1e-7,
    mf = 10 #If the problem is nonstiff, use method flag mf = 10, which selects a nonstiff (Adams) method, no Jacobian used..
  )}, error = function(e) {an.error.occured <<- TRUE})
}


if(!an.error.occured){
  
  out            = as.data.frame(t(as.data.frame(out)))
  colnames(out)  = c(focalSpecies_s,focalSpecies_p)
  out[out<(1/popSize)]=0
  out_all_final      = out
  out_all_challenge  = out_all_final
  print('NO errors for the steady state integration, simulations continuing.')
}else{
  print('An error occured for the steady state integration, terminating this set.')
}
