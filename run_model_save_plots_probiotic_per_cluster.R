#convert days to months
preTreatmentStartIn    = preTreatmentStart*dt 
proTreatmentStartIn    = proTreatmentStart*dt 
addOfPreTtreatmentDays = numOfPreTtreatmentDays*dt;
addOfProTtreatmentDays = numOfProTtreatmentDays*dt;
addOfBreakDays         = numOfBreakDays*dt;
fuptimeIn              = fuptime*dt


if(numOfTotalCourses>0){
  numCourses = numOfTotalCourses
  preTreatmentTimes = seq(0,(addOfPreTtreatmentDays-dt),dt) + preTreatmentStartIn
  inoculumTimes     = c()
  for (k in 1:numCourses){
    inoculumTimes = c(inoculumTimes,k*(addOfBreakDays+addOfProTtreatmentDays)+seq(0,(addOfProTtreatmentDays-dt),dt))
  }
  inoculumTimes = inoculumTimes + proTreatmentStartIn + preTreatmentStartIn
  inoculumTimes = setdiff(inoculumTimes,preTreatmentTimes)
  source('run_model_save_probiotic_per_sub_cluster.R')
}else{
  numCourses = 0
  preTreatmentTimes = seq(0,(addOfPreTtreatmentDays-dt),dt) + preTreatmentStartIn
  inoculumTimes     = c()
  source('run_model_save_probiotic_per_sub_cluster.R')
}
