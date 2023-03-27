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
  treatmentTimes     = c()
  for (k in 1:numCourses){
    treatmentTimes = c(treatmentTimes,k*(addOfBreakDays+addOfProTtreatmentDays)+seq(0,(addOfProTtreatmentDays-dt),dt))
  }
  treatmentTimes = treatmentTimes + proTreatmentStartIn + preTreatmentStartIn
  treatmentTimes = setdiff(treatmentTimes,preTreatmentTimes)
  source('run_model_save_antibiotic_per_sub_cluster.R')
  # }
}else{
  numCourses=0
  preTreatmentTimes = seq(0,(addOfPreTtreatmentDays-dt),dt) + preTreatmentStartIn
  treatmentTimes     = c()
  source('run_model_save_antibiotic_per_sub_cluster.R')
}
