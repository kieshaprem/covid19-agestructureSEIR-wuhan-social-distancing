summariseSimulations = function(VAR,CI,SIMS)
{
  temp = lapply(SIMS,FUN = function(x) rowSums(x[[VAR]]))
  temp1 = do.call(cbind.data.frame, temp)
  var_p_median = apply(temp1,1,function(x) quantile(x,0.5))
  var_p_lci = apply(temp1,1,function(x) quantile(x,(1-CI/100)/2))
  var_p_uci = apply(temp1,1,function(x) quantile(x,1-(1-CI/100)/2))
  SUMMARY = list(median = var_p_median,lci = var_p_lci,uci=var_p_uci)
  rm(temp,temp1,var_p_median,var_p_lci,var_p_uci)
  
  
  results = list(summary=SUMMARY, Sim1 = SIMS[[1]])
  return(results)
  
}

summarisePeakTimePeakSize= function(VAR = 'incidence',SIMS)
{
  time = SIMS[[1]]$time
  temp = lapply(SIMS,FUN = function(x) rowSums(x[[VAR]]))
  temp1 = do.call(cbind.data.frame, temp)
  peaksize = as.numeric(apply(temp1,2,max))
  peaktime = time[as.numeric(apply(temp1,2,function(x) which.max(x)))]
  results = list(time = time, peaktime = peaktime,peaksize = peaksize)
  return(results)
}


summariseSimulations_mid = function(CI,SIMS)
{
  temp = lapply(SIMS,FUN = function(x) rowSums(x[['S']]))
  temp1 = do.call(cbind.data.frame, temp)
  i = which.min(abs(as.numeric(temp1[nrow(temp1),])-as.numeric(quantile(temp1[nrow(temp1),],0.5))))
  j = which.min(abs(as.numeric(temp1[nrow(temp1),])-as.numeric(quantile(temp1[nrow(temp1),],0.25))))
  k = which.min(abs(as.numeric(temp1[nrow(temp1),])-as.numeric(quantile(temp1[nrow(temp1),],0.75))))
  
  S = data.frame(med = temp[[i]],lci = temp[[k]],uci = temp[[j]],time = SIMS[[1]]$time)
  S_age = list(med = SIMS[[i]]$S,lci = SIMS[[k]]$S,uci = SIMS[[j]]$S)
  inc = list(med = SIMS[[i]]$incidence,lci = SIMS[[k]]$incidence,uci = SIMS[[j]]$incidence)
  time = SIMS[[1]]$time
  N_age = SIMS[[1]]$N_age
  results = list(S=S,inc = inc, time = time,N_age=N_age,S_age = S_age)
  return(results)
}

summariseSimulationsAGE = function(VAR,CI,SIMS)
{
  var_p_median = var_p_lci = var_p_uci = array(NA,c(length(SIMS[[1]]$time),16))
  for(age in 1:16)
  {
    temp = lapply(SIMS,FUN = function(x) (x[[VAR]][,age]))
    temp1 = do.call(cbind.data.frame, temp)
    var_p_median[,age] = apply(temp1,1,function(x) quantile(x,0.5))
    var_p_lci[,age] = apply(temp1,1,function(x) quantile(x,(1-CI/100)/2))
    var_p_uci[,age] = apply(temp1,1,function(x) quantile(x,1-(1-CI/100)/2))
    rm(temp,temp1)
  }
  
  SUMMARY = list(median = var_p_median,lci = var_p_lci,uci=var_p_uci)
  rm(var_p_median,var_p_lci,var_p_uci)
  
  
  results = list(summary=SUMMARY, Sim1 = SIMS[[1]])
  return(results)
  
}