### functions to be sused in the estimation 
### of a joint OS/PFS distribution
require(mvtnorm)
CDF <- function(V,sigma){
  return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
}

## joint survival (non-vectorized) 
S_joi <- function(t1,t2, parms, obj_pfs, obj_os){
  # add here a later distinction   
  n_pfs <- length(obj_pfs$res[,"est"])
  n_os <- length(obj_os$res[,"est"])
  parms_pfs <- parms[1:n_pfs]
  parms_os <- parms[(n_pfs+1):(n_pfs + n_os)]
  
  S1_pfs <- do.call(obj_pfs$dfns$p, args = c(as.list(parms_pfs),
                                             list(q = t1, lower.tail = FALSE))
  )
  S1_os <- do.call(obj_os$dfns$p, args = c(as.list(parms_os),
                                           list(q = t2, lower.tail = FALSE))
  )
  z1 <- qnorm(1-S1_pfs)
  z2 <- qnorm(1-S1_os)
  sigma <- diag(2)
  sigma[1,2] <- sigma[2,1] <- parms[length(parms)]
  #if (is.na(z1)){
  #  print(paste0("parms  = ",parms))
  #  print(paste0("t1 = ",t1))
  #}
  #print(c(z1,z2))
  
  return(pmvnorm(lower = c(z1, z2), upper = c(Inf, Inf), sigma = sigma))
}

## example of S_joi(nt)

#S_joi_val <- S_joi(t1 = PFS[[1]][1], t2 = OS[[1]][1], parms = parms, obj_os = obj_os, obj_pfs = obj_pfs)


#S_joi_marPFS <- rep(0, length(t_pred))
#S_joi_marOS <- rep(0, length(t_pred))

#for (i in 1:length(t_pred)){
#  S_joi_marPFS[i] <- S_joi(t1 = t_pred[i], t2 = 0, parms = parms, obj_os = obj_os, obj_pfs = obj_pfs)
#  S_joi_marOS[i] <- S_joi(t1 = 0, t2 = t_pred[i], parms = parms, obj_os = obj_os, obj_pfs = obj_pfs)

#}
# test the marginals 
## to be added 

f_joi <- function(t1,t2, parms, obj_pfs, obj_os){
  n_pfs <- length(obj_pfs$res[,"est"])
  n_os <- length(obj_os$res[,"est"])
  parms_pfs <- parms[1:n_pfs]
  parms_os <- parms[(n_pfs+1):(n_pfs + n_os)]
  
  S1_pfs <- do.call(obj_pfs$dfns$p, args = c(as.list(parms_pfs),
                                             list(q = t1, lower.tail = FALSE))
  )
  S1_os <- do.call(obj_os$dfns$p, args = c(as.list(parms_os),
                                           list(q = t2, lower.tail = FALSE))
  )
  z1 <- qnorm(1-S1_pfs)
  z2 <- qnorm(1-S1_os)
  sigma <- diag(2)
  sigma[1,2] <- sigma[2,1] <- parms[length(parms)]
  num1 <- dmvnorm(c(z1, z2),  sigma =  sigma)
  den1 <- dnorm(z1)*dnorm(z2)
  
  f1_pfs <-  do.call(obj_pfs$dfns$d, args = c(as.list(parms_pfs),
                                              list (x = t1, log = FALSE))
  )
  
  f1_os <-  do.call(obj_os$dfns$d, args = c(as.list(parms_os),
                                            list (x = t2, log = FALSE))
  )
  num2 <- f1_os * f1_pfs
  return(num1*num2/den1)
}

f_joiA <- function(t1,t2, parms, obj_pfs, obj_os){
  
  n_pfs <- length(obj_pfs$res[,"est"])
  n_os <- length(obj_os$res[,"est"])
  parms_pfs <- parms[1:n_pfs]
  parms_os <- parms[(n_pfs+1):(n_pfs + n_os)]
  
  f1_pfs <-  do.call(obj_pfs$dfns$d, args = c(as.list(parms_pfs),
                                              list (x = t1, log = FALSE))
  )
  
  f1_os <-  do.call(obj_os$dfns$d, args = c(as.list(parms_os),
                                            list (x = t2, log = FALSE))
  )
  S1_pfs <- do.call(obj_pfs$dfns$p, args = c(as.list(parms_pfs),
                                             list(q = t1, lower.tail = FALSE))
  )
  S1_os <- do.call(obj_os$dfns$p, args = c(as.list(parms_os),
                                           list(q = t2, lower.tail = FALSE))
  )
  z1 <- qnorm(1-S1_pfs)
  z2 <- qnorm(1-S1_os)
  rho <- parms[length(parms)]
  
  g <- -0.5*log(1 - rho*rho) 
  g <- g - 0.5*(rho*rho*(z1*z1 + z2*z2) - 2*rho*z1*z2)/(1-rho*rho)
  return(f1_pfs*f1_os*exp(g))
  
}

#
#f_joi_val <- f_joi(t1 = PFS[[1]][1], t2 = OS[[1]][1], parms = parms, obj_os = obj_os, obj_pfs = obj_pfs)
#f_joi_valA <- f_joiA(t1 = PFS[[1]][1], t2 = OS[[1]][1], parms = parms, obj_os = obj_os, obj_pfs = obj_pfs)


## computing the terms with partial derivativesl
S1_p <- function(t1,t2, parms, obj_pfs, obj_os){
  n_pfs <- length(obj_pfs$res[,"est"])
  n_os <- length(obj_os$res[,"est"])
  parms_pfs <- parms[1:n_pfs]
  parms_os <- parms[(n_pfs+1):(n_pfs + n_os)]
  
  f1_pfs <-  do.call(obj_pfs$dfns$d, args = c(as.list(parms_pfs),
                                              list (x = t1, log = FALSE))
  )
  
  f1_os <-  do.call(obj_os$dfns$d, args = c(as.list(parms_os),
                                            list (x = t2, log = FALSE))
  )
  S1_pfs <- do.call(obj_pfs$dfns$p, args = c(as.list(parms_pfs),
                                             list(q = t1, lower.tail = FALSE))
  )
  S1_os <- do.call(obj_os$dfns$p, args = c(as.list(parms_os),
                                           list(q = t2, lower.tail = FALSE))
  )
  z1 <- qnorm(1-S1_pfs)
  z2 <- qnorm(1-S1_os)
  rho <- parms[length(parms)]
  
  Psi_pr <- dnorm(z1)* pnorm(z2, mean = rho*z1, sd = sqrt(1-rho*rho), lower.tail = FALSE)
  den1 <- dnorm(z1)
  ret_ <- Psi_pr * f1_pfs / den1 
  return(ret_)
}


## computing the terms with partial derivativesl
S2_p <- function(t1,t2, parms, obj_pfs, obj_os){
  n_pfs <- length(obj_pfs$res[,"est"])
  n_os <- length(obj_os$res[,"est"])
  parms_pfs <- parms[1:n_pfs]
  parms_os <- parms[(n_pfs+1):(n_pfs + n_os)]
  
  f1_pfs <-  do.call(obj_pfs$dfns$d, args = c(as.list(parms_pfs),
                                              list (x = t1, log = FALSE))
  )
  
  f1_os <-  do.call(obj_os$dfns$d, args = c(as.list(parms_os),
                                            list (x = t2, log = FALSE))
  )
  S1_pfs <- do.call(obj_pfs$dfns$p, args = c(as.list(parms_pfs),
                                             list(q = t1, lower.tail = FALSE))
  )
  S1_os <- do.call(obj_os$dfns$p, args = c(as.list(parms_os),
                                           list(q = t2, lower.tail = FALSE))
  )
  z1 <- qnorm(1-S1_pfs)
  z2 <- qnorm(1-S1_os)
  rho <- parms[length(parms)]
  
  Psi_pr <- dnorm(z2)* pnorm(z1, mean = rho*z2, sd = sqrt(1-rho*rho), lower.tail = FALSE)
  den1 <- dnorm(z2)
  ret_<- Psi_pr * f1_os / den1 
  return(ret_)
}


S_joi_vec <- function(t1,t2, parms, obj_pfs, obj_os){ ## t1 and t2 being vectors 
  # add here a later distinction   
  N <- length(t1)
  n_pfs <- length(obj_pfs$res[,"est"])
  n_os <- length(obj_os$res[,"est"])
  parms_pfs <- parms[1:n_pfs]
  parms_os <- parms[(n_pfs+1):(n_pfs + n_os)]
  
  S1_pfs <- do.call(obj_pfs$dfns$p, args = c(as.list(parms_pfs),
                                             list(q = t1, lower.tail = FALSE))
  )
  S1_os <- do.call(obj_os$dfns$p, args = c(as.list(parms_os),
                                           list(q = t2, lower.tail = FALSE))
  )
  z1 <- qnorm(1-S1_pfs)
  z2 <- qnorm(1-S1_os)

  #print(paste0(t2, " ", parms_os))
  
  
  sigma <- diag(2)
  sigma[1,2] <- sigma[2,1] <- parms[length(parms)]
  #if (is.na(z1)){
  #  print(paste0("parms  = ",parms))
  #  print(paste0("t1 = ",t1))
  #}
  #print(c(z1,z2))
  # create matrices for z1 and z2 
  # use of the biviariate normal for the computation 
  #joi1 <- pbivnorm(z1,z2,rho = parms[length(parms)])
  #joi_par1 <- pbivnorm(rep(1e120,N),z2,rho = parms[length(parms)])
  #joi_par2 <- pbivnorm(z1,rep(1e120,N),rho = parms[length(parms)])
  return(apply(cbind(z1,z2),1,CDF,sigma))
  #return(1 - joi_par1 - joi_par2 + joi1 )
}


likCop_vec <- function(parms, obj_pfs, obj_os, data_){
  N <- length(data_$USUBJID)
  logLik <- 0
  #print(parms)
  id00 <- which(data_$OSC == 0 & data_$PFSC == 0)
  id01 <- which(data_$OSC == 1 & data_$PFSC == 0)
  id10 <- which(data_$OSC == 0 & data_$PFSC == 1)
  id11 <- which(data_$OSC == 1 & data_$PFSC == 1)
  logLik <- logLik + sum(log(f_joiA(t1 = data_$PFS[id00], t2 =data_$OS[id00], parms = parms, obj_os = obj_os, obj_pfs = obj_pfs)))
  logLik <- logLik + sum(log(S1_p(t1 = data_$PFS[id01], t2 =data_$OS[id01], parms = parms, obj_os = obj_os, obj_pfs = obj_pfs)))
  logLik <- logLik + sum(log(S2_p(t1 = data_$PFS[id10], t2 =data_$OS[id10], parms = parms, obj_os = obj_os, obj_pfs = obj_pfs)))
  logLik <- logLik + sum(log(S_joi_vec(t1 = data_$PFS[id11], t2 =data_$OS[id11], parms = parms, obj_os = obj_os, obj_pfs = obj_pfs)))
  
  return(-logLik)
}


f_joiAres <- function(t1,t2, parms, obj_pfs, obj_os){
  
  n_pfs <- length(obj_pfs$res[,"est"])
  n_os <- length(obj_os$res[,"est"])
  parms_pfs <- parms[1:n_pfs]
  parms_os <- parms[(n_pfs+1):(n_pfs + n_os)]
  
  f1_pfs <-  do.call(obj_pfs$dfns$d, args = c(as.list(parms_pfs),
                                              list (x = t1, log = FALSE))
  )
  
  f1_os <-  do.call(obj_os$dfns$d, args = c(as.list(parms_os),
                                            list (x = t2, log = FALSE))
  )
  S1_pfs <- do.call(obj_pfs$dfns$p, args = c(as.list(parms_pfs),
                                             list(q = t1, lower.tail = FALSE))
  )
  S1_os <- do.call(obj_os$dfns$p, args = c(as.list(parms_os),
                                           list(q = t2, lower.tail = FALSE))
  )
  z1 <- qnorm(1-S1_pfs)
  z2 <- qnorm(1-S1_os)
  rho <- parms[length(parms)]
  
  g <- -0.5*log(1 - rho*rho) 
  g <- g - 0.5*(rho*rho*(z1*z1 + z2*z2) - 2*rho*z1*z2)/(1-rho*rho)
  
  if (t1 > t2) { 
    return(0)
  }else{
    return(f1_pfs*f1_os*exp(g))
  }
  
}


### addition of an inverse given a matrix representing the joint distribution 
margS_sample <- function(S_dir, time_point , index_ , N_sampl){
  
  # create a return object 
  ret_obj <-c()
  
  ret_ <- c()
  
  if (index_ == 1)
    S_use <- S_dir[time_point,]
  
  
  if (index_ == 2)
    S_use <- S_dir[,time_point]
  

    S_use  <- S_use/S_use[1]
    r_uni <- runif(N_sampl, 0, 1)
  
  for (i in 1:N_sampl){
    ret_ <- c(ret_, min( which(S_use < r_uni[i])))
  }
  ret_obj$mean <- mean(ret_[is.finite(ret_)])
  ret_obj$median <- median(ret_[is.finite(ret_)])
  ret_obj$ret_ <- ret_
  ret_obj$r_uni <- r_uni
  return(ret_obj)
} 