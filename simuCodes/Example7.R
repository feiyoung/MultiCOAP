## d is high-dimensional
rm(list=ls())
source("helpfunc.R")
source('MSFR_main_R_MSFR_V1.R')
library(MultiCOAP)
library(dplyr)

# Estimation performance of parameters: change p from (20, 100, 250) ------------------------------------
pvec <- c(20, 100, 250)
i <- commandArgs(TRUE) %>% as.integer()
#i <- 1
p <- pvec[i] #250;
nvec <- c(150,200); d <- 50; q <- 3; rank0 <- 3; qs <- c(2,2)
N <- 100
method_run <- 'All'#"MultiCOAP_COAP"
methodNames <- c("MultiCOAP", "COAP", "MSFRlog", "MSFR")
n_methods <- length(methodNames)
metricList <- list(A_tr=matrix(NA,N, n_methods), B_tr=matrix(NA, N, n_methods),
                   beta_er=matrix(NA, N, n_methods),
                   labmda_er = matrix(NA, N, n_methods),
                   F_tr =matrix(NA,N, n_methods),H_tr =matrix(NA,N, n_methods), 
                   timeMat = matrix(NA, N, n_methods))
for(i in 1:N){
  # i <- 1
  
  message("i = ", i)
  datList <- gendata_simu_multi1(seed=i, nvec=nvec, p=p, q=q, d= d, 
                                rank0=rank0,rho=c(1,1,1), sigma2_eps=1)
  str(datList)
  XcList <- datList$Xlist; 
  ZList <- datList$Zlist; rank_use=d;
  
  
  res <- MultiCOAP(XcList, ZList, q=q, qs=qs,rank_use = rank0, verbose=F)
  metricList$timeMat[i,1] <- res$time.use
  metricList$A_tr[i,1] <- trace_statistic_fun(res$A, datList$A0)
  metricList$B_tr[i,1] <- trace_list_fun(res$B, datList$Blist0)
  metricList$F_tr[i,1] <- trace_list_fun(res$F, datList$Flist)
  metricList$H_tr[i,1] <- trace_list_fun(res$H, datList$Hlist)
  metricList$beta_er[i, 1] <- normvec(datList$bbeta0-res$bbeta)
  metricList$labmda_er[i, 1] <- normvec(datList$lambdavec -res$Lambda)



  ## COAP
  res_coap <- COAP_run(XcList, ZList, q=q,rank_use = rank0, verbose=F)
  str(res_coap)
  metricList$timeMat[i,2] <- res_coap$time.use
  metricList$A_tr[i,2] <- trace_statistic_fun(res_coap$B, datList$A0)
  metricList$F_tr[i,2] <- trace_list_fun(res_coap$Flist, datList$Flist)
  metricList$beta_er[i, 2] <- normvec(datList$bbeta0-res_coap$bbeta)
  
  
  # MSFR-log
  res_msfa <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=TRUE, load.source = T)
  metricList$timeMat[i,3] <- res_msfa$time.use
  metricList$A_tr[i,3] <- trace_statistic_fun(res_msfa$Phi, datList$A0)
  metricList$B_tr[i,3] <- trace_list_fun(res_msfa$Lambda_s, datList$Blist0)
  metricList$F_tr[i,3] <- trace_list_fun(res_msfa$Flist, datList$Flist)
  metricList$H_tr[i,3] <- trace_list_fun(res_msfa$Hlist, datList$Hlist)
  metricList$beta_er[i, 3] <- normvec(res_msfa$beta-datList$bbeta0)
  metricList$labmda_er[i, 3] <- normvec(datList$lambdavec -res_msfa$lambdavec)

  # MSFR
  try({
    res_msfa2 <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=F, load.source = T)
    metricList$timeMat[i,4] <- res_msfa2$time.use
    metricList$A_tr[i,4] <- trace_statistic_fun(res_msfa2$Phi, datList$A0)
    metricList$B_tr[i,4] <- trace_list_fun(res_msfa2$Lambda_s, datList$Blist0)
    metricList$F_tr[i,4] <- trace_list_fun(res_msfa2$Flist, datList$Flist)
    metricList$H_tr[i,4] <- trace_list_fun(res_msfa2$Hlist, datList$Hlist)
    metricList$beta_er[i, 4] <- normvec(res_msfa2$beta-datList$bbeta0)
    metricList$labmda_er[i, 4] <- normvec(datList$lambdavec -res_msfa2$lambdavec)
  },silent=T)

  save(metricList, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec), "p",p, "d", d,'_metricList.rds'))

}
str(metricList)
save(metricList, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec), "p",p, "d", d,'_metricList.rds'))

head(metricList$A_tr)
head(metricList$F_tr)
head(metricList$beta_er)



# Add p = 200 case --------------------------------------------------------
p <- 200
  nvec <- c(150,200); d <- 50; q <- 3; rank0 <- 3; qs <- c(2,2)
  N <- 100
  method_run <- 'All'#"MultiCOAP_COAP"
  methodNames <- c("MultiCOAP", "COAP", "MSFRlog", "MSFR")
  n_methods <- length(methodNames)
  metricList <- list(A_tr=matrix(NA,N, n_methods), B_tr=matrix(NA, N, n_methods),
                     beta_er=matrix(NA, N, n_methods),
                     labmda_er = matrix(NA, N, n_methods),
                     F_tr =matrix(NA,N, n_methods),H_tr =matrix(NA,N, n_methods), 
                     timeMat = matrix(NA, N, n_methods))
  for(i in 1:N){
    # i <- 1
    
    message("i = ", i)
    datList <- gendata_simu_multi1(seed=i, nvec=nvec, p=p, q=q, d= d, 
                                   rank0=rank0,rho=c(1,1,1), sigma2_eps=1)
    str(datList)
    XcList <- datList$Xlist; 
    ZList <- datList$Zlist; rank_use=d;
    
    
    res <- MultiCOAP(XcList, ZList, q=q, qs=qs,rank_use = rank0, verbose=F, init="MSFRVI")
    metricList$timeMat[i,1] <- res$time.use
    metricList$A_tr[i,1] <- trace_statistic_fun(res$A, datList$A0)
    metricList$B_tr[i,1] <- trace_list_fun(res$B, datList$Blist0)
    metricList$F_tr[i,1] <- trace_list_fun(res$F, datList$Flist)
    metricList$H_tr[i,1] <- trace_list_fun(res$H, datList$Hlist)
    metricList$beta_er[i, 1] <- normvec(datList$bbeta0-res$bbeta)
    metricList$labmda_er[i, 1] <- normvec(datList$lambdavec -res$Lambda)
    
    
    
    ## COAP
    res_coap <- COAP_run(XcList, ZList, q=q,rank_use = rank0, verbose=F)
    str(res_coap)
    metricList$timeMat[i,2] <- res_coap$time.use
    metricList$A_tr[i,2] <- trace_statistic_fun(res_coap$B, datList$A0)
    metricList$F_tr[i,2] <- trace_list_fun(res_coap$Flist, datList$Flist)
    metricList$beta_er[i, 2] <- normvec(datList$bbeta0-res_coap$bbeta)
    
    
    # MSFR-log
    res_msfa <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=TRUE, load.source = T)
    metricList$timeMat[i,3] <- res_msfa$time.use
    metricList$A_tr[i,3] <- trace_statistic_fun(res_msfa$Phi, datList$A0)
    metricList$B_tr[i,3] <- trace_list_fun(res_msfa$Lambda_s, datList$Blist0)
    metricList$F_tr[i,3] <- trace_list_fun(res_msfa$Flist, datList$Flist)
    metricList$H_tr[i,3] <- trace_list_fun(res_msfa$Hlist, datList$Hlist)
    metricList$beta_er[i, 3] <- normvec(res_msfa$beta-datList$bbeta0)
    metricList$labmda_er[i, 3] <- normvec(datList$lambdavec -res_msfa$lambdavec)
  
    # MSFR
    try({
      res_msfa2 <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=F, load.source = T)
      metricList$timeMat[i,4] <- res_msfa2$time.use
      metricList$A_tr[i,4] <- trace_statistic_fun(res_msfa2$Phi, datList$A0)
      metricList$B_tr[i,4] <- trace_list_fun(res_msfa2$Lambda_s, datList$Blist0)
      metricList$F_tr[i,4] <- trace_list_fun(res_msfa2$Flist, datList$Flist)
      metricList$H_tr[i,4] <- trace_list_fun(res_msfa2$Hlist, datList$Hlist)
      metricList$beta_er[i, 4] <- normvec(res_msfa2$beta-datList$bbeta0)
      metricList$labmda_er[i, 4] <- normvec(datList$lambdavec -res_msfa2$lambdavec)
    },silent=T)
  
    save(metricList, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec), "p",p, "d", d,'_metricList.rds'))
  
  }
  str(metricList)



# Estimation performance of parameters: p=50, change n from (100,150) ------------------------------------
p <- 50
nvec.list <- list(c(100, 150), c(150,200), c(200, 300))
i <- commandArgs(TRUE) %>% as.integer()
#i <- 3

nvec <- nvec.list[[i]]
d <- 50; q <- 3; rank0 <- 3; qs <- c(2,2)
N <- 100
method_run <- 'All'#"MultiCOAP_COAP"
methodNames <- c("MultiCOAP", "COAP", "MSFRlog", "MSFR")
n_methods <- length(methodNames)
metricList <- list(A_tr=matrix(NA,N, n_methods), B_tr=matrix(NA, N, n_methods),
                   beta_er=matrix(NA, N, n_methods),
                   labmda_er = matrix(NA, N, n_methods),
                   F_tr =matrix(NA,N, n_methods),H_tr =matrix(NA,N, n_methods), 
                   timeMat = matrix(NA, N, n_methods))
for(i in 1:N){
  # i <- 1
  
  message("i = ", i)
  datList <- gendata_simu_multi1(seed=i, nvec=nvec, p=p, q=q, d= d, 
                                 rank0=rank0,rho=c(1,1,1), sigma2_eps=1)
  str(datList)
  XcList <- datList$Xlist; 
  ZList <- datList$Zlist; rank_use=d;
  
  
  res <- MultiCOAP(XcList, ZList, q=q, qs=qs,rank_use = rank0, verbose=F)
  metricList$timeMat[i,1] <- res$time.use
  metricList$A_tr[i,1] <- trace_statistic_fun(res$A, datList$A0)
  metricList$B_tr[i,1] <- trace_list_fun(res$B, datList$Blist0)
  metricList$F_tr[i,1] <- trace_list_fun(res$F, datList$Flist)
  metricList$H_tr[i,1] <- trace_list_fun(res$H, datList$Hlist)
  metricList$beta_er[i, 1] <- normvec(datList$bbeta0-res$bbeta)
  metricList$labmda_er[i, 1] <- normvec(datList$lambdavec -res$Lambda)
  
  
  
  ## COAP
  res_coap <- COAP_run(XcList, ZList, q=q,rank_use = rank0, verbose=F)
  str(res_coap)
  metricList$timeMat[i,2] <- res_coap$time.use
  metricList$A_tr[i,2] <- trace_statistic_fun(res_coap$B, datList$A0)
  metricList$F_tr[i,2] <- trace_list_fun(res_coap$Flist, datList$Flist)
  metricList$beta_er[i, 2] <- normvec(datList$bbeta0-res_coap$bbeta)
  
  
  # # # MSFR-log
  res_msfa <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=TRUE, load.source = T)
  metricList$timeMat[i,3] <- res_msfa$time.use
  metricList$A_tr[i,3] <- trace_statistic_fun(res_msfa$Phi, datList$A0)
  metricList$B_tr[i,3] <- trace_list_fun(res_msfa$Lambda_s, datList$Blist0)
  metricList$F_tr[i,3] <- trace_list_fun(res_msfa$Flist, datList$Flist)
  metricList$H_tr[i,3] <- trace_list_fun(res_msfa$Hlist, datList$Hlist)
  metricList$beta_er[i, 3] <- normvec(res_msfa$beta-datList$bbeta0)
  metricList$labmda_er[i, 3] <- normvec(datList$lambdavec -res_msfa$lambdavec)
  # 
  # MSFR
  try({
    res_msfa2 <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=F)
    metricList$timeMat[i,4] <- res_msfa2$time.use
    metricList$A_tr[i,4] <- trace_statistic_fun(res_msfa2$Phi, datList$A0)
    metricList$B_tr[i,4] <- trace_list_fun(res_msfa2$Lambda_s, datList$Blist0)
    metricList$F_tr[i,4] <- trace_list_fun(res_msfa2$Flist, datList$Flist)
    metricList$H_tr[i,4] <- trace_list_fun(res_msfa2$Hlist, datList$Hlist)
    metricList$beta_er[i, 4] <- normvec(res_msfa2$beta-datList$bbeta0)
    metricList$labmda_er[i, 4] <- normvec(datList$lambdavec -res_msfa2$lambdavec)
  },silent=T)
  
  save(metricList, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec), "p",p, "d", d,'_metricList.rds'))
  
}
str(metricList)
save(metricList, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec), "p",p, "d", d,'_metricList.rds'))

head(metricList$A_tr)
head(metricList$F_tr)
head(metricList$beta_er)


# Select the number of factors and ranks ----------------------------------
selectQ.COAP <- function(reslist,threshold=c(1e-1, 1e-2)){
  
  thre1 <- threshold[1]
  beta_svalues <- svd(reslist$bbeta)$d
  beta_svalues <- beta_svalues[beta_svalues>thre1]
  ratio1 <- beta_svalues[-length(beta_svalues)] / beta_svalues[-1]
  hr <- which.max(ratio1[-length(ratio1)])
  
  
  thre2 <- threshold[2]
  B_svalues <- svd(reslist$B)$d
  B_svalues <- B_svalues[B_svalues>thre2]
  ratio_fac <- B_svalues[-length(B_svalues)] / B_svalues[-1]
  hq <- which.max(ratio_fac)
  
  return(c(hr=hr, hq=hq))
}

selectQR <- function(res, tune.coef.rank=FALSE,method=c("SVR", "var.prop"), 
                     threshold=c(1e-1, 1e-2, 1e-2), upper.var.prop=0.80){
  
  method <- match.arg(method)
  qlist <- list()
  
  d_svdA <- svd(res$A)$d
  # 
  if(method=='SVR'){
    d_svdA <- d_svdA[d_svdA>threshold[1]]
    qq <- length(d_svdA)
    if(qq>1){
      rsigs <- d_svdA[-qq]/d_svdA[-1]
      qlist$q <-  which.max(rsigs)
    }else{
      qlist$q <- 1
    }
  }else if(method=='var.prop'){
    prop.vec <- cumsum(d_svdA^2)/sum(d_svdA^2)
    qlist$q <- which(prop.vec> upper.var.prop)[1]
  }
  
  n_qs <- length(res$B)
  qvec <- rep(NA, n_qs)
  names(qvec) <- paste0("qs", 1:n_qs)
  for(i in 1:n_qs){
    # i <-1 
    d_svdB1 <- svd(res$B[[i]])$d
    # 
    if(method=='SVR'){
      d_svdB1 <- d_svdB1[d_svdB1>threshold[2]]
      qq1 <- length(d_svdB1)
      if(qq1>1){
        rsigB1 <- d_svdB1[-qq1]/d_svdB1[-1]
        qvec[i] <-  which.max(rsigB1)
      }else{
        qvec[i] <-  1
      }
    }else if(method=='var.prop'){
      prop.vec <- cumsum(d_svdB1^2)/sum(d_svdB1^2)
      qvec[i] <- which(prop.vec> upper.var.prop)[1]
    }
  }
  qlist$qs <- qvec
  if(tune.coef.rank){
    d_svdA <- svd(res$bbeta)$d[1:res$qrlist$r]
    if(method=='SVR'){
    d_svdA <- d_svdA[d_svdA>threshold[3]]
    qq <- length(d_svdA)
    if(qq>2){
      d_svdA <- d_svdA[-1] # remove the intercept term.
      rsigs <- d_svdA[-(qq-1)]/d_svdA[-1]
      qlist$r <-  which.max(rsigs) 
    }else{
      qlist$r <- 1
    }
    }else if(method=='var.prop'){
      prop.vec <- cumsum(d_svdA^2)/sum(d_svdA^2)
      qlist$r <- which(prop.vec> upper.var.prop)[1]
    }
  }
  return(qlist)
}

p <- 100
nvec <- c(150,200); d <- 50; q <- 3; rank0 <- 3; qs <- c(2,2)
N <- 100
method.use <- c("var.prop", "two.var.prop", "SVR", "twoSVR")
qArray <- array(NA, dim=c(4, length(method.use), N))
row.names(qArray) <- c("q", "qs1", "qs2", "r")
colnames(qArray) <- method.use
sigma2_eps <- 2
for(i in 1:N){
  # i<- 1
  message("i = ", i)
  datList <- gendata_simu_multi1(seed=i, nvec=nvec, p=p, q=q, d= d, 
                                 rank0=rank0,rho=c(2,3,1), sigma2_eps=sigma2_eps)
  str(datList)
  XcList <- datList$Xlist; 
  ZList <- datList$Zlist; rank_use=d;
  q_max <- 6; qs_max <- 4; r_max <- 10
  res <- MultiCOAP(XcList, ZList, q= q_max, qs=rep(qs_max,2), rank_use = r_max, verbose = F)
  qlist <- selectQR(res, tune.coef.rank = TRUE, method='var.prop')
  selectQR(res, tune.coef.rank = TRUE, method='SVR')
  res_tmp <- MultiCOAP(XcList, ZList, q=qlist$q, qs = qlist$qs ,rank_use = r_max, verbose=F)
  selectQR(res_tmp, tune.coef.rank = TRUE, method='var.prop')
  qArray[,1,i] <- unlist(qlist)
  res2 <- MultiCOAP(XcList, ZList, q=qlist$q, qs = qlist$qs ,rank_use = qlist$r, verbose=F)
  qlist2 <- selectQR(res2, tune.coef.rank = T, method='var.prop')
  qArray[,2,i] <- c(qlist$q, qlist2$qs, qlist2$r)
  qArray[,3,i] <- unlist(selectQR(res, tune.coef.rank = T, method='SVR'))
  qlist.svr <- selectQR(res, tune.coef.rank = T, method='SVR')
  res3 <- MultiCOAP(XcList, ZList, q=qlist.svr$q, qs = rep(qs_max,2),rank_use = qlist.svr$r)
  qArray[,4,i] <- c(unlist(qlist.svr)[1],
                    unlist(selectQR(res3, tune.coef.rank = T, method='SVR'))[-1])
  
  save(qArray, file=paste0(dir.name,"_","nsum", sum(nvec), "p",p, "d", d,'_selectRank_sigma2_',sigma2_eps,'.rds'))
}
apply(qArray, c(1,2), mean)
apply(qArray, c(1,2), sd)

# res_coap <- COAP_run(XcList, ZList, q=q_max,rank_use = r_max)
# selectQ.COAP(res_coap)[2:1]



# The influence of misselection of r --------------------------------------
## we fix q and qs
p <- 100
nvec <- c(150,200); d <- 50; q <- 3; rank0 <- 3; qs <- c(2,2)
N <- 100
r.vec <- c(1,3,5, 9)
ii <- commandArgs(TRUE) %>% as.integer()
# ii <- 1
N <- 100
method_run <- 'All'#"MultiCOAP_COAP"
methodNames <- c("MultiCOAP", "COAP")
n_methods <- length(methodNames)
metricList <- list(A_tr=matrix(NA,N, n_methods), B_tr=matrix(NA, N, n_methods),
                   beta_er=matrix(NA, N, n_methods),
                   labmda_er = matrix(NA, N, n_methods),
                   F_tr =matrix(NA,N, n_methods),H_tr =matrix(NA,N, n_methods), 
                   timeMat = matrix(NA, N, n_methods))
for(i in 1:N){
  # i <- 2
  
  message("i = ", i)
  datList <- gendata_simu_multi1(seed=i, nvec=nvec, p=p, q=q, d= d, 
                                 rank0=rank0,rho=c(2,3,1), sigma2_eps=1)
  str(datList)
  XcList <- datList$Xlist; 
  ZList <- datList$Zlist; rank_use= r.vec[ii];
  
  
  res <- MultiCOAP(XcList, ZList, q=q, qs=qs,rank_use = rank_use, verbose=F, init='MSFRVI')
  metricList$timeMat[i,1] <- res$time.use
  metricList$A_tr[i,1] <- trace_statistic_fun(res$A, datList$A0)
  metricList$B_tr[i,1] <- trace_list_fun(res$B, datList$Blist0)
  metricList$F_tr[i,1] <- trace_list_fun(res$F, datList$Flist)
  metricList$H_tr[i,1] <- trace_list_fun(res$H, datList$Hlist)
  metricList$beta_er[i, 1] <- normvec(datList$bbeta0-res$bbeta)
  metricList$labmda_er[i, 1] <- normvec(datList$lambdavec -res$Lambda)
  
  
  
  ## COAP
  res_coap <- COAP_run(XcList, ZList, q=q,rank_use = rank_use, verbose=F)
  str(res_coap)
  metricList$timeMat[i,2] <- res_coap$time.use
  metricList$A_tr[i,2] <- trace_statistic_fun(res_coap$B, datList$A0)
  metricList$F_tr[i,2] <- trace_list_fun(res_coap$Flist, datList$Flist)
  metricList$beta_er[i, 2] <- normvec(datList$bbeta0-res_coap$bbeta)
  
  
  save(metricList, file=paste0(dir.name,"_misselect_r_nsum", sum(nvec), "p",p, "d", d,"_rankuse",rank_use,'_metricList.rds'))
  
}
metricList


