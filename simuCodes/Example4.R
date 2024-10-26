
### Fix p=100, nvec <- c(150,200); d <- 3; q <- 3
### Change (rho_A, rho_B) \in {(1,1), (1,3), (2,1)}
rm(list=ls())

dir.old <- "/home/liuwei1/Rfiles/OtherDataAnalysis/MultiCOAP"
setwd(dir.old)
source("helpfunc.R")
source('MSFR_main_R_MSFR_V1.R')
library(MultiCOAP)
library(dplyr)

# Estimation performance of parameters: Fix rhoA=1, change rhoB ------------------------------------
## Fix rhoA=1, change rhoB
rho.ABlist <- list(c(0.8,1),  c(2,1), c(2,3))
## (2,1) (2,5), the results are good!
# i <- commandArgs(TRUE) %>% as.integer()
i <- 1  
rhoAB <- rho.ABlist[[i]]
p <- 100
nvec <- c(150,200); d <- 3; q <- 3; rank0 <- 3; qs <- c(2,2)
N <- 5
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
                                 rank0=rank0,rho=c(rhoAB[1],rhoAB[2],1), sigma2_eps=1)
  # str(datList)
  XcList <- datList$Xlist; 
  ZList <- datList$Zlist; rank_use=d;
  
  
  res <- MultiCOAP(XcList, ZList, q=q, qs=qs,rank_use = rank0, verbose=F,init="LFM") # "MSFRVI"
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
  str(metricList)
  
  # MSFR-log
  res_msfa <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=TRUE, load.source=TRUE)
  str(res_msfa)
  res_msfa$Flist <- lapply(res_msfa$E_f, t)
  res_msfa$Hlist <- lapply(res_msfa$E_l, t)
  metricList$timeMat[i,3] <- res_msfa$time.use
  metricList$A_tr[i,3] <- trace_statistic_fun(res_msfa$Phi, datList$A0)
  metricList$B_tr[i,3] <- trace_list_fun(res_msfa$Lambda_s, datList$Blist0)
  metricList$F_tr[i,3] <- trace_list_fun(res_msfa$Flist, datList$Flist)
  metricList$H_tr[i,3] <- trace_list_fun(res_msfa$Hlist, datList$Hlist)
  metricList$beta_er[i, 3] <- normvec(res_msfa$beta-datList$bbeta0)
  metricList$labmda_er[i, 3] <- normvec(datList$lambdavec -res_msfa$lambdavec)

  # MSFR
  try({
    res_msfa2 <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=F, load.source=TRUE)
    metricList$timeMat[i,4] <- res_msfa2$time.use
    metricList$A_tr[i,4] <- trace_statistic_fun(res_msfa2$Phi, datList$A0)
    metricList$B_tr[i,4] <- trace_list_fun(res_msfa2$Lambda_s, datList$Blist0)
    metricList$F_tr[i,4] <- trace_list_fun(res_msfa2$Flist, datList$Flist)
    metricList$H_tr[i,4] <- trace_list_fun(res_msfa2$Hlist, datList$Hlist)
    metricList$beta_er[i, 4] <- normvec(res_msfa2$beta-datList$bbeta0)
    metricList$labmda_er[i, 4] <- normvec(datList$lambdavec -res_msfa2$lambdavec)
  },silent=T)
  save(metricList, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec),"rhoAB",paste0(rhoAB, collapse = "_"),
                               "q", q ,"p",p, "d", d,'_metricList.rds'))
  
}
str(metricList)
save(metricList, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec),"rhoAB",paste0(rhoAB, collapse = "_"),
                             "q", q ,"p",p, "d", d,'_metricList.rds'))


head(metricList$A_tr)
head(metricList$F_tr)
head(metricList$beta_er)
