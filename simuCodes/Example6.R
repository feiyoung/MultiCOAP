#### Compare the computational time by change nvec.
rm(list=ls())
setwd("/home/liuwei1/Rfiles/OtherDataAnalysis/MultiCOAP")
source("helpfunc.R")
source('MSFR_main_R_MSFR_V1.R')
library(MultiCOAP)
library(dplyr)


# Fix p=800, change nvec ---------------------------------------------------------------
p <- 800
nvec.base <- c(200, 300)
cseq <- seq(1, 20, by=3)

d <- 3; q <- 3; rank0 <- 3; qs <- c(2,2)
N <- 5
method_run <- 'All'#"MultiCOAP_COAP"
methodNames <- c("MultiCOAP", "COAP", "MSFRlog", "MSFR")
n_methods <- length(methodNames)
timeArray = array(NA, dim=c(N, n_methods, length(cseq)))
colnames(timeArray) <- methodNames
for(jj in seq_along(cseq)){
  ## jj <-1 
  for(i in 1:N){
    # i <- 1
    
    message("jj= ", jj,", i = ", i)
    nvec <- nvec.base * cseq[jj]
    datList <- gendata_simu_multi1(seed=i, nvec=nvec, p=p, q=q, d= d, 
                                   rank0=rank0,rho=c(1,1,1), sigma2_eps=1)
    str(datList)
    XcList <- datList$Xlist; 
    ZList <- datList$Zlist; rank_use=d;
    
    
    res <- MultiCOAP(XcList, ZList, q=q, qs=qs,rank_use = rank0, verbose=T, epsELBO=1e-10)
    timeArray[i,1,jj] <- res$time.use
    
    
    
    ## COAP
    res_coap <- COAP_run(XcList, ZList, q=q,rank_use = rank0, verbose=T, epsELBO=1e-10)
    timeArray[i,2,jj] <- res_coap$time.use
    
    # # # MSFR-log
    try({
    res_msfa <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=TRUE, load.source = T)
    timeArray[i,3,jj] <- res_msfa$time.use
    },silent=T)
    # 
    # MSFR
    try({
      res_msfa2 <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=F, load.source=T)
      timeArray[i,4,jj] <- res_msfa2$time.use
    },silent=T)
    
    
  }
  save(timeArray, file=paste0(dir.name,"_",method_run,"nseq_",  "p",p, "d", d,'_timeArray.rds'))
  
}









# #Fix nvec=c(1000, 2000); change p ---------------------------------------


nvec <- c(1000, 2000)
p.seq <- c(200, 500, 800,  900,  2000, 5000, 8000, 1e4)

d <- 3; q <- 3; rank0 <- 3; qs <- c(2,2)
N <- 5
method_run <- 'All'#"MultiCOAP_COAP"
methodNames <- c("MultiCOAP", "COAP", "MSFRlog", "MSFR")
n_methods <- length(methodNames)
timeArray = array(NA, dim=c(N, n_methods, length(p.seq)))
colnames(timeArray) <- methodNames
for(jj in 5:length(p.seq)){
  ## jj <-1 
  for(i in 1:N){
    # i <- 1
    p <- p.seq[jj]
    message("jj= ", jj,", i = ", i)
    datList <- gendata_simu_multi1(seed=i, nvec=nvec, p=p, q=q, d= d, 
                                   rank0=rank0,rho=c(1,1,1), sigma2_eps=1)
    str(datList)
    XcList <- datList$Xlist; 
    ZList <- datList$Zlist; rank_use=d;
    
    
    res <- MultiCOAP(XcList, ZList, q=q, qs=qs,rank_use = rank0, verbose=T, epsELBO=1e-10)
    timeArray[i,1,jj] <- res$time.use
    
    
    
    ## COAP
    res_coap <- COAP_run(XcList, ZList, q=q,rank_use = rank0, verbose=T, epsELBO=1e-10)
    timeArray[i,2,jj] <- res_coap$time.use
    
    # # # MSFR-log
    try({
    res_msfa <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=TRUE, load.source = T)
    timeArray[i,3,jj] <- res_msfa$time.use
    },silent=T)
    # 
    # MSFR
    try({
      res_msfa2 <- MSFR_run(XcList, ZList, q=q, qs=qs, maxIter=1e3, log.transform=F, load.source=T)
      timeArray[i,4,jj] <- res_msfa2$time.use
    },silent=T)
    
    
  }
  save(timeArray, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec), "pseq_","d", d,'_timeArray.rds'))
  
}


