
#setwd("D:\\LearnFiles\\Research paper\\idea\\MultiPoisFactor\\simu")
#setwd("/share/analysisdata/liuw/multicoap_res")
rm(list=ls())
source("helpfunc.R")
source('MSFR_main_R_MSFR_V1.R')
library(MultiCOAP)
library(dplyr)
### stepwise selection for q

selectQ <- function(res, tune.coef.rank=FALSE,  threshold=c(1e-1, 1e-2, 1e-2), 
                    method=c("var.prop", "SVR"), upper.var.prop=0.80){
  
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
    
  }
  return(qlist)
}
getq.from.mat <- function(mat, qmax){
  
  if(all(is.na(mat))) return(c(NA, NA, NA))
  mat[is.na(mat)] <- Inf
  id.min <- which.min(mat)
  hqs <- ceiling(id.min/qmax)
  hq <-  id.min%% qmax
  if(hq==0) hq <- qmax
  return(c(q=hq, qs1=hqs, qs2=hqs))
}
# Select the number of factors under fixed d ------------------------------

library(MultiCOAP)
p <- 100; nvec <- c(150,200); d <- 3; q <- 3
N <- 100
method_run <- "All"
nMethod <- 5 # ours, MSFR-log: AIC, BIC; MSFR: AIC, BIC
qArray <- array(NA, dim=c(3, nMethod, N))
row.names(qArray) <- c("q", "qs1", "qs2")
colnames(qArray) <- c("MultiCOAP", "AIC-log", "BIC-log", "AIC", "BIC")
timeMat <- matrix(NA, nrow=N, ncol=3)
colnames(timeMat) <- c("MultiCOAP", "MSFR-log", "MSFR")

for(i in 1:N){
  # i <- 1
  message("i = ", i)
  datList <- gendata_simu_multi(seed=i, nvec=nvec, p=p, q=q, d= d, rho=c(2,1), sigma2_eps=1)
  str(datList)
  XcList <- datList$Xlist; 
  ZList <- datList$Zlist; rank_use=d;
  q_max <- 6; qs_max <- 4
  
  aicMat1 <- array(NA, dim=c(q_max, qs_max))
  bicMat1 <- aicMat1
  tic <- proc.time()
  for(q1 in 1:q_max){
    message("q1 = ", q1, "/", q_max)
    for(qs1 in 1:qs_max){
      res_msfa <- MSFR_run(XcList, ZList, q=q1, qs=c(qs1, qs1), maxIter=1e3, log.transform=TRUE, load.source = T)
      aicMat1[q1, qs1] <- res_msfa$AIC
      bicMat1[q1, qs1] <- res_msfa$BIC
    }
  }
  toc <- proc.time()
  time.msfa.log <- toc[3] - tic[3]
  timeMat[i,2] <- time.msfa.log
  qArray[,2,i] <- getq.from.mat(aicMat1, qmax = q_max)
  qArray[,3,i] <- getq.from.mat(bicMat1, qmax = q_max)
  
  tic <- proc.time()
  aicMat2 <- array(NA, dim=c(q_max, qs_max))
  bicMat2 <- aicMat2
  for(q1 in 1:q_max){
    message("q1 = ", q1, "/", q_max)
    for(qs1 in 1:qs_max){
      try({
        res_msfa <- MSFR_run(XcList, ZList, q=q1, qs=c(qs1, qs1), maxIter=1e3, log.transform=FALSE, load.source = T)
        aicMat2[q1, qs1] <- res_msfa$AIC
        bicMat2[q1, qs1] <- res_msfa$BIC 
      }, silent=TRUE)
      
    }
  }
  toc <- proc.time()
  time.msfa <- toc[3] - tic[3]
  timeMat[i,3] <- time.msfa
  qArray[,4,i] <- getq.from.mat(aicMat2, qmax = q_max)
  qArray[,5,i] <- getq.from.mat(bicMat2, qmax = q_max)
  
  save(qArray, timeMat, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec),"q", q ,"p",p, "d", d,'_qArray.rds'))
  
}

save(qArray, timeMat, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec),"q", q ,"p",p, "d", d,'_qArray.rds'))



# ### Run MultiCOAP -----------------------------------------------------

N <- 100
method.use <- c("var.prop", "two.var.prop", "SVR", "twoSVR")
qArray <- array(NA, dim=c(3, length(method.use), N))
row.names(qArray) <- c("q", "qs1", "qs2")
colnames(qArray) <- method.use
timeMat <- matrix(NA, nrow=N, ncol=2)
colnames(timeMat) <- c("var.prop", "SVR")

for(i in 1:N){
  # i <- 2
  message("i = ", i)
  datList <- gendata_simu_multi(seed=i, nvec=nvec, p=p, q=q, d= d, rho=c(2,1), sigma2_eps=1)
  str(datList)
  XcList <- datList$Xlist; 
  ZList <- datList$Zlist; rank_use=d;
  q_max <- 6; qs_max <- 4
  res <- MultiCOAP(XcList, ZList, q=q_max, qs=rep(qs_max, length(XcList)),rank_use = d, verbose=F)
  # str(res)
  # trace_statistic_fun(res$A, datList$A0)
  #  trace_list_fun(res$B, datList$Blist0)
  ## select q first
  hq <- selectQ(res, method='var.prop')$q
  hq.svr <- selectQ(res, method='SVR')$q
  print(c(hq, hq.svr))
  res1 <- MultiCOAP(XcList, ZList, q=hq, qs=rep(qs_max, length(XcList)),rank_use = d, verbose=F)
  ### select qs next
  qArray[,1,i] <-   c(hq, unlist(selectQ(res1, method='var.prop'))[-1])
  qArray[,2,i] <-   unlist(selectQ(res1, method='var.prop'))
  
  res2 <- MultiCOAP(XcList, ZList, q=hq.svr, qs=rep(qs_max, length(XcList)),rank_use = d, verbose=F)
  qArray[,3,i] <-  c(hq.svr, unlist( selectQ(res2, method='SVR'))[-1])
  qArray[,4,i] <-  c(unlist( selectQ(res2, method='SVR')))
}
apply(qArray, c(1,2), mean, na.rm=T)
apply(qArray, c(1,2), sd, na.rm=T)
method_run <- "MultiCOAP4"
save(qArray, timeMat, file=paste0(dir.name,"_",method_run,"nsum", sum(nvec),"q", q ,"p",p, "d", d,'_qArray.rds'))




