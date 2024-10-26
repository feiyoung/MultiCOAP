


# Generat data ------------------------------------------------------------
# 
# gendata_simu_multi1 <- function (seed = 1, nvec = c(100,300), p = 50, d=3, q = 3,
#                                  qs= rep(2, length(nvec)),
#                                  rank0=3, rho = c(rhoA=1, rhoB= 1, rhoZ=1), sigma2_eps=1){
#   # seed = 1; nvec = c(100,300); p = 50; d=3; q = 3
#   # qs= rep(2, length(nvec))
#   # rank0=3; rho = c(1, 1); sigma2_eps=0.1;
#   
#   
#   if(length(nvec)<2) stop("nvec must have at least two elements!")
#   S <- length(nvec)
#   
#   require(MASS)
#   if(rank0<=1) stop("rank0 must be greater than 1!")
#   cor.mat<-function (p, rho, type = "toeplitz") {
#     if (p == 1) 
#       return(matrix(1, 1, 1))
#     mat <- diag(p)
#     if (type == "toeplitz") {
#       for (i in 2:p) {
#         for (j in 1:i) {
#           mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
#         }
#       }
#     }
#     if (type == "identity") {
#       mat[mat == 0] <- rho
#     }
#     return(mat)
#   }
#   Diag<-function (vec){
#     q <- length(vec)
#     if (q > 1) {
#       y <- diag(vec)
#     }
#     else {
#       y <- matrix(vec, 1, 1)
#     }
#     return(y)
#   }
#   
#   if(length(rho)<2) stop("rho must be a numeric vector of length 2!")
#   
#   factor_term_A<- rho[1]
#   factor_term_B <- rho[2]
#   factor_term_z <- rho[3]
#   set.seed(1) 
#   rank_true <- rank0
#   bbeta0 <- t(matrix(rnorm(d*rank_true), d, rank_true) %*% matrix(rnorm(rank_true* p), rank_true, p)) / p *4 * factor_term_z
#   
#   Blist <- list()
#   IC <- 'Orth'
#   if(IC=='Orth'){
#     set.seed(1)
#     Ztmp <- matrix(rnorm(p * (q+qs[1])), p, (q+qs[1]))
#     A <- qr(Ztmp)
#     A1 <- qr.Q(A) %*% Diag(seq(q+qs[1], 1, length=q+qs[1]))
#     A1 <- A1 %*% Diag(sign(A1[1, ]))*factor_term_A ## Fixed B0 and mu0 for each repeat.
#     A0 <- A1[,1:q]; B1 <- A1[,(q+1):(q+qs[1])]
#     t(A0) %*% A0; t(B1) %*% B1
#     Blist[[1]] <- B1
#     for(r in 2:S){
#       set.seed(r)
#       Ztmp <- matrix(rnorm(p * (qs[2])), p, (qs[2]))
#       A <- qr(Ztmp)
#       A1 <- qr.Q(A) %*% Diag(seq(qs[2], 1, length=qs[1]))  * factor_term_B
#       t(A1) %*% A1
#       Blist[[r]]   <- A1 %*% Diag(sign(A1[1, ])) 
#     }
#   }
#   
#   set.seed(seed)
#   
#   if(d<2) stop("d must be greater than 1!")
#   Zlist <- list()
#   Xlist <- list()
#   Flist <- list()
#   Hlist <- list()
#   for(s in 1:S){
#     
#     n <- nvec[s]
#     Z <- MASS::mvrnorm(n, mu=rep(0, d-1), Sigma = cor.mat(d-1, rho=0.5))
#     Z <- cbind(1, Z)
#     Zlist[[s]] <- Z
#     epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
#     
#     
#     FH <- mvrnorm(n, mu = rep(0, q+qs[s]), Diag(rep(1,q+qs[s])))
#     Flist[[s]] <- FH[,1:q]
#     Hlist[[s]] <- FH[,(q+1):(q+qs[s])]
#     AB1 <- cbind(A0, Blist[[s]])
#     mu <- exp(Z %*% t(bbeta0) + FH %*% t(AB1) + epsi) # + matrix(mu0, n, p, byrow = T) 
#     Xlist[[s]] <- matrix(rpois(n * p, lambda = mu), n, p)
#     
#   }
#   lambdavec <- rep(sigma2_eps, S)
#   
#   return(list(Xlist = Xlist, Zlist=Zlist, bbeta0=bbeta0, A0 = A0, Blist0 = Blist,
#               lambdavec= lambdavec, Flist = Flist, Hlist = Hlist, rank=rank0, q=q, qs=qs))
# }
# 


# Compared methods --------------------------------------------------------
mat2list <-function(z_int, nvec){
  
  zList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){
    
    zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
    istart <- istart + nvec[i]
  }
  return(zList_int)
}

COAP_run <- function(XcList, ZList, q, rank_use, aList=NULL, ...){
  
  require(COAP)
  nvec <- sapply(XcList, nrow)
  if(is.null(aList)){
    aList <- lapply(nvec, function(n1) rep(1, n1))
  }
  
  Xmat <- Reduce(rbind, XcList); Zmat <- Reduce( rbind, ZList)
  avec <- unlist(aList)
  tic <- proc.time()
  res_coap <- RR_COAP(Xmat, multiFac= avec, Z = Zmat,  q=q, rank_use=rank_use, ...)
  toc <- proc.time()
  res_coap$Flist <- mat2list(res_coap$H, nvec)
  res_coap$time.use <- toc[3] - tic[3]
  return(res_coap)
}

MSFR_run <- function(XcList, ZList, q, qs, maxIter=1e4, log.transform=FALSE,load.source=FALSE, dir.source=NULL){
  
  # require(MSFA)
  require(psych)
  if(!load.source){
    source(paste0(dir.source, "MSFR_main_R_MSFR_V1.R"))
  }
  #fa <- psych::fa
  B_s <- ZList
  if(log.transform){
    X_s <- lapply(XcList, function(x) log(1+x))
  }else{
    X_s <- XcList #
  }
  t1<- proc.time()
  test <- start_msfa(X_s, B_s, 5, k=q, j_s=qs, constraint = "block_lower2", method = "adhoc")
  # EM_beta <- ecm_msfa(X_s, B_s, start=test, trace = FALSE, nIt=maxIter, constraint = "block_lower1")
  EM_beta <- ecm_msfa(X_s, B_s, start=test, trace = FALSE, nIt=maxIter)
  t2<- proc.time()
  EM_beta$Flist <- lapply(EM_beta$E_f, t)
  EM_beta$Hlist <- lapply(EM_beta$E_l, t)
  EM_beta$time.use <- t2[3] - t1[3]
  EM_beta$lambdavec <- sapply(EM_beta$psi_s,mean)
  return(EM_beta)
}


# Metrics -----------------------------------------------------------------
# mean.Fnorm <- function(x) sum(x^2)/ length(x)
normvec <- function(x) sqrt(sum(x^2)/ length(x))

trace_statistic_fun <- function(H, H0){
  
  tr_fun <- function(x) sum(diag(x))
  mat1 <- t(H0) %*% H %*% qr.solve(t(H) %*% H) %*% t(H) %*% H0
  
  tr_fun(mat1) / tr_fun(t(H0) %*% H0)
  
}
trace_list_fun <- function(Hlist, H0list){
  trvec <- rep(NA, length(Hlist))
  for(i in seq_along(trvec)){
    trvec[i] <- trace_statistic_fun(Hlist[[i]], H0list[[i]])
  }
  return(mean(trvec))
}

