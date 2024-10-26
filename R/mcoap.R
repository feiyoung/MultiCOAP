# generate man files
# devtools::document()
# R CMD check --as-cran MultiCOAP_1.1.tar.gz
## usethis::use_data(dat_r2_mac)
# pkgdown::build_site()
# pkgdown::build_home()
# pkgdown::build_reference()
# pkgdown::build_article("COAPsimu")
# pkgdown::build_article("ProFASTdlpfc2")

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('html_document'))

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('pdf_document'), clean = F)

#' Fit the multi-study covariate-augmented overdispersed Poisson factor model via variational inference
#' @description Fit the high-dimensional multi-study covariate-augmented overdispersed Poisson factor model via variational inference.
#' @param XcList a length-M list with each component a count matrix, which is the observed count matrix from each source/study.
#' @param ZList a length-M list with each component a matrix that is the covariate matrix from each study.
#' @param q an optional integer, specify the number of study-shared factors; default as 15.
#' @param qs a integer vector with length M, specify the number of study-specifed factors; default as 2.
#' @param rank_use an optional integer, specify the rank of the regression coefficient matrix; default as NULL, which means that rank is the dimension of covariates in Z.
#' @param aList an optional length-M list with each component a vector, the normalization factors of each study; default as full-one vector.
#' @param init an optional string, specify the initialization method, default as "MSFRVI".
#' @param epsELBO  an optional positive vlaue, tolerance of relative variation rate of the envidence lower bound value, defualt as '1e-5'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed an optional integer, specify the random seed for reproducibility in initialization.
#' @return return a list including the following components: (1) F, a list composed by the posterior estimation of study-shared factor matrix for each study; (2) H,  a list composed by the posterior estimation of study-specified factor matrix for each study; 
#' (3) Sf, a list consisting of the posterior estimation of covariance matrix of study-shared factors for each study; (4) Sh, a list consisting of the posterior estimation of covariance matrix of study-specified factors for each study;
#' (5) A, the loading matrix corresponding to study-shared factors; (6) B, a list composed by the loading matrices corresponding to the study-specified factors;
#' (7) bbeta, the estimated regression  coefficient matrix; (8) invLambda, the inverse of the estimated variances of error;  (9) ELBO: the ELBO value when algorithm stops; (7) ELBO_seq: the sequence of ELBO values.
#' (11) qrlist, the number of factors and rank of regression coefficient matrix used in fitting;  (12) time.use, the elapsed time for model fitting.
#' @details If \code{init="MSFRVI"}, it will use the results from multi-study linear factor model as initial values; If \code{init="LFM"}, it will use the results from linear factor model by combing data from all studies as initials.
#' @seealso \code{\link{MSFRVI}}
#' @references None
#' @export
#' @useDynLib MultiCOAP, .registration = TRUE
#' @importFrom  irlba irlba
#' @importFrom  Rcpp evalCpp
#'
#' @examples 
#' seed <- 1; nvec <- c(100,300); p<- 300;
#' d <- 3; q<- 3; qs <- rep(2,2)
#' datlist <- gendata_simu_multi2(seed=seed, nvec=nvec, p=p, d=d, q=3, qs=qs)
#' fit_mcoap <- MultiCOAP(datlist$Xlist, ZList = datlist$Zlist, q=3, qs=qs, rank_use = d)
#' str(fit_mcoap)

## Add MSFRVI for initialization
MultiCOAP <- function(XcList, ZList, q=15, qs=rep(2, length(XcList)), rank_use = NULL, aList = NULL,
                      init=c("MSFRVI", "LFM"),epsELBO = 1e-5, maxIter= 30, verbose= TRUE,seed=1){
  
  # init=c("MSFRVI")
  # epsELBO = 1e-5; maxIter= 30; verbose= TRUE;seed=1
  Diag<-function (vec){
    q <- length(vec)
    if (q > 1) {
      y <- diag(vec)
    }
    else {
      y <- matrix(vec, 1, 1)
    }
    return(y)
  }
  normlize <- function(Z){
    nc <- ncol(Z)
    A <- qr(Z)
    A1 <- qr.Q(A)
    A1 <- A1 %*% Diag(sign(A1[1,])) 
    return(A1)
  }
  mat1fun <- function(X){
    matrix(1, nrow(X), ncol(X))
  }
  dimrandomfun_norm <- function(q, n){
    matrix(rnorm(prod(q*n)), nrow=n, ncol= q) ## sometimes unform is better, sometimes normal is better!!!
  }
  dimrandomfun_unif <- function(q, n){
    matrix(runif(prod(q*n), max=0.5), nrow=n, ncol= q) ## uniform distribution is better
  }
  approxPCA <- function(X, q){ ## speed the computation for initial values.
    # require(irlba) 
    n <- nrow(X)
    svdX  <- irlba(A =X, nv = q)
    PCs <- svdX$u * sqrt(n)
    loadings <- svdX$v %*% diag(svdX$d[1:q]) /sqrt(n)
    errMat <- X - PCs %*% t(loadings)
    return(list(PCs = PCs, loadings = loadings, errMat=errMat))
  }
  
  init <- match.arg(init)
  ## Basic info.
  p <- ncol(XcList[[1]])
  nvec <- sapply(XcList, nrow)
  
  
  if(is.null(rank_use)) rank_use <- ncol(ZList[[1]])
  if(is.null(aList)){
    aList <- lapply(nvec, function(n1) rep(1, n1))
  }
  
  ## Initialize the variational parameters and model parameters
  MuList_y_int <- lapply(XcList, function(x) log(x+1))
  ## consider the different a_(si) effects
  MuList_y_int <- lapply(seq_along(XcList), function(r){
    MuList_y_int[[r]] - log(aList[[r]]+1e-20) ## difference by column
  })
  
  SList_y_int <- lapply(XcList, mat1fun); invLambda_int<- rep(1, length(nvec));
  Xc <- Reduce(rbind, MuList_y_int);
  lm1 <- lm(Xc ~ 0+Reduce(rbind, ZList))
  rm(Xc)
  bbeta_int <- t(coef(lm1))
  if(init %in% c('LFM', "MSFRVI")){
    E <- resid(lm1)
    A_int <- approxPCA(E, q)$loadings
    rm(E, lm1)
    MuList_f_int <- list()
  }
  # else if(init=='COAP'){
  #   # require(COAP)
  #   Xmat <- Reduce(rbind, XcList); Zmat <- Reduce( rbind, ZList)
  #   tic <- proc.time()
  #   res_coap <- RR_COAP(Xmat, Z = Zmat,  q=q, rank_use=rank_use, verbose = FALSE)
  #   toc <- proc.time()
  #   A_int <- res_coap$B
  #   MuList_f_int <- mat2list(res_coap$H, nvec)
  # }
  
  
  Flist <- list(matrix(0, nvec[1], q), matrix(0, nvec[2], q))
  SList_f_int <- lapply(nvec, function(x) diag(rep(1,q)))
  Hlist <- list(matrix(0, nvec[1], qs[1]), matrix(0, nvec[2], qs[2]))
  SList_h_int <- lapply(qs, function(x) diag(rep(1, x)))
  set.seed(seed) # random initialization
  BList_int <- lapply(qs, dimrandomfun_norm, n=p)
  BList_int <- lapply(BList_int, normlize)
  MuList_h_int <-  list()
  
  for(i in seq_along(nvec)){
    if(init %in% c('LFM', "MSFRVI")){
      MuList_f_int[[i]] <- dimrandomfun_norm(q, nvec[i])
    }
    MuList_h_int[[i]] <- dimrandomfun_norm(qs[i], nvec[i])
  }
  ic <- 1 ## select the identifiable method
  if(init=="MSFRVI"){
    #res1 <- MSFRVI(XcList, ZList, q=q, qs=qs, rank_use=rank_use, verbose=F)
    res1 <- MSFRVI_cpp(MuList_y_int, ZList, rank_use, 
                       invLambda_int, A_int, BList_int, bbeta_int, MuList_f_int, 
                       SList_f_int, MuList_h_int, SList_h_int, ic, epsELBO, maxIter, 
                       verbose=FALSE, loop_ic=TRUE) 
    A_int <- res1$A; BList_int <- res1$B
    MuList_f_int <- res1$F; MuList_h_int <- res1$H
    SList_f_int <- res1$Sf; SList_h_int <- res1$Sh
  }
  
  tic <- proc.time()
  res <- MuCOAP_cpp(XcList, aList, ZList, rank_use, MuList_y_int, SList_y_int, 
                    invLambda_int, A_int, BList_int, bbeta_int, MuList_f_int, 
                    SList_f_int, MuList_h_int, SList_h_int, ic, epsELBO, maxIter, 
                    verbose, loop_ic=TRUE) 
  toc <- proc.time()
  res$Lambda <- 1.0/ res$invLambda
  res$qrlist <- list(q=q, qs=qs, r= rank_use)
  res$time.use <- toc[3] - tic[3]
  
  return(res)
}


#' Fit the multi-study covariate-augmented linear factor model via variational inference
#' @description Fit the multi-study covariate-augmented linear factor model via variational inference
#' @param XList A length-M list, where each component represents a matrix and is the observed response matrix from each source/study. Ideally, each matrix should be continuous.
#' @param ZList a length-M list with each component a matrix that is the covariate matrix from each study.
#' @param q an optional integer, specify the number of study-shared factors; default as 15.
#' @param qs a integer vector with length M, specify the number of study-specifed factors; default as 2.
#' @param rank_use an optional integer, specify the rank of the regression coefficient matrix; default as NULL, which means that rank is the dimension of covariates in Z.
#' @param aList an optional length-M list with each component a vector, the normalization factors of each study; default as full-one vector.
#' @param epsELBO  an optional positive vlaue, tolerance of relative variation rate of the envidence lower bound value, defualt as '1e-5'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed an optional integer, specify the random seed for reproducibility in initialization.
#' @return return a list including the following components: (1) F, a list composed by the posterior estimation of study-shared factor matrix for each study; (2) H,  a list composed by the posterior estimation of study-specified factor matrix for each study; 
#' (3) Sf, a list consisting of the posterior estimation of covariance matrix of study-shared factors for each study; (4) Sh, a list consisting of the posterior estimation of covariance matrix of study-specified factors for each study;
#' (5) A, the loading matrix corresponding to study-shared factors; (6) B, a list composed by the loading matrices corresponding to the study-specified factors;
#' (7) bbeta, the estimated regression  coefficient matrix; (8) invLambda, the inverse of the estimated variances of error;  (9) ELBO: the ELBO value when algorithm stops; (7) ELBO_seq: the sequence of ELBO values.
#' (11) qrlist, the number of factors and rank of regression coefficient matrix used in fitting;  (12) time.use, the elapsed time for model fitting.
#' @details None
#' @seealso \code{\link{MultiCOAP}}
#' @references None
#' @export
#' @useDynLib MultiCOAP, .registration = TRUE
#' @importFrom  irlba irlba
#' @importFrom  Rcpp evalCpp
#'
#' @examples 
#' seed <- 1; nvec <- c(100,300); p<- 300;
#' d <- 3; q<- 3; qs <- rep(2,2)
#' datlist <- gendata_simu_multi2(seed=seed, nvec=nvec, p=p, d=d, q=3, qs=qs)
#' XList <- lapply(datlist$Xlist,  function(x) log(1+x))
#' fit_msfavi <- MSFRVI(XList, ZList = datlist$Zlist, q=3, qs=qs, rank_use = d)
#' str(fit_msfavi)
### MSFR using variational inference
MSFRVI <- function(XList, ZList, q=15, qs=rep(2, length(XList)), rank_use = NULL, aList=NULL, epsELBO = 1e-5, maxIter= 30, verbose= TRUE,seed=1){
  
  Diag<-function (vec){
    q <- length(vec)
    if (q > 1) {
      y <- diag(vec)
    }
    else {
      y <- matrix(vec, 1, 1)
    }
    return(y)
  }
  normlize <- function(Z){
    nc <- ncol(Z)
    A <- qr(Z)
    A1 <- qr.Q(A)
    A1 <- A1 %*% Diag(sign(A1[1,])) 
    return(A1)
  }
  mat1fun <- function(X){
    matrix(1, nrow(X), ncol(X))
  }
  dimrandomfun_norm <- function(q, n){
    matrix(rnorm(prod(q*n)), nrow=n, ncol= q) ## sometimes unform is better, sometimes normal is better!!!
  }
  dimrandomfun_unif <- function(q, n){
    matrix(runif(prod(q*n), max=0.5), nrow=n, ncol= q) ## uniform distribution is better
  }
  approxPCA <- function(X, q){ ## speed the computation for initial values.
    # require(irlba) 
    n <- nrow(X)
    svdX  <- irlba(A =X, nv = q)
    PCs <- svdX$u * sqrt(n)
    loadings <- svdX$v %*% diag(svdX$d[1:q]) /sqrt(n)
    errMat <- X - PCs %*% t(loadings)
    return(list(PCs = PCs, loadings = loadings, errMat=errMat))
  }
  
  
  ## Basic info.
  p <- ncol(XList[[1]])
  nvec <- sapply(XList, nrow)
  
  
  if(is.null(rank_use)) rank_use <- ncol(ZList[[1]])
  if(is.null(aList)){
    aList <- lapply(nvec, function(n1) rep(1, n1))
  }
  
  
  ## Initialize the variational parameters and model parameters
  MuList_y_int <- XList #lapply(XList, function(x) log(x+1))
  ## consider the different a_(si) effects
  MuList_y_int <- lapply(seq_along(XList), function(r){
    #MuList_y_int[[r]] - log(aList[[r]]+1e-20) ## difference by column
    MuList_y_int[[r]] - aList[[r]] ## difference by column
  })
  
  SList_y_int <- lapply(XList, mat1fun); invLambda_int<- rep(1, length(nvec));
  Xc <- Reduce(rbind, MuList_y_int);
  lm1 <- lm(Xc ~ 0+Reduce(rbind, ZList))
  rm(Xc)
  bbeta_int <- t(coef(lm1))
  init <- "LFM"
  if(init=='LFM'){
    E <- resid(lm1)
    A_int <- approxPCA(E, q)$loadings
    rm(E, lm1)
    MuList_f_int <- list()
  }
  
  
  Flist <- list(matrix(0, nvec[1], q), matrix(0, nvec[2], q))
  
  SList_f_int <- lapply(nvec, function(x) diag(rep(1,q)))
  Hlist <- list(matrix(0, nvec[1], qs[1]), matrix(0, nvec[2], qs[2]))
  SList_h_int <- lapply(qs, function(x) diag(rep(1, x)))
  set.seed(seed) # random initialization
  BList_int <- lapply(qs, dimrandomfun_norm, n=p)
  BList_int <- lapply(BList_int, normlize)
  MuList_h_int <-  list()
  
  for(i in seq_along(nvec)){
    if(init=='LFM'){
      MuList_f_int[[i]] <- dimrandomfun_norm(q, nvec[i])
    }
    MuList_h_int[[i]] <- dimrandomfun_norm(qs[i], nvec[i])
  }
  
  
  ic <- 1 ## select the identifiable method
  tic <- proc.time()
  res <- MSFRVI_cpp(MuList_y_int, ZList, rank_use, 
                    invLambda_int, A_int, BList_int, bbeta_int, MuList_f_int, 
                    SList_f_int, MuList_h_int, SList_h_int, ic, epsELBO, maxIter, 
                    verbose, loop_ic=TRUE) 
  toc <- proc.time()
  res$Lambda <- 1.0/ res$invLambda
  res$qrlist <- list(q=q, qs=qs, r= rank_use)
  res$time.use <- toc[3] - tic[3]
  
  return(res)
}



# ## Change the method to get initial values
# MultiCOAPv2 <- function(XcList, ZList, q=15, qs=rep(2, length(XcList)), rank_use = NULL, aList = NULL,
#                       init=c("LFM", "COAP"),epsELBO = 1e-5, maxIter= 30, verbose= TRUE,seed=1){
#   
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
#   normlize <- function(Z){
#     nc <- ncol(Z)
#     A <- qr(Z)
#     A1 <- qr.Q(A)
#     A1 <- A1 %*% Diag(sign(A1[1,])) 
#     return(A1)
#   }
#   mat1fun <- function(X){
#     matrix(1, nrow(X), ncol(X))
#   }
#   dimrandomfun_norm <- function(q, n){
#     matrix(rnorm(prod(q*n)), nrow=n, ncol= q) ## sometimes unform is better, sometimes normal is better!!!
#   }
#   dimrandomfun_unif <- function(q, n){
#     matrix(runif(prod(q*n), max=0.5), nrow=n, ncol= q) ## uniform distribution is better
#   }
#   approxPCA <- function(X, q){ ## speed the computation for initial values.
#     require(irlba) 
#     n <- nrow(X)
#     svdX  <- irlba(A =X, nv = q)
#     PCs <- svdX$u * sqrt(n)
#     loadings <- svdX$v %*% diag(svdX$d[1:q]) /sqrt(n)
#     errMat <- X - PCs %*% t(loadings)
#     return(list(PCs = PCs, loadings = loadings, errMat=errMat))
#   }
#   
#   init <- match.arg(init)
#   ## Basic info.
#   p <- ncol(XcList[[1]])
#   nvec <- sapply(XcList, nrow)
#   
#   
#   if(is.null(rank_use)) rank_use <- ncol(ZList[[1]])
#   if(is.null(aList)){
#     aList <- lapply(nvec, function(n1) rep(1, n1))
#   }
#   
#   ## Initialize the variational parameters and model parameters
#   MuList_y_int <- lapply(XcList, function(x) log(x+1))
#   SList_y_int <- lapply(XcList, mat1fun); invLambda_int<- rep(1, length(nvec));
#   Xc <- Reduce(rbind, MuList_y_int);
#   lm1 <- lm(Xc ~ 0+Reduce(rbind, ZList))
#   rm(Xc)
#   bbeta_int <- t(coef(lm1))
#   if(init=='LFM'){
#     E <- resid(lm1)
#     A_int <- approxPCA(E, q)$loadings
#     rm(E, lm1)
#     MuList_f_int <- list()
#   }else if(init=='COAP'){
#     require(COAP)
#     Xmat <- Reduce(rbind, XcList); Zmat <- Reduce( rbind, ZList)
#     tic <- proc.time()
#     res_coap <- RR_COAP(Xmat, Z = Zmat,  q=q, rank_use=rank_use, verbose = FALSE)
#     toc <- proc.time()
#     A_int <- res_coap$B
#     MuList_f_int <- mat2list(res_coap$H, nvec)
#   }
#   
#   
#   Flist <- list(matrix(0, nvec[1], q), matrix(0, nvec[2], q))
#   
#   SList_f_int <- lapply(nvec, function(x) diag(rep(1,q)))
#   Hlist <- list(matrix(0, nvec[1], qs[1]), matrix(0, nvec[2], qs[2]))
#   SList_h_int <- lapply(qs, function(x) diag(rep(1, x)))
#   set.seed(seed) # random initialization
#   BList_int <- lapply(qs, dimrandomfun_norm, n=p)
#   BList_int <- lapply(BList_int, normlize)
#   MuList_h_int <-  list()
#   
#   for(i in seq_along(nvec)){
#     if(init=='LFM'){
#       MuList_f_int[[i]] <- dimrandomfun_norm(q, nvec[i])
#     }
#     MuList_h_int[[i]] <- dimrandomfun_norm(qs[i], nvec[i])
#   }
#   
#   
#   ic <- 1 ## select the identifiable method
#   tic <- proc.time()
#   res <- MuCOAP_cpp(XcList, aList, ZList, rank_use, MuList_y_int, SList_y_int, 
#                     invLambda_int, A_int, BList_int, bbeta_int, MuList_f_int, 
#                     SList_f_int, MuList_h_int, SList_h_int, ic, epsELBO, maxIter, 
#                     verbose, loop_ic=TRUE) 
#   toc <- proc.time()
#   res$Lambda <- 1.0/ res$invLambda
#   res$qrlist <- list(q=q, qs=qs, r= rank_use)
#   res$time.use <- toc[3] - tic[3]
#   
#   return(res)
# }


## Version v1
# MultiCOAP <- function(XcList, ZList, q=15, qs=rep(2, length(XcList)), rank_use = NULL, aList = NULL,
#                        epsELBO = 1e-5, maxIter= 30, verbose= TRUE,seed=1){
#   
#   mat1fun <- function(X){
#     matrix(1, nrow(X), ncol(X))
#   }
#   dimrandomfun <- function(q, n){
#     matrix(rnorm(prod(q*n)), nrow=n, ncol= q)
#   }
#   approxPCA <- function(X, q){ ## speed the computation for initial values.
#     require(irlba) 
#     n <- nrow(X)
#     svdX  <- irlba(A =X, nv = q)
#     PCs <- svdX$u * sqrt(n)
#     loadings <- svdX$v %*% diag(svdX$d[1:q]) /sqrt(n)
#     errMat <- X - PCs %*% t(loadings)
#     return(list(PCs = PCs, loadings = loadings, errMat=errMat))
#   }
#   
#   ## Basic info.
#   p <- ncol(XcList[[1]])
#   nvec <- sapply(XcList, nrow)
#   
#   
#   if(is.null(rank_use)) rank_use <- ncol(ZList[[1]])
#   if(is.null(aList)){
#     aList <- lapply(nvec, function(n1) rep(1, n1))
#   }
#   
#   ## Initialize the variational parameters and model parameters
#   MuList_y_int <- lapply(XcList, function(x) log(x+1))
#   SList_y_int <- lapply(XcList, mat1fun); invLambda_int<- rep(1, length(nvec));
#   Xc <- Reduce(rbind, MuList_y_int);
#   lm1 <- lm(Xc ~ 0+Reduce(rbind, ZList))
#   rm(Xc)
#   bbeta_int <- t(coef(lm1))
#   E <- resid(lm1)
#   A_int <- approxPCA(E, q)$loadings
#   rm(E, lm1)
#   
#   Flist <- list(matrix(0, nvec[1], q), matrix(0, nvec[2], q))
#   
#   SList_f_int <- lapply(nvec, function(x) diag(rep(1,q)))
#   Hlist <- list(matrix(0, nvec[1], qs[1]), matrix(0, nvec[2], qs[2]))
#   SList_h_int <- lapply(qs, function(x) diag(rep(1, x)))
#   set.seed(seed) # random initialization
#   BList_int <- lapply(qs, dimrandomfun, n=p)
#   MuList_h_int <- MuList_f_int <- list()
#   for(i in seq_along(nvec)){
#     
#     MuList_f_int[[i]] <- dimrandomfun(q, nvec[i])
#     MuList_h_int[[i]] <- dimrandomfun(qs[i], nvec[i])
#   }
#   
#   
#   ic <- 1 ## select the identifiable method
#   tic <- proc.time()
#   res <- MuCOAP_cpp(XcList, aList, ZList, rank_use, MuList_y_int, SList_y_int, 
#                     invLambda_int, A_int, BList_int, bbeta_int, MuList_f_int, 
#                     SList_f_int, MuList_h_int, SList_h_int, ic, epsELBO, maxIter, 
#                     verbose, loop_ic=TRUE) 
#   toc <- proc.time()
#   res$Lambda <- 1.0/ res$invLambda
#   res$qrlist <- list(q=q, qs=qs, r= rank_use)
#   res$time.use <- toc[3] - tic[3]
#   
#   return(res)
# }
# 
