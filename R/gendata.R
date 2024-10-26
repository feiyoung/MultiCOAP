#' Generate simulated data
#' @description Generate simulated data from covariate-augmented Poisson factor models
#' @param seed a postive integer, the random seed for reproducibility of data generation process.
#' @param nvec a  vector with postive integers, specify the sample size in each study/source.
#' @param a_interval a numeric vector with two elements, specify the range of offset term values in each study. 
#' @param p a postive integer, specify the dimension of count variables.
#' @param d a postive integer,  specify the dimension of covariate matrix.
#' @param q a postive integer,  specify the number of study-shared factors.
#' @param qs a  vector with postive integers, specify the number of study-specified factors.
#' @param rank0 a postive integer, specify the rank of the coefficient matrix.
#' @param rho a numeric vector with length 3 and positive elements, specify the signal strength of regression coefficient and loading matrices, respectively. 
#' @param sigma2_eps a positive real, the variance of overdispersion error.
#' @param seed.beta a postive integer, the random seed for fixing the regression coefficient matrix and loading matrix generation.
#' @return return a list including the following components: 
#' (1) Xlist, the list consisting of high-dimensional count matrices from multiple studies; (2) aList: the known normalization term (offset) for each study; (3) Zlist, the list consisting of covariate matrix;
#' (4) bbeta0, the true regression coefficient matrix; (5) A0, the loading matrix of study-shared factors;  (6) Blist, the list consisting of loading matrices of study-specified factors;
#' (7)lambdavec, the variance vector of the random error vector; (8)Flist, the list composed by study-shared factor matrices; (9) Hlist, the list composed by study-specified factor matrices; 
#' (10) rank0, the rank of underlying regression coefficient matrix;  (11) q, the number of study-shared factors; (12)qs, the numbers of study-specified factors.
#' @details None
#' @seealso None
#' @references None
#' @export
#' @importFrom  MASS mvrnorm
#' @importFrom stats coef lm resid rnorm rpois runif
#' @examples 
#' seed <- 1; nvec <- c(100,300); p<- 300;
#' d <- 3; q<- 3; qs <- rep(2,2)
#' datlist <- gendata_simu_multi2(seed=seed, nvec=nvec, p=p, d=d, q=3, qs=qs)
#' str(datlist)


### Add the different a_(si) in the model
gendata_simu_multi2 <- function (seed = 1, nvec = c(100,300), a_interval=c(0,1), p = 50, d=3, q = 3,
                                 qs= rep(2, length(nvec)),
                                 rank0=3, rho = c(rhoA=1, rhoB= 1, rhoZ=1), sigma2_eps=1, seed.beta=1){
  # seed = 1; nvec = c(100,300); p = 50; d=3; q = 3
  # qs= rep(2, length(nvec))
  # rank0=3; rho = c(1, 1, 0.5); sigma2_eps=0.1;a_interval=c(0,1)
  
  
  if(length(nvec)<2) stop("nvec must have at least two elements!")
  S <- length(nvec)
  
  #require(MASS)
  if(rank0<=1) stop("rank0 must be greater than 1!")
  cor.mat<-function (p, rho, type = "toeplitz") {
    if (p == 1) 
      return(matrix(1, 1, 1))
    mat <- diag(p)
    if (type == "toeplitz") {
      for (i in 2:p) {
        for (j in 1:i) {
          mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
        }
      }
    }
    if (type == "identity") {
      mat[mat == 0] <- rho
    }
    return(mat)
  }
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
  
  if(length(rho)<2) stop("rho must be a numeric vector of length 2!")
  
  factor_term_A<- rho[1]
  factor_term_B <- rho[2]
  factor_term_z <- rho[3]
  set.seed(seed.beta) 
  rank_true <- rank0
  bbeta0 <- t(matrix(rnorm(d*rank_true), d, rank_true) %*% matrix(rnorm(rank_true* p), rank_true, p)) / p *4 * factor_term_z
  
  Blist <- list()
  IC <- 'Orth'
  if(IC=='Orth'){
    set.seed(seed.beta)
    Ztmp <- matrix(rnorm(p * (q+qs[1])), p, (q+qs[1]))
    A <- qr(Ztmp)
    A1 <- qr.Q(A) %*% Diag(seq(q+qs[1], 1, length=q+qs[1]))
    A1 <- A1 %*% Diag(sign(A1[1, ]))*factor_term_A ## Fixed B0 and mu0 for each repeat.
    A0 <- A1[,1:q]; B1 <- A1[,(q+1):(q+qs[1])]
    t(A0) %*% A0; t(B1) %*% B1
    Blist[[1]] <- B1
    for(r in 2:S){
     
      Ztmp <- matrix(rnorm(p * (qs[r])), p, (qs[r]))
      A <- qr(Ztmp)
      A1 <- qr.Q(A) %*% Diag(seq(qs[r], 1, length=qs[r]))  * factor_term_B
      t(A1) %*% A1
      Blist[[r]]   <- A1 %*% Diag(sign(A1[1, ])) 
    }
  }
  
  set.seed(seed)
  ### randomly draw a_(si) from a_interval arguments
  if(max(a_interval)<=1){
    aList <- lapply(nvec, function(n1) rep(1, n1))
  }else{
    aList <- lapply(nvec, function(n1) {
      a1 <- round(runif(n1, min=a_interval[1], max=a_interval[2]))
      a1[a1==0] <- 1
      return(a1)
    })
  }
  
  if(d<2) stop("d must be greater than 1!")
  Zlist <- list()
  Xlist <- list()
  Flist <- list()
  Hlist <- list()
  for(s in 1:S){
    # s<- 1
    n <- nvec[s]
    Z <- MASS::mvrnorm(n, mu=rep(0, d-1), Sigma = cor.mat(d-1, rho=0.5))
    Z <- cbind(1, Z)
    Zlist[[s]] <- Z
    epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
    
    
    FH <- mvrnorm(n, mu = rep(0, q+qs[s]), Diag(rep(1,q+qs[s])))
    Flist[[s]] <- FH[,1:q]
    Hlist[[s]] <- FH[,(q+1):(q+qs[s])]
    AB1 <- cbind(A0, Blist[[s]])
    mu <- exp(Z %*% t(bbeta0) + FH %*% t(AB1) + epsi) # + matrix(mu0, n, p, byrow = T)
    mu <- mu * aList[[s]] ## add the normalized factor a_si
    Xlist[[s]] <- matrix(rpois(n * p, lambda = mu), n, p)
    
  }
  lambdavec <- rep(sigma2_eps, S)
  
  return(list(Xlist = Xlist,aList=aList, Zlist=Zlist, bbeta0=bbeta0, A0 = A0, Blist0 = Blist,
              lambdavec= lambdavec, Flist = Flist, Hlist = Hlist, rank=rank0, q=q, qs=qs))
}


### Conder the absolute value of counts
gendata_simu_multi1 <- function (seed = 1, nvec = c(100,300), p = 50, d=3, q = 3,
                                 qs= rep(2, length(nvec)),
                                 rank0=3, rho = c(rhoA=1, rhoB= 1, rhoZ=1), sigma2_eps=1,seed.beta=1){
  ## rho = c(rhoA=1, rhoB= 1, rhoZ=1) control the signal strength of A, B and Z
  # seed = 1; nvec = c(100,300); p = 50; d=3; q = 3
  # qs= rep(2, length(nvec))
  # rank0=3; rho = c(1, 1); sigma2_eps=0.1;
  
  
  if(length(nvec)<2) stop("nvec must have at least two elements!")
  S <- length(nvec)
  
  #require(MASS)
  if(rank0<=1) stop("rank0 must be greater than 1!")
  cor.mat<-function (p, rho, type = "toeplitz") {
    if (p == 1) 
      return(matrix(1, 1, 1))
    mat <- diag(p)
    if (type == "toeplitz") {
      for (i in 2:p) {
        for (j in 1:i) {
          mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
        }
      }
    }
    if (type == "identity") {
      mat[mat == 0] <- rho
    }
    return(mat)
  }
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
  
  if(length(rho)<2) stop("rho must be a numeric vector of length 2!")
  
  factor_term_A<- rho[1]
  factor_term_B <- rho[2]
  factor_term_z <- rho[3]
  set.seed(seed.beta) 
  rank_true <- rank0
  bbeta0 <- t(matrix(rnorm(d*rank_true), d, rank_true) %*% matrix(rnorm(rank_true* p), rank_true, p)) / p *4 * factor_term_z
  
  Blist <- list()
  IC <- 'Orth'
  if(IC=='Orth'){
    set.seed(seed.beta)
    Ztmp <- matrix(rnorm(p * (q+qs[1])), p, (q+qs[1]))
    A <- qr(Ztmp)
    A1 <- qr.Q(A) %*% Diag(seq(q+qs[1], 1, length=q+qs[1]))
    A1 <- A1 %*% Diag(sign(A1[1, ]))*factor_term_A ## Fixed B0 and mu0 for each repeat.
    A0 <- A1[,1:q]; B1 <- A1[,(q+1):(q+qs[1])]
    t(A0) %*% A0; t(B1) %*% B1
    Blist[[1]] <- B1
    for(r in 2:S){
      #set.seed(seed.beta+r)
      Ztmp <- matrix(rnorm(p * (qs[r])), p, (qs[r]))
      A <- qr(Ztmp)
      A1 <- qr.Q(A) %*% Diag(seq(qs[r], 1, length=qs[r]))  * factor_term_B
      t(A1) %*% A1
      Blist[[r]]   <- A1 %*% Diag(sign(A1[1, ])) 
    }
  }
  
  set.seed(seed)
  
  if(d<2) stop("d must be greater than 1!")
  Zlist <- list()
  Xlist <- list()
  Flist <- list()
  Hlist <- list()
  for(s in 1:S){
    
    n <- nvec[s]
    Z <- MASS::mvrnorm(n, mu=rep(0, d-1), Sigma = cor.mat(d-1, rho=0.5))
    Z <- cbind(1, Z)
    Zlist[[s]] <- Z
    epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
    
    
    FH <- mvrnorm(n, mu = rep(0, q+qs[s]), Diag(rep(1,q+qs[s])))
    Flist[[s]] <- FH[,1:q]
    Hlist[[s]] <- FH[,(q+1):(q+qs[s])]
    AB1 <- cbind(A0, Blist[[s]])
    mu <- exp(Z %*% t(bbeta0) + FH %*% t(AB1) + epsi) # + matrix(mu0, n, p, byrow = T) 
    Xlist[[s]] <- matrix(rpois(n * p, lambda = mu), n, p)
    
  }
  lambdavec <- rep(sigma2_eps, S)
  
  return(list(Xlist = Xlist, Zlist=Zlist, bbeta0=bbeta0, A0 = A0, Blist0 = Blist,
              lambdavec= lambdavec, Flist = Flist, Hlist = Hlist, rank=rank0, q=q, qs=qs))
}



gendata_simu_multi <-function (seed = 1, nvec = c(100,300), p = 50, d=3, q = 3,
                                       qs= rep(2, length(nvec)),
                                       rank0=3, rho = c(1, 1), sigma2_eps=1,seed.beta=1){
  # seed = 1; nvec = c(100,300); p = 50; d=3; q = 3
  # qs= rep(2, length(nvec))
  # rank0=3; rho = c(1, 1); sigma2_eps=0.1;
  
  
  if(length(nvec)<2) stop("nvec must have at least two elements!")
  S <- length(nvec)
  
  # require(MASS)
  if(rank0<=1) stop("rank0 must be greater than 1!")
  cor.mat<-function (p, rho, type = "toeplitz") {
    if (p == 1) 
      return(matrix(1, 1, 1))
    mat <- diag(p)
    if (type == "toeplitz") {
      for (i in 2:p) {
        for (j in 1:i) {
          mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
        }
      }
    }
    if (type == "identity") {
      mat[mat == 0] <- rho
    }
    return(mat)
  }
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
  
  if(length(rho)<2) stop("rho must be a numeric vector of length 2!")
  
  
  factor_term_B <- rho[1]
  factor_term_z <- rho[2]
  set.seed(seed.beta) 
  rank_true <- rank0
  bbeta0 <- t(matrix(rnorm(d*rank_true), d, rank_true) %*% matrix(rnorm(rank_true* p), rank_true, p)) / p *4 * factor_term_z
  
  Blist <- list()
  IC <- 'Orth'
  if(IC=='Orth'){
    set.seed(seed.beta)
    Ztmp <- matrix(rnorm(p * (q+qs[1])), p, (q+qs[1]))
    A <- qr(Ztmp)
    A1 <- qr.Q(A) %*% Diag(seq(q+qs[1], 1, length=q+qs[1]))
    A1 <- A1 %*% Diag(sign(A1[1, ]))*factor_term_B ## Fixed B0 and mu0 for each repeat.
    A0 <- A1[,1:q]; B1 <- A1[,(q+1):(q+qs[1])]
    t(A0) %*% A0; t(B1) %*% B1
    Blist[[1]] <- B1
    for(r in 2:S){
      
      Ztmp <- matrix(rnorm(p * (qs[2])), p, (qs[2]))
      A <- qr(Ztmp)
      A1 <- qr.Q(A) %*% Diag(seq(qs[2], 1, length=qs[1]))  * factor_term_B
      t(A1) %*% A1
      Blist[[r]]   <- A1 %*% Diag(sign(A1[1, ])) 
    }
  }
  
  set.seed(seed)
  
  if(d<2) stop("d must be greater than 1!")
  Zlist <- list()
  Xlist <- list()
  Flist <- list()
  Hlist <- list()
  for(s in 1:S){
    
    n <- nvec[s]
    Z <- MASS::mvrnorm(n, mu=rep(0, d-1), Sigma = cor.mat(d-1, rho=0.5))
    Z <- cbind(1, Z)
    Zlist[[s]] <- Z
    epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
    
    
    FH <- mvrnorm(n, mu = rep(0, q+qs[s]), Diag(rep(1,q+qs[s])))
    Flist[[s]] <- FH[,1:q]
    Hlist[[s]] <- FH[,(q+1):(q+qs[s])]
    AB1 <- cbind(A0, Blist[[s]])
    mu <- exp(Z %*% t(bbeta0) + FH %*% t(AB1) + epsi) # + matrix(mu0, n, p, byrow = T) 
    Xlist[[s]] <- matrix(rpois(n * p, lambda = mu), n, p)
    
  }
  lambdavec <- rep(sigma2_eps, S)
  
  return(list(Xlist = Xlist, Zlist=Zlist, bbeta0=bbeta0, A0 = A0, Blist0 = Blist,
              lambdavec= lambdavec, Flist = Flist, Hlist = Hlist, rank=rank0, q=q, qs=qs))
}
