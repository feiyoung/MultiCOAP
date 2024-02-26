// This script implement multi-study covariate-augumented Poisson factor model.
// Date: 2023-08-17

// Revised log:


#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>

#define INT_MIN1 (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;


// Define global variables
//double X_count_ij, a_i, invLambda_j,Mu_x_ij;
/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//' 

// diag(W0* Cki * W0.t())
vec decomp(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  svd(U, s, V, Cki);
  WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}  


List irlbaCpp(const mat& X, const int& q){
  Rcpp::Environment irlba("package:irlba");
  Rcpp::Function f = irlba["irlba"];
  
  return f(Named("A") = X, Named("nv") = q);
  
}  
  
// update A matrix
double calELBO(const field<mat>& Xf, const field<vec>& af, const field<mat>& Zf, 
               const field<mat>& Muf_y, const field<mat>& Sf_y,
               const mat& A,  const field<mat>& Bf,const mat& bbeta, 
               const field<mat>& Muf_f, const field<mat>& Sf_f,
               const field<mat>& Muf_h, const field<mat>& Sf_h, const vec& invLambda){
  
  int r, r_max = Xf.n_elem, p = Xf(0).n_cols;
  
  double pois_term=0.0, entropy=0.0, ELBO=0.0;
  for(r=0; r< r_max; ++r){  
    pois_term+= accu(Xf(r) % Muf_y(r)) -accu(repmat(af(r), 1, p) % exp(Muf_y(r) + Sf_y(r)/2));
  }
  double dimreduction_term1 = 0.0, dimreduction_term2=0.0;
  int n;
  for(r=0; r< r_max; ++r){
    n = Xf(r).n_rows;
    mat Sf_bar = n* Sf_f(r);
    mat Sh_bar = n* Sf_h(r);
    mat dX = (Muf_y(r) - Zf(r)*bbeta.t() - Muf_f(r) * A.t() - Muf_h(r) * Bf(r).t()) * sqrt(invLambda(r));
    vec rowSum_ASA = decomp(Sf_bar, A)* invLambda(r);
    vec rowSum_BSB = decomp(Sh_bar, Bf(r)) * invLambda(r);
    dimreduction_term1 += -0.5*(accu(dX % dX)+ 
                                      accu(Sf_y(r))*invLambda(r)+ accu(rowSum_BSB+rowSum_ASA) - n*p* accu(log(invLambda(r))) );
    
    dimreduction_term2 += -0.5*(accu(Muf_f(r)% Muf_f(r))+accu(Muf_h(r)% Muf_h(r))+
                                      trace(Sf_bar)+trace(Sh_bar));
  }
  
  for(r=0; r< r_max; ++r){
    n = Xf(r).n_rows;
    entropy += 0.5*accu(log(Sf_y(r))) + 
      0.5 * n* (log_det_sympd(Sf_f(r)) + log_det_sympd(Sf_h(r)));
  }
 
    
  ELBO = pois_term +dimreduction_term1 + dimreduction_term2 + entropy;
  return ELBO;
}


arma::mat update_A(const field<mat>& Muf_y, const field<mat>& Zf, const mat& bbeta,
                   const field<mat>& Muf_h, const field<mat>& Bf, 
                   const field<mat>& Muf_f, const field<mat>& Sf_f){
  int r, r_max = Muf_y.n_elem, n, p = Muf_y(0).n_cols, q=Muf_f(0).n_cols;
  mat A1(q, q, fill::zeros), A2(p, q, fill::zeros);
  
  for(r=0; r<r_max; ++r){
    n = Muf_f(r).n_rows;
    A1 += Muf_f(r).t() * Muf_f(r) + n*Sf_f(r);
    A2 += trans(Muf_y(r) - Zf(r) *bbeta.t() - Muf_h(r) * Bf(r).t()) * Muf_f(r);
  }
  
  return (A2 * A1.i());
}

// update B_s
arma::field<mat> update_B(const field<mat>& Muf_y, const field<mat>& Zf, const mat& bbeta,
                   const field<mat>& Muf_f, const mat& A, 
                   const field<mat>& Muf_h, const field<mat>& Sf_h){
  int r, n, r_max = Muf_y.n_elem;
  field<mat> Bf(r_max);
  
  for(r=0; r<r_max; ++r){
    n = Muf_h(r).n_rows;
    mat A1_tmp, A2_tmp;
    A1_tmp = Muf_h(r).t() * Muf_h(r) + n*Sf_h(r);
    A2_tmp = trans(Muf_y(r) - Zf(r) *bbeta.t() - Muf_f(r) *A.t()) * Muf_h(r);
    Bf(r) = A2_tmp * A1_tmp.i();
  }
  
  return Bf;
}


// update lambda_s
arma::vec update_Lambda( const field<mat>& Zf, const field<mat>& Muf_y, const field<mat>& Sf_y,
                         const mat& A,  const field<mat>& Bf,const mat& bbeta, 
                         const field<mat>& Muf_f, const field<mat>& Sf_f,
                         const field<mat>& Muf_h, const field<mat>& Sf_h,  field<mat>& Yf_tilde){

  int r, n, r_max = Zf.n_elem,  p = A.n_rows;
  vec Lambda(r_max);
  for(r = 0; r<r_max; ++r){
    n = Muf_y(r).n_rows;
    vec a_vec = decomp(Sf_f(r), A);
    vec b_vec = decomp(Sf_h(r), Bf(r));
    mat tmpX= Muf_y(r) - Muf_f(r)*A.t() - Muf_h(r)* Bf(r).t();
    mat dX =  tmpX - Zf(r) * bbeta.t(); // speed up using Yf_tilde
    Lambda(r) = accu(dX % dX)/n/p  + accu(Sf_y(r))/n/p + accu(a_vec)/p + accu(b_vec)/p;
    Yf_tilde(r) = tmpX/sqrt(Lambda(r));
  }
  
  return Lambda;
}

// update bbeta
// separate method
mat update_bbeta_sep(const field<mat>& Zf, const field<mat>&Yf_tilde,
                     const vec& invLambda, const int& rank_use, const bool& fast_svd=true){ 
  
  // Perform singular value decomposition of X
  int r, r_max = Zf.n_elem, N=0, p=Yf_tilde(0).n_cols, d=Zf(0).n_cols;
  vec n_vec(r_max+1);
  n_vec(0) = 0;
  for(r=0; r<r_max; ++r){
    N += Zf(r).n_rows;
    n_vec(r+1) = N;
  }
  mat Y_bar(N, p);
  mat Z_bar(N, d);
  N = 0;
  for(r=0; r<r_max; ++r){
    Y_bar.rows(n_vec(r), n_vec(r+1)-1) = Yf_tilde(r);
    Z_bar.rows(n_vec(r), n_vec(r+1)-1) = Zf(r)  * sqrt(invLambda(r));
  }
  
  
  vec s;
  mat C_ls = pinv(Z_bar.t() * Z_bar) * Z_bar.t()*Y_bar; // d*p
  mat C_rr;
  //Rprintf("good2\n");
  if(d>rank_use){ // If conduct reduced rank regression, we add rank constraint.
    mat ZC = Z_bar * C_ls; // n*p
   
    mat V;
    if(fast_svd){
      Rcpp::List svdX = irlbaCpp(ZC, rank_use);
      mat V1 = svdX["v"];
      V = V1;
      V1.reset();
    }else{
      mat U, V1;
      vec s;
      svd(U, s, V1, ZC);  // how to speed up using approximate SVD
      U.reset();
      V1.reset();
    }
    C_rr = C_ls* (V.cols(0, rank_use-1)) * trans(V.cols(0, rank_use-1)) ; // d*p
  }else{
    C_rr = C_ls;
  }
  
  return C_rr.t();
  //return C_ls.t();
  
}


void add_IC_Orth(mat& A, field<mat>& Bf){
  // Try the orthogonal matrix method
  int r, q = A.n_cols, qs1 = Bf(0).n_cols, r_max=Bf.n_elem;
  mat U, V;
  vec s;
  mat AB1 = join_rows(A,Bf(0));
  svd(U, s, V, AB1);
  vec signU = sign(U.row(0).t()); // can be speed up by only caculate the first columns
  A = U.cols(0, q-1) * diagmat(s.subvec(0,q-1) % signU.subvec(0, q-1));
  Bf(0) = U.cols(q, q+qs1-1) * diagmat(s.subvec(q,q+qs1-1) % signU.subvec(q,q+qs1-1));
  
  for(r=1; r<r_max; ++r){
    qs1 = Bf(r).n_cols;
    mat U1, V1;
    vec s1;
    svd(U1, s1, V1, Bf(r));
    vec signU1 = sign(U1.row(0).t());
    Bf(r) = U1.cols(0,qs1-1) * diagmat(s1 % signU1.subvec(0,qs1-1));
  }
  
  
}

void add_IC_LT(mat& A, field<mat>& Bf){
  // Try the lower triangular matrix method
  int r, q = A.n_cols, qs1 = Bf(0).n_cols, r_max=Bf.n_elem;
  mat AB1 = join_rows(A,Bf(0));
  mat Ls = AB1.submat(0, 0, q+qs1-1, q+qs1-1 );
  mat L  = trimatl(Ls);
  AB1.submat(0, 0, q+qs1-1, q+qs1-1 ) = L;
  AB1 = AB1 * diagmat(sign(diagvec(L)));
  A = AB1.cols(0, q-1);
  Bf(0) = AB1.cols(q, q+qs1-1);
  
  for(r=1; r<r_max; ++r){
    qs1 = Bf(r).n_cols;
    mat Ls1 = Bf(r).submat(0, 0, qs1-1, qs1-1 );
    mat L1  = trimatl(Ls1);
    Bf(r).submat(0, 0, qs1-1, qs1-1 ) = L1;
    Bf(r) = Bf(r)*diagmat(sign(diagvec(L1)));
  }
  
}

void VB_Estep(const field<mat>& Xf, const field<vec>& af, const field<mat>& Zf, field<mat>& Muf_y, 
              field<mat>& Sf_y,
              const vec& invLambda, const mat& A, const field<mat>& Bf,const mat& bbeta, 
              field<mat>& Muf_f,field<mat>& Sf_f,field<mat>& Muf_h,field<mat>& Sf_h){
 
  int r, r_max = Xf.n_elem, q= A.n_cols, p = Xf(0).n_cols;
  //  ## VB E-step
  // update posterior variance of y: Sf_y
  // update posterior mean of y: Muf_y
  for(r=0; r<r_max; ++r){
    // Matrix operation for this method:
    mat tmp_mat = Zf(r) * bbeta.t() + Muf_f(r) * A.t() +Muf_h(r) * Bf(r).t();
    Muf_y(r) = (Xf(r) - repmat(af(r), 1, p) % exp(Muf_y(r)) % (1- Muf_y(r)) +  invLambda(r)* tmp_mat) / 
      (repmat(af(r), 1, p) % exp(Muf_y(r)) + invLambda(r));
    
    Sf_y(r) = 1.0 / (repmat(af(r), 1, p) % exp(Muf_y(r)) + invLambda(r));
  }
    
  // double elbo1 = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, Mu_h, S_h, Sigma_h);
  // Rprintf("VB dY= %4f\n", elbo1 - elbo0);
  
  // update posterior variance of f_i, Sf_f
  // update posterior mean of f_i, Muf_f
  vec Si_inv;
  vec si = diagvec(A.t() * A);
  for(r=0; r<r_max; ++r){
    Si_inv = invLambda(r) * si + 1.0; // diagmat(invLambda) * B
    Sf_f(r) = diagmat(1.0/Si_inv);
    Muf_f(r) = invLambda(r) * (Muf_y(r) - Zf(r) *bbeta.t() - Muf_h(r) * Bf(r).t()) * A * Sf_f(r);
  }
  // double elbo2 = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, Mu_h, S_h, Sigma_h);
  // Rprintf("VB dH= %4f\n", elbo2 - elbo1);
  
  // update posterior variance of h_i, S_i
  // update posterior mean of h_i, Mu_h
  
  for(r=0; r<r_max; ++r){
    Sf_h(r) = inv(invLambda(r) * Bf(r).t() * Bf(r) + eye(Bf(r).n_cols, Bf(r).n_cols));
    Sf_h(r) = diagmat(Sf_h(r));
    Muf_h(r) = invLambda(r) * (Muf_y(r) - Zf(r) *bbeta.t() - Muf_f(r) * A.t()) * Bf(r) * Sf_h(r);
  }
  
}



// [[Rcpp::export]]
Rcpp::List MuCOAP_cpp(const Rcpp::List& XcList, const Rcpp::List& aList, const  Rcpp::List& ZList,
                      const int& rank_use,
                       const Rcpp::List& MuList_y_int, 
                       const Rcpp::List& SList_y_int,
                       const arma::vec& invLambda_int, const arma::mat& A_int, 
                       const Rcpp::List& BList_int,
                       const arma::mat& bbeta_int,
                       const Rcpp::List& MuList_f_int,
                       const Rcpp::List& SList_f_int,
                      const Rcpp::List& MuList_h_int,
                       const Rcpp::List& SList_h_int, const int& ic,
                       const double& epsELBO, const int& maxIter, const bool& verbose,
                       const bool& loop_ic = false){
  
  
  int r,r_max = XcList.length(); // get the number of data source
  
  // Initialize
  // transfer list to field.
  field<mat> Xf(r_max), Muf_y(r_max), Sf_y(r_max), Muf_f(r_max), Sf_f(r_max);
  field<mat> Muf_h(r_max), Sf_h(r_max), Bf(r_max), Zf(r_max);
  field<vec> af(r_max);
  //Rprintf("good starting!\n");
  for(r=0; r < r_max; ++r){ 
    mat Xtmp = XcList[r]; // enforce to become a matrix.
    Xf(r) = Xtmp;
    mat Xtmp1 = MuList_y_int[r];
    Muf_y(r) = Xtmp1;
    mat Xtmp2 = SList_y_int[r];
    Sf_y(r) = Xtmp2;
    mat Xtmp3 = MuList_f_int[r];
    Muf_f(r) = Xtmp3;
    mat Xtmp4 = SList_f_int[r];
    Sf_f(r) = Xtmp4;
    mat Xtmp5 = MuList_h_int[r];
    Muf_h(r) = Xtmp5;
    mat Xtmp6 = SList_h_int[r];
    Sf_h(r) = Xtmp6;
    mat tmp = BList_int[r];
    Bf(r) = tmp;
    
    
    vec tmp_vec = aList[r];
    af(r) = tmp_vec;
    mat tmp_mat = ZList[r];
    Zf(r) = tmp_mat;
  }  
  
  vec invLambda(invLambda_int);
  mat A(A_int), bbeta(bbeta_int);
  
  //Rprintf("good starting2!\n");
  field<mat> Yf_tilde(Muf_y);
  vec ELBO_vec(maxIter), bsb, Lambda;
  ELBO_vec(0) = INT_MIN1;
  int iter;
  
  
  for(iter = 1; iter < maxIter; ++iter){
    
    
    //Rprintf("E step starting!\n");
    // VB E-step
    VB_Estep(Xf, af, Zf, Muf_y, Sf_y, invLambda, A, Bf, bbeta,Muf_f, Sf_f, Muf_h, Sf_h);
    //Rprintf("E step Finished!\n");
    
    //VB M-step
    // double elbo1 = calELBO( Xf, af, Zf, Muf_y, Sf_y, A, Bf, bbeta, Muf_f, Sf_f, Muf_h, Sf_h, invLambda);
    // Rprintf("ELOB1 = %4f!\n", elbo1);
    //update A
    A = update_A(Muf_y, Zf, bbeta, Muf_h, Bf, Muf_f, Sf_f);
    // double elbo2 = calELBO( Xf, af, Zf, Muf_y, Sf_y, A, Bf, bbeta, Muf_f, Sf_f, Muf_h, Sf_h, invLambda);
    // Rprintf("dA= %4f \n", elbo2 - elbo1);
    // update B
    Bf = update_B( Muf_y, Zf, bbeta, Muf_f, A, Muf_h, Sf_h) ;
    // double elbo3 = calELBO( Xf, af, Zf, Muf_y, Sf_y, A, Bf, bbeta, Muf_f, Sf_f, Muf_h, Sf_h, invLambda);
    // Rprintf("dB= %4f \n", elbo3 - elbo2);
    // Add identifiability conditions
    if(loop_ic){
      if(ic == 1){
        add_IC_Orth(A, Bf);
      }else if(ic==2){
        add_IC_LT(A, Bf);
      }
    }
    
    
    
    // update Lambda
    Lambda =  update_Lambda( Zf, Muf_y, Sf_y, A, Bf, bbeta, Muf_f, Sf_f, Muf_h, Sf_h, Yf_tilde);
    invLambda = 1.0 / Lambda;
    // double elbo4 = calELBO( Xf, af, Zf, Muf_y, Sf_y, A, Bf, bbeta, Muf_f, Sf_f, Muf_h, Sf_h, invLambda);
    // Rprintf("dLambda= %4f\n", elbo4 - elbo3);
    
    // update bbeta
    bbeta = update_bbeta_sep(Zf, Yf_tilde, invLambda, rank_use) ;
    
    
    // double elbo5 = calELBO( Xf, af, Zf, Muf_y, Sf_y, A, Bf, bbeta, Muf_f, Sf_f, Muf_h, Sf_h, invLambda);
    // Rprintf("dbbeta= %4f \n", elbo5 - elbo4);
  
    ELBO_vec(iter) = calELBO( Xf, af, Zf, Muf_y, Sf_y, A, Bf, bbeta, Muf_f, Sf_f, Muf_h, Sf_h, invLambda);
    
    
    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n", 
              iter +1, ELBO_vec(iter), abs(ELBO_vec(iter)  - ELBO_vec(iter-1))/ abs(ELBO_vec(iter-1)));
    }
    if(abs((ELBO_vec(iter)  - ELBO_vec(iter-1))/ ELBO_vec(iter-1)) < epsELBO) break;
    
    
  }
  
  // Add identifiability conditions
  if(!loop_ic){
    if(ic == 1){
      add_IC_Orth(A, Bf);
    }else if(ic==2){
      add_IC_LT(A, Bf);
    }
    // re-estimate the parameters.
    Lambda =  update_Lambda( Zf, Muf_y, Sf_y, A, Bf, bbeta, Muf_f, Sf_f, Muf_h, Sf_h, Yf_tilde);
    invLambda = 1.0 / Lambda;
    // update bbeta
    bbeta = update_bbeta_sep(Zf, Yf_tilde, invLambda, rank_use) ;
    VB_Estep(Xf, af, Zf, Muf_y, Sf_y, invLambda, A, Bf, bbeta,Muf_f, Sf_f, Muf_h, Sf_h);
  }
  
  
  // output return value
  List resList = List::create(
    Rcpp::Named("F") = Muf_f,
    Rcpp::Named("H") = Muf_h,
    Rcpp::Named("Sf") =Sf_f,
    Rcpp::Named("Sh") = Sf_h,
    Rcpp::Named("A") = A,
    Rcpp::Named("B") = Bf,
    Rcpp::Named("bbeta") = bbeta,
    Rcpp::Named("invLambda") = invLambda,
    Rcpp::Named("Mu_y") = Muf_y,
    Rcpp::Named("ELBO") = ELBO_vec(iter-1),
    Rcpp::Named("dELBO") = ELBO_vec(iter-1)  - ELBO_vec(iter-2),
    Rcpp::Named("ELBO_seq") = ELBO_vec.subvec(0, iter-1)
  );
  return(resList);
  
}