// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace arma;
using namespace std;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::vec rowsum_Mat(arma::mat M) {
  int nr=M.n_rows;
  arma::vec out(nr);
  for(int i=0;i<nr;i++){
    out(i)=sum(M.row(i));
  }
  return out;
}

// [[Rcpp::export]]
arma::vec colsum_Mat(arma::mat M) {
  int nc=M.n_cols;
  arma::vec out(nc);
  for(int i=0;i<nc;i++){
    out(i)=sum(M.col(i));
  }
  return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// gamma, Tau update function, gradient and hessian functions and ELBO convergence function

// [[Rcpp::export]]
arma::cube gamma_update_HMM_stat_undir(arma::mat gamma, arma::vec pi, arma::vec theta, arma::mat network, int N, int K){
  arma::cube quad_lin_coeff(N,K,2);
  for(int i = 0; i < N; i++){
    if(i!=(N-1)){
      for(int k = 0; k < K; k++){
        float t1=0;
        for(int j = i+1; j < N; j++){
          for(int l = 0; l < K; l++){
            float exp_val=exp(theta(k)+theta(l));
            t1+=(gamma(j,l)/(2*gamma(i,k)))*((network(i,j)*(theta(k)+theta(l)))-log(1+exp_val));
          }
        }
        quad_lin_coeff(i,k,0)=t1-(1/gamma(i,k));
        quad_lin_coeff(i,k,1)=log(pi(k))-log(gamma(i,k))+1;
      }
    } else if(i==(N-1)){
      for(int k = 0; k < K; k++){
        quad_lin_coeff(i,k,0)=-(1/gamma((N-1),k));
        quad_lin_coeff(i,k,1)=log(pi(k))-log(gamma((N-1),k))+1;
      }
    }
  }
  return quad_lin_coeff;
}

// [[Rcpp::export]]
arma::mat grad_HMM_stat_undir(arma::vec theta, arma::mat gamma, arma::mat network, int N, int K){
  arma::vec grad_vector(K);
  arma::mat grad_mat(K,K,fill::zeros);
  for(int i = 0; i < (N-1); i++){
    for(int j = i+1; j < N; j++){
      arma::mat grad_matsub(K,K);
      for(int k = 0; k < K; k++){
        for(int l = 0; l < K; l++){
          float exp_val=exp(theta(k)+theta(l));
          grad_matsub(k,l)=gamma(i,k)*gamma(j,l)*(network(i,j)-(exp_val/(1+exp_val)));
        }
      }
      grad_mat+=grad_matsub;
    }
  }
  arma::vec rsum=rowsum_Mat(grad_mat);
  arma::vec csum=colsum_Mat(grad_mat);
  for(int k = 0; k < K; k++){
    grad_vector(k)=rsum(k)+csum(k);
  }
  return std::move(grad_vector);
}

// [[Rcpp::export]]
arma::mat hess_HMM_stat_undir(arma::vec theta, arma::mat gamma, int N, int K){
  arma::mat t1(K,K);
  arma::mat hess_mat(K,K,fill::zeros);
  for(int i = 0; i < (N-1); i++){
    for(int j = i+1; j < N; j++){
      arma::mat hess_matsub(K,K);
      for(int k = 0; k < K; k++){
        for(int l = 0; l < K; l++){
          float exp_val=exp(theta(k)+theta(l));
          hess_matsub(k,l)=-(gamma(i,k)*gamma(j,l)*(exp_val/pow((1+exp_val),2)));
        }
      }
      hess_mat+=hess_matsub;
    }
  }
  for(int k = 0; k < K; k++){
    for(int l = 0; l < K; l++){
      if(k!=l){
        t1(k,l)=hess_mat(k,l)+hess_mat(l,k);
      }
    }
  }
  arma::vec rsum=rowsum_Mat(hess_mat);
  arma::vec csum=colsum_Mat(hess_mat);
  for(int k = 0; k < K; k++){
    t1(k,k)=(csum(k)+rsum(k)+(2*hess_mat(k,k)));
  }
  return t1;
}

// [[Rcpp::export]]
float ELBO_conv_HMM_stat_undir(arma::mat gamma, arma::vec pi, arma::vec theta, arma::mat network, int N, int K){
  float t1=0;
  for(int i = 0; i < (N-1); i++){
    for(int j = i+1; j < N; j++){
      for(int k = 0; k < K; k++){
        for(int l = 0; l < K; l++){
          float exp_val=exp(theta(k)+theta(l));
          t1+=(gamma(i,k)*gamma(j,l)*(network(i,j)*(theta(k)+theta(l))-log(1+exp_val)));
        }
      }
    }
  }
  float t2=0;
  for(int i = 0; i < N; i++){
    for(int k = 0; k < K; k++){
      if((pi(k)>=(pow(10,(-100))))&(gamma(i,k)>=(pow(10,(-100))))){
        t2+=gamma(i,k)*(log(pi(k))-log(gamma(i,k)));
      }
    }
  }
  
  float ELBO_val=t1+t2;
  return ELBO_val;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Defining functions for K=1

// [[Rcpp::export]]
float grad_HMM_stat_undir_K1(float theta, arma::mat network, int N){
  float exp_val=exp(2*theta);
  float exp_val_ratio=exp_val/(1+exp_val);
  float grad_val=0;
  for(int i = 0; i < (N-1); i++){
    for(int j = i+1; j < N; j++){
      grad_val+=network(i,j)-exp_val_ratio;
    }
  }
  float grad_val_final=2*grad_val;
  return grad_val_final;
}

// [[Rcpp::export]]
float hess_HMM_stat_undir_K1(float theta, int N){
  float exp_val=exp(2*theta);
  float exp_val_ratio=exp_val/(pow((1+exp_val),2));
  float hess_val=0;
  for(int i = 0; i < (N-1); i++){
    for(int j = i+1; j < N; j++){
      hess_val+=(-exp_val_ratio);
    }
  }
  float hess_val_final=4*hess_val;
  return hess_val_final;
}

// [[Rcpp::export]]
float ELBO_conv_HMM_stat_undir_K1(float theta, arma::mat network, int N){
  float exp_val=exp(2*theta);
  float log_exp_val=log(1+exp_val);
  float ELBO_val=0;
  for(int i = 0; i < (N-1); i++){
    for(int j = i+1; j < N; j++){
      ELBO_val+=(network(i,j)*(2*theta)-log_exp_val);
    }
  }
  return ELBO_val;
}
///////////////////////////////////////////////////////////////////////////////////////////////
