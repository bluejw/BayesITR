//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp; 


const double log2pi = log(2.0*M_PI);


// [[Rcpp::export]]
double dmvn_rcpp(rowvec& x, rowvec& mean, mat& sigma, bool logd = false){ 
  
  // calculate density of multivariate normal distribution
  // args: x: row vector data
  //      mean: row vector mean, sigma: covariance matrix  
  //      logd: true for taking log
  // returns: out: pdf (or log pdf) of multivariate normal distribution
  
  int xdim = x.size(); 
  mat rooti = trans(inv(trimatu(chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0)*log2pi;
  
  vec z = rooti*trans(x-mean);
  double out = constants-0.5*sum(z%z)+rootisum;
  
  if (logd == false){ out = exp(out); }
  return(out);
}


// [[Rcpp::export]]
mat rmvn_rcpp(const int n, vec& mean, mat& sigma){
  
  // randomly generate samples from multivariate normal distribution
  // args: n: number of data 
  //      mean: row vector mean, sigma: covariance matrix  
  // returns: out: random samples from multivariate normal distribution
  
  int k = sigma.n_cols; // dimension of the multivariate normal distribution
  mat z = randn(n, k);
  mat out = repmat(mean,1,n).t()+z*chol(sigma);
  return(out);
}


// [[Rcpp::export]]
double rinvgamma_rcpp(const double a, const double b){
  
  // generate random samples from inverse-gamma distribution
  // args: inverse-gamma(a, b)
  // returns: random sample from inverse-gamma distribution
  
  return(1/R::rgamma(a, 1/b));
}


// [[Rcpp::export]]
double rinvgaussian_rcpp(const double mu, const double lambda){
  
  // generate random samples from inverse-gaussian distribution
  // args: inverse-gaussian(mu, lambda)
  // returns: random sample from inverse-gaussian distribution
  
  double out;
  double nu = R::rnorm(0, 1);
  double y = pow(nu, 2);
  double x = mu + ((pow(mu,2)*y)/(2*lambda)) - (mu/(2*lambda))*
    sqrt(4*mu*lambda*y + pow(mu,2)*pow(y,2));
  double z = R::runif(0, 1);
  if (z > (mu/(mu + x))){ out = mu*mu/x; 
  }else{ out = x; }
  return(out);
}


// [[Rcpp::export]]
double update_eta_gamma_rcpp(vec& zeta_gamma, double xi_gamma, double nu_gamma, 
                             cube& alpha, mat& beta, vec& U, vec& lambda, vec& mu, mat& tau, vec& sigma2, 
                             const int I, const int S_sum, const int B, const int P, 
                             vec& Y, mat& A_bs, cube& X_bs, mat& Z){
  
  double eta_gamma_update;
  
  double A_sum = 0; double Ay_sum = 0;
  for (int i=0; i<I; i++){
    // calculate A_sum
    rowvec A_bs_i = A_bs.submat(i, 0, i, B-1);
    double A_zeta_gamma = conv_to<double>::from(A_bs_i*zeta_gamma);
    A_sum += pow(A_zeta_gamma, 2)*tau(i,0);
    
    // calculate Ay_sum
    double Y_i_tilde = Y(i) - mu(0) - U(i)*lambda(0);
    rowvec Z_i = Z.submat(i, 0, i, P-1);
    vec beta_q = conv_to<vec>::from(beta.row(0));
    Y_i_tilde -= conv_to<double>::from(Z_i*beta_q);
    for (int s=0; s<S_sum; s++){
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      vec alpha_s = alpha.subcube(0, s, 0, 0, s, B-1);
      Y_i_tilde -= conv_to<double>::from(X_bs_is*alpha_s);
    }
    Ay_sum += A_zeta_gamma*tau(i,0)*Y_i_tilde;
  }
  
  // Gaussian posterior
  double V_n = 1/(A_sum/sigma2(0) + 1/(xi_gamma*nu_gamma));
  double mu_n = V_n*Ay_sum/sigma2(0);
  eta_gamma_update = R::rnorm(mu_n, sqrt(V_n));
  
  return(eta_gamma_update);
}


// [[Rcpp::export]]
vec update_zeta_gamma_rcpp(double eta_gamma, vec& m_gamma,
                           cube& alpha, mat& beta, vec& U, vec& lambda, vec& mu, mat& tau, vec& sigma2, 
                           const int I, const int S_sum, const int B, const int P, 
                           vec& Y, mat& A_bs, cube& X_bs, mat& Z){
  
  vec zeta_gamma_update(B);
  
  mat A_sum(B, B); A_sum.fill(0);
  vec Ay_sum(B); Ay_sum.fill(0);
  
  for (int i=0; i<I; i++){
    // calculate A_sum
    rowvec A_bs_i = A_bs.submat(i, 0, i, B-1);
    rowvec A_eta_gamma = A_bs_i*eta_gamma;
    A_sum += A_eta_gamma.t()*A_eta_gamma*tau(i,0);
    
    // calculate Ay_sum
    double Y_i_tilde = Y(i) - mu(0) - U(i)*lambda(0);
    rowvec Z_i = Z.submat(i, 0, i, P-1);
    vec beta_q = conv_to<vec>::from(beta.row(0));
    Y_i_tilde -= conv_to<double>::from(Z_i*beta_q);
    for (int s=0; s<S_sum; s++){
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      vec alpha_s = alpha.subcube(0, s, 0, 0, s, B-1);
      Y_i_tilde -= conv_to<double>::from(X_bs_is*alpha_s);
    }
    Ay_sum += conv_to<vec>::from(A_eta_gamma*tau(i,0)*Y_i_tilde);
  }
  
  // Gaussian posterior
  mat I_B = eye(B, B);
  mat V_n = inv_sympd(A_sum/sigma2(0) + I_B);
  vec mu_n = V_n*(Ay_sum/sigma2(0) + m_gamma);
  zeta_gamma_update = conv_to<vec>::from(rmvn_rcpp(1, mu_n, V_n));
  
  return(zeta_gamma_update);
}


// [[Rcpp::export]]
mat update_eta_alpha_rcpp(vec& gamma, mat& eta_alpha, cube& zeta_alpha, mat& xi_alpha, mat& nu_alpha, 
                          mat& beta, vec& U, vec& lambda, vec& mu, mat& tau, vec& sigma2, 
                          const int I, const int Q, const int S, const int S_sum, const int B, const int P,
                          vec& Y, vec& A, mat& A_bs, cube& X_bs, mat& Z){
  
  mat eta_alpha_update = eta_alpha;
  
  // calcualte alpha_update matrix 
  cube alpha_update(Q, S_sum, B);
  for (int q=0; q<Q; q++){
    for (int s=0; s<S_sum; s++){
      double eta_alpha_update_qs = eta_alpha_update(q,s);
      vec zeta_alpha_qs = zeta_alpha.subcube(q, s, 0, q, s, B-1);
      alpha_update.subcube(q, s, 0, q, s, B-1) = eta_alpha_update_qs*zeta_alpha_qs;
    }
  }
  
  for (int s=0; s<S_sum; s++){
    
    vec zeta_alpha_s = zeta_alpha.subcube(0, s, 0, 0, s, B-1);
    
    double X_sum = 0; double Xy_sum = 0;
    for (int i=0; i<I; i++){
      // calculate X_sum
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      double X_zeta_alpha = conv_to<double>::from(X_bs_is*zeta_alpha_s);
      X_sum += pow(X_zeta_alpha, 2)*tau(i,0);
      
      // calculate Xy_sum
      double Y_i_tilde = Y(i) - mu(0) - U(i)*lambda(0);
      rowvec A_bs_i = A_bs.submat(i, 0, i, B-1);
      Y_i_tilde -= conv_to<double>::from(A_bs_i*gamma);
      rowvec Z_i = Z.submat(i, 0, i, P-1);
      vec beta_q = conv_to<vec>::from(beta.row(0));
      Y_i_tilde -= conv_to<double>::from(Z_i*beta_q);
      for (int s2=0; s2<S_sum; s2++){
        rowvec X_bs_is2 = X_bs.subcube(i, s2, 0, i, s2, B-1);
        vec alpha_update_s2 = alpha_update.subcube(0, s2, 0, 0, s2, B-1);
        Y_i_tilde -= conv_to<double>::from(X_bs_is2*alpha_update_s2);
      }
      vec alpha_update_s = alpha_update.subcube(0, s, 0, 0, s, B-1);
      Y_i_tilde += conv_to<double>::from(X_bs_is*alpha_update_s);
      Xy_sum += X_zeta_alpha*tau(i,0)*Y_i_tilde;
    }
    
    // Gaussian posterior
    double V_n = 1/(X_sum/sigma2(0) + 1/(xi_alpha(0,s)*nu_alpha(0,s)));
    double mu_n = V_n*Xy_sum/sigma2(0);
    eta_alpha_update(0,s) = R::rnorm(mu_n, sqrt(V_n));
    
    // update alpha
    alpha_update.subcube(0, s, 0, 0, s, B-1) = eta_alpha_update(0,s)*zeta_alpha_s;
  }
  
  for (int s=0; s<S; s++){
    
    vec zeta_alpha_s = zeta_alpha.subcube(1, s, 0, 1, s, B-1);
    
    double X_sum = 0; double Xa_sum = 0;
    for (int i=0; i<I; i++){
      // calculate X_sum
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      double X_zeta_alpha = conv_to<double>::from(X_bs_is*zeta_alpha_s);
      X_sum += pow(X_zeta_alpha, 2)*tau(i,1);
      
      // calculate Xa_sum
      double A_i_tilde = A(i) - mu(1) - U(i)*lambda(1);
      rowvec Z_i = Z.submat(i, 0, i, P-1);
      vec beta_q = conv_to<vec>::from(beta.row(1));
      A_i_tilde -= conv_to<double>::from(Z_i*beta_q);
      for (int s2=0; s2<S; s2++){
        rowvec X_bs_is2 = X_bs.subcube(i, s2, 0, i, s2, B-1);
        vec alpha_update_s2 = alpha_update.subcube(1, s2, 0, 1, s2, B-1);
        A_i_tilde -= conv_to<double>::from(X_bs_is2*alpha_update_s2);
      }
      vec alpha_update_s = alpha_update.subcube(1, s, 0, 1, s, B-1);
      A_i_tilde += conv_to<double>::from(X_bs_is*alpha_update_s);
      Xa_sum += X_zeta_alpha*tau(i,1)*A_i_tilde;
    }
    
    // Gaussian posterior
    double V_n = 1/(X_sum/sigma2(1) + 1/(xi_alpha(1,s)*nu_alpha(1,s)));
    double mu_n = V_n*Xa_sum/sigma2(1);
    eta_alpha_update(1,s) = R::rnorm(mu_n, sqrt(V_n));
    
    // update alpha
    alpha_update.subcube(1, s, 0, 1, s, B-1) = eta_alpha_update(1,s)*zeta_alpha_s;
  }
  
  for (int s=S; s<S_sum; s++){
    eta_alpha_update(1,s) = 0;
  }
  
  return(eta_alpha_update);
}


// [[Rcpp::export]]
cube update_zeta_alpha_rcpp(vec& gamma, mat& eta_alpha, cube& zeta_alpha, cube& m_alpha, 
                            mat& beta, vec& U, vec& lambda, vec& mu, mat& tau, vec& sigma2, 
                            const int I, const int Q, const int S, const int S_sum, const int B, const int P,
                            vec& Y, vec& A, mat& A_bs, cube& X_bs, mat& Z){
  
  cube zeta_alpha_update = zeta_alpha;
  
  // calcualte alpha_update matrix 
  cube alpha_update(Q, S_sum, B);
  for (int q=0; q<Q; q++){
    for (int s=0; s<S_sum; s++){
      double eta_alpha_qs = eta_alpha(q,s);
      vec zeta_alpha_update_qs = zeta_alpha_update.subcube(q, s, 0, q, s, B-1);
      alpha_update.subcube(q, s, 0, q, s, B-1) = eta_alpha_qs*zeta_alpha_update_qs;
    }
  }
  
  for (int s=0; s<S_sum; s++){
    
    double eta_alpha_s = eta_alpha(0,s);
    
    mat X_sum(B, B); X_sum.fill(0);
    vec Xy_sum(B); Xy_sum.fill(0);
    for (int i=0; i<I; i++){
      // calculate X_sum
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      rowvec X_eta_alpha = X_bs_is*eta_alpha_s;
      X_sum += X_eta_alpha.t()*X_eta_alpha*tau(i,0);
      
      // calculate Xy_sum
      double Y_i_tilde = Y(i) - mu(0) - U(i)*lambda(0);
      rowvec A_bs_i = A_bs.submat(i, 0, i, B-1);
      Y_i_tilde -= conv_to<double>::from(A_bs_i*gamma);
      rowvec Z_i = Z.submat(i, 0, i, P-1);
      vec beta_q = conv_to<vec>::from(beta.row(0));
      Y_i_tilde -= conv_to<double>::from(Z_i*beta_q);
      for (int s2=0; s2<S_sum; s2++){
        rowvec X_bs_is2 = X_bs.subcube(i, s2, 0, i, s2, B-1);
        vec alpha_update_s2 = alpha_update.subcube(0, s2, 0, 0, s2, B-1);
        Y_i_tilde -= conv_to<double>::from(X_bs_is2*alpha_update_s2);
      }
      vec alpha_update_s = alpha_update.subcube(0, s, 0, 0, s, B-1);
      Y_i_tilde += conv_to<double>::from(X_bs_is*alpha_update_s);
      Xy_sum += conv_to<vec>::from(X_eta_alpha*tau(i,0)*Y_i_tilde);
    }
    
    // Gaussian posterior
    mat I = eye(B, B);
    mat V_n = inv_sympd(X_sum/sigma2(0) + I);
    vec m_alpha_s = m_alpha.subcube(0, s, 0, 0, s, B-1);
    vec mu_n = V_n*(Xy_sum/sigma2(0) + m_alpha_s);
    zeta_alpha_update.subcube(0, s, 0, 0, s, B-1) = rmvn_rcpp(1, mu_n, V_n);
    
    // update alpha
    alpha_update.subcube(0, s, 0, 0, s, B-1) = eta_alpha_s*zeta_alpha_update.subcube(0, s, 0, 0, s, B-1);
  }
  
  for (int s=0; s<S; s++){
    
    double eta_alpha_s = eta_alpha(1, s);
    
    mat X_sum(B, B); X_sum.fill(0);
    vec Xa_sum(B); Xa_sum.fill(0);
    for (int i=0; i<I; i++){
      // calculate X_sum
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      rowvec X_eta_alpha = X_bs_is*eta_alpha_s;
      X_sum += X_eta_alpha.t()*X_eta_alpha*tau(i,1);
      
      // calculate Xa_sum
      double A_i_tilde = A(i) - mu(1) - U(i)*lambda(1);
      rowvec Z_i = Z.submat(i, 0, i, P-1);
      vec beta_q = conv_to<vec>::from(beta.row(1));
      A_i_tilde -= conv_to<double>::from(Z_i*beta_q);
      for (int s2=0; s2<S; s2++){
        rowvec X_bs_is2 = X_bs.subcube(i, s2, 0, i, s2, B-1);
        vec alpha_update_s2 = alpha_update.subcube(1, s2, 0, 1, s2, B-1);
        A_i_tilde -= conv_to<double>::from(X_bs_is2*alpha_update_s2);
      }
      vec alpha_update_s = alpha_update.subcube(1, s, 0, 1, s, B-1);
      A_i_tilde += conv_to<double>::from(X_bs_is*alpha_update_s);
      Xa_sum += conv_to<vec>::from(X_eta_alpha*tau(i,1)*A_i_tilde);
    }
    
    // Gaussian posterior
    mat I = eye(B, B);
    mat V_n = inv_sympd(X_sum/sigma2(1) + I);
    vec m_alpha_s = m_alpha.subcube(1, s, 0, 1, s, B-1);
    vec mu_n = V_n*(Xa_sum/sigma2(1) + m_alpha_s);
    zeta_alpha_update.subcube(1, s, 0, 1, s, B-1) = rmvn_rcpp(1, mu_n, V_n);
    
    // update alpha
    alpha_update.subcube(1, s, 0, 1, s, B-1) = eta_alpha_s*zeta_alpha_update.subcube(1, s, 0, 1, s, B-1);
  }
  
  vec zero_vec(B); zero_vec.fill(0); 
  for (int s=S; s<S_sum; s++){
    zeta_alpha_update.subcube(1, s, 0, 1, s, B-1) = zero_vec;
  }
  
  return(zeta_alpha_update);
}


// [[Rcpp::export]]
mat update_beta_rcpp(mat& xi_beta, mat& nu_beta, vec& gamma, cube& alpha, 
                     vec& U, vec& lambda, vec& mu, mat& tau, vec& sigma2,
                     const int I, const int Q, const int S, const int S_sum, const int B, const int P,
                     vec& Y, vec& A, mat& A_bs, cube& X_bs, mat& Z){
  
  mat beta_update(Q, P);
  
  mat Z_sum(P, P); Z_sum.fill(0);
  vec Zy_sum(P); Zy_sum.fill(0);
  
  for (int i=0; i<I; i++){
    // calculate Z_sum
    rowvec Z_i = Z.row(i);
    Z_sum += Z_i.t()*Z_i*tau(i,0);
    
    // calculate Zy_sum
    double Y_i_tilde = Y(i) - mu(0) - U(i)*lambda(0);
    rowvec A_bs_i = A_bs.submat(i, 0, i, B-1);
    Y_i_tilde -= conv_to<double>::from(A_bs_i*gamma);
    for (int s=0; s<S_sum; s++){
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      vec alpha_s = alpha.subcube(0, s, 0, 0, s, B-1);
      Y_i_tilde -= conv_to<double>::from(X_bs_is*alpha_s);
    }
    Zy_sum += conv_to<vec>::from(Z_i*tau(i,0)*Y_i_tilde);
  }
  
  // Gaussian posterior
  mat Sigma_beta_inv(P, P); Sigma_beta_inv.fill(0);
  for (int p=0; p<P; p++){ Sigma_beta_inv(p,p) = 1/(xi_beta(0,p)*nu_beta(0,p)); }
  mat V_n = inv_sympd(Z_sum/sigma2(0) + Sigma_beta_inv);
  vec mu_n = V_n*(Zy_sum/sigma2(0));
  beta_update.row(0) = rmvn_rcpp(1, mu_n, V_n);
  
  Z_sum.fill(0);
  vec Za_sum(P); Za_sum.fill(0);
  
  for (int i=0; i<I; i++){
    // calculate Z_sum
    rowvec Z_i = Z.row(i);
    Z_sum += Z_i.t()*Z_i*tau(i,1);
    // calculate Za_sum
    double A_i_tilde = A(i) - mu(1) - U(i)*lambda(1);
    for (int s=0; s<S; s++){
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      vec alpha_s = alpha.subcube(1, s, 0, 1, s, B-1);
      A_i_tilde -= conv_to<double>::from(X_bs_is*alpha_s);
    }
    Za_sum += conv_to<vec>::from(Z_i*tau(i,1)*A_i_tilde);
  }
  
  // Gaussian posterior
  Sigma_beta_inv.fill(0);
  for (int p=0; p<P; p++){ Sigma_beta_inv(p,p) = 1/(xi_beta(1,p)*nu_beta(1,p)); }
  V_n = inv_sympd(Z_sum/sigma2(1) + Sigma_beta_inv);
  mu_n = V_n*(Za_sum/sigma2(1));
  beta_update.row(1) = rmvn_rcpp(1, mu_n, V_n);
  
  return(beta_update);
}


// [[Rcpp::export]]
vec update_U_rcpp(vec& gamma, cube& alpha, mat& beta, vec& lambda, vec& mu, mat& tau, vec& sigma2, 
                  const int I, const int Q, const int S, const int S_sum, const int B, const int P,
                  vec& Y, vec& A, mat& A_bs, cube& X_bs, mat& Z){
  
  vec U_update(I);
  
  for (int i=0; i<I; i++){
    
    // calculate YA_tilde_i
    vec YA_tilde_i(Q);
    YA_tilde_i(0) = Y(i) - mu(0);
    rowvec A_bs_i = A_bs.submat(i, 0, i, B-1);
    YA_tilde_i(0) -= conv_to<double>::from(A_bs_i*gamma);
    rowvec Z_i = Z.submat(i, 0, i, P-1);
    vec beta_q = conv_to<vec>::from(beta.row(0));
    YA_tilde_i(0) -= conv_to<double>::from(Z_i*beta_q);
    for (int s=0; s<S_sum; s++){
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      vec alpha_s = alpha.subcube(0, s, 0, 0, s, B-1);
      YA_tilde_i(0) -= conv_to<double>::from(X_bs_is*alpha_s);
    }
    
    YA_tilde_i(1) = A(i) - mu(1);
    Z_i = Z.submat(i, 0, i, P-1);
    beta_q = conv_to<vec>::from(beta.row(1));
    YA_tilde_i(1) -= conv_to<double>::from(Z_i*beta_q);
    for (int s=0; s<S_sum; s++){
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      vec alpha_s = alpha.subcube(1, s, 0, 1, s, B-1);
      YA_tilde_i(1) -= conv_to<double>::from(X_bs_is*alpha_s);
    }
    
    // calculate mu_u and V_u
    mat Sigma_i_inv(Q,Q); Sigma_i_inv.fill(0);
    for (int q=0; q<Q; q++){ Sigma_i_inv(q,q) = tau(i,q)/sigma2(q); }
    double V_u = 1/(conv_to<double>::from(lambda.t()*Sigma_i_inv*lambda) + 1);
    double mu_u = V_u*(conv_to<double>::from(lambda.t()*Sigma_i_inv*YA_tilde_i));
    
    // Gaussian posterior distribution
    U_update(i) = R::rnorm(mu_u, sqrt(V_u));
  }
  
  return(U_update);
}


// [[Rcpp::export]]
vec update_lambda_rcpp(vec& gamma, cube& alpha, mat& beta, vec& U, vec& mu, mat& tau, vec& sigma2, 
                       const int I, const int Q, const int S, const int S_sum, const int B, const int P,
                       vec& Y, vec& A, mat& A_bs, cube& X_bs, mat& Z){
  
  vec lambda_update(Q);
  
  double U_sum = 0; double Uy_sum = 0;
  
  for (int i=0; i<I; i++){
    // calculate U_sum
    U_sum += pow(U(i), 2)*tau(i,0);
    
    // calculate Uy_sum
    double Y_i_tilde = Y(i) - mu(0);
    rowvec A_bs_i = A_bs.submat(i, 0, i, B-1);
    Y_i_tilde -= conv_to<double>::from(A_bs_i*gamma);
    rowvec Z_i = Z.submat(i, 0, i, P-1);
    vec beta_q = conv_to<vec>::from(beta.row(0));
    Y_i_tilde -= conv_to<double>::from(Z_i*beta_q);
    for (int s=0; s<S_sum; s++){
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      vec alpha_s = alpha.subcube(0, s, 0, 0, s, B-1);
      Y_i_tilde -= conv_to<double>::from(X_bs_is*alpha_s);
    }
    Uy_sum += U(i)*tau(i,0)*Y_i_tilde;
  }
  
  // Gaussian posterior
  double V_m = 1/(U_sum/sigma2(0) + 1/100);
  double mu_m = V_m*Uy_sum/sigma2(0);
  lambda_update(0) = R::rnorm(mu_m, sqrt(V_m));
  
  U_sum = 0; double Ua_sum = 0;
  
  for (int i=0; i<I; i++){
    // calculate U_sum
    U_sum += pow(U(i), 2)*tau(i,1);
    
    // calculate Ua_sum
    double A_i_tilde = A(i) - mu(1);
    rowvec Z_i = Z.submat(i, 0, i, P-1);
    vec beta_q = conv_to<vec>::from(beta.row(1));
    A_i_tilde -= conv_to<double>::from(Z_i*beta_q);
    for (int s=0; s<S; s++){
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      vec alpha_s = alpha.subcube(1, s, 0, 1, s, B-1);
      A_i_tilde -= conv_to<double>::from(X_bs_is*alpha_s);
    }
    Ua_sum += U(i)*tau(i,1)*A_i_tilde;
  }
  
  // Gaussian posterior
  double V_n = 1/(U_sum/sigma2(1) + 1/100);
  double mu_n = V_n*Ua_sum/sigma2(1);
  lambda_update(1) = R::rnorm(mu_n, sqrt(V_n));
  
  return(lambda_update);
}


// [[Rcpp::export]]
mat update_tau_rcpp(vec& gamma, cube& alpha, mat& beta, vec& U, vec& lambda, vec& mu, vec& sigma2,
                    const int I, const int Q, const int S, const int S_sum, const int B, const int P,
                    vec& Y, vec& A, mat& A_bs, cube& X_bs, mat& Z){
  
  mat tau_update(I,Q);
  double eps = 1e-3; // avoid numerical issues
  
  for (int i=0; i<I; i++){
    // calculate Y_tilde_i
    double Y_tilde_i = Y(i) - mu(0) - U(i)*lambda(0);
    rowvec A_bs_i = A_bs.submat(i, 0, i, B-1);
    Y_tilde_i -= conv_to<double>::from(A_bs_i*gamma);
    rowvec Z_i = Z.submat(i, 0, i, P-1);
    vec beta_q = conv_to<vec>::from(beta.row(0));
    Y_tilde_i -= conv_to<double>::from(Z_i*beta_q);
    for (int s=0; s<S_sum; s++){
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      vec alpha_s = alpha.subcube(0, s, 0, 0, s, B-1);
      Y_tilde_i -= conv_to<double>::from(X_bs_is*alpha_s);
    }
    // Inverse-Gaussian posterior distribution
    double mean = sqrt(sigma2(0))/(2*abs(Y_tilde_i));
    tau_update(i,0) = max(rinvgaussian_rcpp(mean, 0.25), eps);
    
    // calculate A_tilde_i
    double A_tilde_i = A(i) - mu(1) - U(i)*lambda(1);
    Z_i = Z.submat(i, 0, i, P-1);
    beta_q = conv_to<vec>::from(beta.row(1));
    A_tilde_i -= conv_to<double>::from(Z_i*beta_q);
    for (int s=0; s<S; s++){
      rowvec X_bs_is = X_bs.subcube(i, s, 0, i, s, B-1);
      vec alpha_s = alpha.subcube(1, s, 0, 1, s, B-1);
      A_tilde_i -= conv_to<double>::from(X_bs_is*alpha_s);
    }
    // Inverse-Gaussian posterior distribution
    mean = sqrt(sigma2(1))/(2*abs(A_tilde_i));
    tau_update(i,1) = max(rinvgaussian_rcpp(mean, 0.25), eps);
  }
  
  return(tau_update);
}