//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp; 


const double log2pi = log(2.0*M_PI);


// [[Rcpp::export]]
double dmvn_rcpp(rowvec& x, rowvec& mean, mat& sigma, bool logd = false){ 
  
  // Calculate density of the multivariate normal distribution.
  // Input: 1. x (row vector, dim=c(xdim)): row vector data;
  //        2. mean (row vector, dim=c(xdim)): mean of multivariate normal distribution;
  //        3. sigma (matrix, dim=c(xdim,xdim)): covariance matrix of multivariate normal distribution;
  //        4. logd (scalar, boolean): true for taking log.
  // Output: out (scalar, double): pdf (or log pdf) of multivariate normal distribution.
  
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
  
  // Randomly generate samples from the multivariate normal distribution.
  // Input: 1. n (scalar, integer): number of data points;
  //        2. mean (row vector, dim=c(k)): mean of multivariate normal distribution;
  //        3. sigma (matrix, dim=c(k,k)): covariance matrix of multivariate normal distribution.
  // Output: out (matrix, dim=c(n,k)): n random samples from the k-dimensional multivariate normal distribution.
  
  int k = sigma.n_cols; // dimension of the multivariate normal distribution
  mat z = randn(n, k);
  mat out = repmat(mean,1,n).t()+z*chol(sigma);
  return(out);
}


// [[Rcpp::export]]
double rinvgamma_rcpp(const double a, const double b){
  
  // Generate a random sample from the Inverse-Gamma distribution.
  // Input: 1. a (scalar, double): the first parameter in Inverse-Gamma(a, b);
  //        2. b (scalar, double): the second parameter in Inverse-Gamma(a, b).
  // Output: a random sample from Inverse-Gamma(a,b).
  
  return(1/R::rgamma(a, 1/b));
}


// [[Rcpp::export]]
double rinvgaussian_rcpp(const double mu, const double lambda){
  
  // Generate a random sample from the Inverse-Gaussian distribution.
  // Input: 1. mu (scalar, double): the first parameter in Inverse-Gaussian(mu, lambda);
  //        2. lambda (scalar, double): the second parameter in Inverse-Gaussian(mu, lambda).
  // Output: out (scalar, double): a random sample from Inverse-Gaussian(mu, lambda).
  
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
    
  // Update the parameter eta_gamma in the MCMC algorithm.
  // Input: 1. zeta_gamma (vector, dim=c(B)): the spike-and-slab prior hyper-parameter that represents the magnitude of gamma if it is selected as a signal;
  //        2. xi_gamma (scalar, double): the hyper-parameter in the prior for eta_gamma, where eta_gamma ~ N(0, xi_gamma*nu_gamma), and xi_gamma ~ rho_gamma*Delta_{1}(xi_gamma) + (1-rho_gamma)*Delta_{nu_0}(xi_gamma), nu_0=2.5e-4;
  //        3. nu_gamma (scalar, double): the hyper-parameter in the prior for eta_gamma, where eta_gamma ~ N(0, xi_gamma*nu_gamma), and nu_gamma ~ Inverse-Gamma(a_nu, b_nu), a_nu=5, b_nu=50;
  //        4. alpha (cube, dim=c(Q=2, S_sum, B)): the estimated effects of the continuous covariates X, and the interactions of both the continuous and discrete covariates X and Z with the continuous treatment A (i.e., A*X and A*Z) on the outcome Y (q=1) and the treatment A (q=2) using cubic B-spline expansions, where S_sum=S*2+P, S and P denote the dimensions of X and Z, respectively;
  //        5. beta (matrix, dim=c(Q=2, P)): the estimated effects of the discrete covariates Z on Y (q=1) and A (q=2);
  //        6. U (vector, dim=c(I)): the unmeasured confounder vector, where U[i] denotes the value for individual i, i=1,2,...,I, and I denotes the total number of individuals;
  //        7. lambda (vector, dim=c(Q=2)): the effects of unmeasured confounder U on the outcome Y (q=1) and the treatment A (q=2);
  //        8. mu (vector, dim=c(Q=2)): the global intercept terms for the outcome Y (q=1) and the treatment A (q=2);
  //        9. tau (matrix, dim=c(I, Q=2)): the auxiliary variable associated with the Laplace error term for Y (q=1) and A (q=2), where tau[i,1] and tau[i,2] ~ Inverse-Gamma(1, 1/8) for individual i, i=1,2,...,I;
  //        10. sigma2 (vector, dim=c(Q=2)): the variance of the Laplace error term for Y (q=1) and A (q=2), where sigma2 ~ Inverse-Gamma(a_sigma, b_sigma), a_sigma=b_sigma=1;
  //        11. I (scalar, integer): the total number of individuals;
  //        12. S_sum (scalar, integer): S_sum=S*2+P, where S and P denote the dimensions of continuous and discrete covariates X and Z, respectively;
  //        13. B (scalar, integer): the degrees of freedom B=B0-1 for the cubic B-spline expansion after imposing the sum-to-zero constraint;
  //        14. P (scalar, integer): the dimension of the discrete covariates Z;
  //        15. Y (vector, dim=c(I)): the outcome vector, where Y[i] denotes the value for individual i, i=1,2,...,I;
  //        16. A_bs (matrix, dim=c(I, B)): the cubic B-spline expansion for the treatment A, where A_bs[i,] denotes the vector value for individual i, i=1,2,...,I;
  //        17. X_bs (cube, dim=c(I, S_sum, B)): the cubic B-spline expansion for the continuous covariates and interactions, where X_bs[i,s,] denotes the vector value for individual i and covariate s, i=1,2,...,I, s=1,2,...,S_sum;
  //        18. Z (matrix, dim=c(I, P)): the discrete covariates Z, where Z[i,p] denotes the value for individual i and discrete covariate p, i=1,2,...,I, p=1,2,...,P.
  // Output: eta_gamma_update (scalar, double): the updated value of eta_gamma for the current MCMC iteration.
  
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
    
  // Update the parameter zeta_gamma in the MCMC algorithm.
  // Input: 1. eta_gamma (scalar, double): the spike-and-slab prior hyper-parameter that indicates whether gamma is selected as a signal (i.e., non-zero values) or not, where gamma=eta_gamma*zeta_gamma;
  //        2. m_gamma (vector, dim=c(B)): the hyper-parameter in the prior for zeta_gamma, where zeta_gamma[b] ~ N(m_gamma[b], 1), for b=1,2,...,B, and m_gamma[b] ~ 0.5*Delta_{1}(m_gamma[b]) + 0.5*Delta_{-1}(m_gamma[b]);
  //        3. alpha (cube, dim=c(Q=2, S_sum, B)): the estimated effects of the continuous covariates X, and the interactions of both the continuous and discrete covariates X and Z with the continuous treatment A (i.e., A*X and A*Z) on the outcome Y (q=1) and the treatment A (q=2) using cubic B-spline expansions, where S_sum=S*2+P, S and P denote the dimensions of X and Z, respectively;
  //        4. beta (matrix, dim=c(Q=2, P)): the estimated effects of the discrete covariates Z on Y (q=1) and A (q=2);
  //        5. U (vector, dim=c(I)): the unmeasured confounder vector, where U[i] denotes the value for individual i, i=1,2,...,I, and I denotes the total number of individuals;
  //        6. lambda (vector, dim=c(Q=2)): the effects of unmeasured confounder U on the outcome Y (q=1) and the treatment A (q=2);
  //        7. mu (vector, dim=c(Q=2)): the global intercept terms for the outcome Y (q=1) and the treatment A (q=2);
  //        8. tau (matrix, dim=c(I, Q=2)): the auxiliary variable associated with the Laplace error term for Y (q=1) and A (q=2), where tau[i,1] and tau[i,2] ~ Inverse-Gamma(1, 1/8) for individual i, i=1,2,...,I;
  //        9. sigma2 (vector, dim=c(Q=2)): the variance of the Laplace error term for Y (q=1) and A (q=2), where sigma2 ~ Inverse-Gamma(a_sigma, b_sigma), a_sigma=b_sigma=1;
  //        10. I (scalar, integer): the total number of individuals;
  //        11. S_sum (scalar, integer): S_sum=S*2+P, where S and P denote the dimensions of continuous and discrete covariates X and Z, respectively;
  //        12. B (scalar, integer): the degrees of freedom B=B0-1 for the cubic B-spline expansion after imposing the sum-to-zero constraint;
  //        13. P (scalar, integer): the dimension of the discrete covariates Z;
  //        14. Y (vector, dim=c(I)): the outcome vector, where Y[i] denotes the value for individual i, i=1,2,...,I;
  //        15. A_bs (matrix, dim=c(I, B)): the cubic B-spline expansion for the treatment A, where A_bs[i,] denotes the vector value for individual i, i=1,2,...,I;
  //        16. X_bs (cube, dim=c(I, S_sum, B)): the cubic B-spline expansion for the continuous covariates and interactions, where X_bs[i,s,] denotes the vector value for individual i and covariate s, i=1,2,...,I, s=1,2,...,S_sum;
  //        17. Z (matrix, dim=c(I, P)): the discrete covariates Z, where Z[i,p] denotes the value for individual i and discrete covariate p, i=1,2,...,I, p=1,2,...,P.
  // Output: zeta_gamma_update (vector, dim=c(B)): the updated value of zeta_gamma for the current MCMC iteration.
  
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
  
  // Update the parameter eta_alpha in the MCMC algorithm.
  // Input: 1. gamma (vector, dim=c(B)): the estimated treatment effect of A on Y using cubic B-spline expansion, where the degrees of freedom B=B0-1 after imposing the sum-to-zero constraint;
  //        2. eta_alpha (matrix, dim=c(Q=2, S_sum)): eta_alpha[q,s] indicates whether alpha[q,s,] is selected as a signal or not, where alpha[q,s,]=eta_alpha[q,s]*zeta_alpha[q,s,], for q=1,2 and s=1,2,...,S_sum;
  //        3. zeta_alpha (cube, dim=c(Q=2, S_sum, B)): zeta_alpha[q,s,] represents the magnitude of alpha[q,s,] if it is selected as a signal;
  //        4. xi_alpha (matrix, dim=c(Q=2, S_sum)): the hyper-parameter in the prior for eta_alpha, where eta_alpha[q,s] ~ N(0, xi_alpha[q,s]*nu_alpha[q,s]);
  //        5. nu_alpha (matrix, dim=c(Q=2, S_sum)): the hyper-parameter in the prior for eta_alpha, where eta_alpha[q,s] ~ N(0, xi_alpha[q,s]*nu_alpha[q,s]), and nu_alpha[q,s] ~ Inverse-Gamma(a_nu, b_nu);
  //        6. beta (matrix, dim=c(Q=2, P)): the estimated effects of the discrete covariates Z on Y (q=1) and A (q=2);
  //        7. U (vector, dim=c(I)): the unmeasured confounder vector, where U[i] denotes the value for individual i, i=1,2,...,I, and I denotes the total number of individuals;
  //        8. lambda (vector, dim=c(Q=2)): the effects of unmeasured confounder U on the outcome Y (q=1) and the treatment A (q=2);
  //        9. mu (vector, dim=c(Q=2)): the global intercept terms for the outcome Y (q=1) and the treatment A (q=2);
  //        10. tau (matrix, dim=c(I, Q=2)): the auxiliary variable associated with the Laplace error term for Y (q=1) and A (q=2), where tau[i,1] and tau[i,2] ~ Inverse-Gamma(1, 1/8) for individual i, i=1,2,...,I;
  //        11. sigma2 (vector, dim=c(Q=2)): the variance of the Laplace error term for Y (q=1) and A (q=2), where sigma2 ~ Inverse-Gamma(a_sigma, b_sigma), a_sigma=b_sigma=1;
  //        12. I (scalar, integer): the total number of individuals;
  //        13. Q (scalar, integer): the number of outcome and treatment, where Q=2;
  //        14. S (scalar, integer): the dimension of continuous covariates X;
  //        15. S_sum (scalar, integer): S_sum=S*2+P, where S and P denote the dimensions of continuous and discrete covariates X and Z, respectively;
  //        16. B (scalar, integer): the degrees of freedom B=B0-1 for the cubic B-spline expansion after imposing the sum-to-zero constraint;
  //        17. P (scalar, integer): the dimension of the discrete covariates Z;
  //        18. Y (vector, dim=c(I)): the outcome vector, where Y[i] denotes the value for individual i, i=1,2,...,I;
  //        19. A (vector, dim=c(I)): the treatment vector, where A[i] denotes the value for individual i, i=1,2,...,I;
  //        20. A_bs (matrix, dim=c(I, B)): the cubic B-spline expansion for the treatment A, where A_bs[i,] denotes the vector value for individual i, i=1,2,...,I;
  //        21. X_bs (cube, dim=c(I, S_sum, B)): the cubic B-spline expansion for the continuous covariates and interactions, where X_bs[i,s,] denotes the vector value for individual i and covariate s, i=1,2,...,I, s=1,2,...,S_sum;
  //        22. Z (matrix, dim=c(I, P)): the discrete covariates Z, where Z[i,p] denotes the value for individual i and discrete covariate p, i=1,2,...,I, p=1,2,...,P.
  // Output: eta_alpha_update (matrix, dim=c(Q=2, S_sum)): the updated value of eta_alpha for the current MCMC iteration.
  
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
    
  // Update the parameter zeta_alpha in the MCMC algorithm.
  // Input: 1. gamma (vector, dim=c(B)): the estimated treatment effect of A on Y using cubic B-spline expansion, where the degrees of freedom B=B0-1 after imposing the sum-to-zero constraint;
  //        2. eta_alpha (matrix, dim=c(Q=2, S_sum)): eta_alpha[q,s] indicates whether alpha[q,s,] is selected as a signal or not, where alpha[q,s,]=eta_alpha[q,s]*zeta_alpha[q,s,], for q=1,2 and s=1,2,...,S_sum;
  //        3. zeta_alpha (cube, dim=c(Q=2, S_sum, B)): zeta_alpha[q,s,] represents the magnitude of alpha[q,s,] if it is selected as a signal;
  //        4. m_alpha (cube, dim=c(Q=2, S_sum, B)): the hyper-parameter in the prior for zeta_alpha, where zeta_alpha[q,s,b] ~ N(m_alpha[q,s,b], 1), for q=1,2, s=1,2,...,S_sum, and b=1,2,...,B, and m_alpha[q,s,b] ~ 0.5*Delta_{1}(m_alpha[q,s,b]) + 0.5*Delta_{-1}(m_alpha[q,s,b]);
  //        5. beta (matrix, dim=c(Q=2, P)): the estimated effects of the discrete covariates Z on Y (q=1) and A (q=2);
  //        6. U (vector, dim=c(I)): the unmeasured confounder vector, where U[i] denotes the value for individual i, i=1,2,...,I, and I denotes the total number of individuals;
  //        7. lambda (vector, dim=c(Q=2)): the effects of unmeasured confounder U on the outcome Y (q=1) and the treatment A (q=2);
  //        8. mu (vector, dim=c(Q=2)): the global intercept terms for the outcome Y (q=1) and the treatment A (q=2);
  //        9. tau (matrix, dim=c(I, Q=2)): the auxiliary variable associated with the Laplace error term for Y (q=1) and A (q=2), where tau[i,1] and tau[i,2] ~ Inverse-Gamma(1, 1/8) for individual i, i=1,2,...,I;
  //        10. sigma2 (vector, dim=c(Q=2)): the variance of the Laplace error term for Y (q=1) and A (q=2), where sigma2 ~ Inverse-Gamma(a_sigma, b_sigma), a_sigma=b_sigma=1;
  //        11. I (scalar, integer): the total number of individuals;
  //        12. Q (scalar, integer): the number of outcome and treatment, where Q=2;
  //        13. S (scalar, integer): the dimension of continuous covariates X;
  //        14. S_sum (scalar, integer): S_sum=S*2+P, where S and P denote the dimensions of continuous and discrete covariates X and Z, respectively;
  //        15. B (scalar, integer): the degrees of freedom B=B0-1 for the cubic B-spline expansion after imposing the sum-to-zero constraint;
  //        16. P (scalar, integer): the dimension of the discrete covariates Z;
  //        17. Y (vector, dim=c(I)): the outcome vector, where Y[i] denotes the value for individual i, i=1,2,...,I;
  //        18. A (vector, dim=c(I)): the treatment vector, where A[i] denotes the value for individual i, i=1,2,...,I;
  //        19. A_bs (matrix, dim=c(I, B)): the cubic B-spline expansion for the treatment A, where A_bs[i,] denotes the vector value for individual i, i=1,2,...,I;
  //        20. X_bs (cube, dim=c(I, S_sum, B)): the cubic B-spline expansion for the continuous covariates and interactions, where X_bs[i,s,] denotes the vector value for individual i and covariate s, i=1,2,...,I, s=1,2,...,S_sum;
  //        21. Z (matrix, dim=c(I, P)): the discrete covariates Z, where Z[i,p] denotes the value for individual i and discrete covariate p, i=1,2,...,I, p=1,2,...,P.
  // Output: zeta_alpha_update (cube, dim=c(Q=2, S_sum, B)): the updated value of zeta_alpha for the current MCMC iteration.
  
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
    
  // Update the parameter beta in the MCMC algorithm.
  // Input: 1. xi_beta (matrix, dim=c(Q=2, P)): the hyper-parameter in the prior for beta, where beta[q,p] ~ N(0, xi_beta[q,p]*nu_beta[q,p]), and xi_beta[q,p] ~ rho_beta*Delta_{1}(xi_beta[q,p]) + (1-rho_beta)*Delta_{nu_0}(xi_beta[q,p]);
  //        2. nu_beta (matrix, dim=c(Q=2, P)): the hyper-parameter in the prior for beta, where beta[q,p] ~ N(0, xi_beta[q,p]*nu_beta[q,p]), and nu_beta[q,p] ~ Inverse-Gamma(a_nu, b_nu), for q=1,2, and p=1,2,...,P;
  //        3. gamma (vector, dim=c(B)): the estimated treatment effect of A on Y using cubic B-spline expansion, where the degrees of freedom B=B0-1 after imposing the sum-to-zero constraint;
  //        4. alpha (cube, dim=c(Q=2, S_sum, B)): the estimated effects of the continuous covariates X, and the interactions of both the continuous and discrete covariates X and Z with the continuous treatment A (i.e., A*X and A*Z) on the outcome Y (q=1) and the treatment A (q=2) using cubic B-spline expansions, where S_sum=S*2+P, S and P denote the dimensions of X and Z, respectively;
  //        5. U (vector, dim=c(I)): the unmeasured confounder vector, where U[i] denotes the value for individual i, i=1,2,...,I, and I denotes the total number of individuals;
  //        6. lambda (vector, dim=c(Q=2)): the effects of unmeasured confounder U on the outcome Y (q=1) and the treatment A (q=2);
  //        7. mu (vector, dim=c(Q=2)): the global intercept terms for the outcome Y (q=1) and the treatment A (q=2);
  //        8. tau (matrix, dim=c(I, Q=2)): the auxiliary variable associated with the Laplace error term for Y (q=1) and A (q=2), where tau[i,1] and tau[i,2] ~ Inverse-Gamma(1, 1/8) for individual i, i=1,2,...,I;
  //        9. sigma2 (vector, dim=c(Q=2)): the variance of the Laplace error term for Y (q=1) and A (q=2), where sigma2 ~ Inverse-Gamma(a_sigma, b_sigma), a_sigma=b_sigma=1;
  //        10. I (scalar, integer): the total number of individuals;
  //        11. Q (scalar, integer): the number of outcome and treatment, where Q=2;
  //        12. S (scalar, integer): the dimension of continuous covariates X;
  //        13. S_sum (scalar, integer): S_sum=S*2+P, where S and P denote the dimensions of continuous and discrete covariates X and Z, respectively;
  //        14. B (scalar, integer): the degrees of freedom B=B0-1 for the cubic B-spline expansion after imposing the sum-to-zero constraint;
  //        15. P (scalar, integer): the dimension of the discrete covariates Z;
  //        16. Y (vector, dim=c(I)): the outcome vector, where Y[i] denotes the value for individual i, i=1,2,...,I;
  //        17. A (vector, dim=c(I)): the treatment vector, where A[i] denotes the value for individual i, i=1,2,...,I;
  //        18. A_bs (matrix, dim=c(I, B)): the cubic B-spline expansion for the treatment A, where A_bs[i,] denotes the vector value for individual i, i=1,2,...,I;
  //        19. X_bs (cube, dim=c(I, S_sum, B)): the cubic B-spline expansion for the continuous covariates and interactions, where X_bs[i,s,] denotes the vector value for individual i and covariate s, i=1,2,...,I, s=1,2,...,S_sum;
  //        20. Z (matrix, dim=c(I, P)): the discrete covariates Z, where Z[i,p] denotes the value for individual i and discrete covariate p, i=1,2,...,I, p=1,2,...,P.
  // Output: beta_update (matrix, dim=c(Q=2, P)): the updated value of beta for the current MCMC iteration.
  
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
    
  // Update the parameter U in the MCMC algorithm.
  // Input: 1. gamma (vector, dim=c(B)): the estimated treatment effect of A on Y using cubic B-spline expansion, where the degrees of freedom B=B0-1 after imposing the sum-to-zero constraint;
  //        2. alpha (cube, dim=c(Q=2, S_sum, B)): the estimated effects of the continuous covariates X, and the interactions of both the continuous and discrete covariates X and Z with the continuous treatment A (i.e., A*X and A*Z) on the outcome Y (q=1) and the treatment A (q=2) using cubic B-spline expansions, where S_sum=S*2+P, S and P denote the dimensions of X and Z, respectively;
  //        3. beta (matrix, dim=c(Q=2, P)): the estimated effects of the discrete covariates Z on Y (q=1) and A (q=2);
  //        4. lambda (vector, dim=c(Q=2)): the effects of unmeasured confounder U on the outcome Y (q=1) and the treatment A (q=2);
  //        5. mu (vector, dim=c(Q=2)): the global intercept terms for the outcome Y (q=1) and the treatment A (q=2);
  //        6. tau (matrix, dim=c(I, Q=2)): the auxiliary variable associated with the Laplace error term for Y (q=1) and A (q=2), where tau[i,1] and tau[i,2] ~ Inverse-Gamma(1, 1/8) for individual i, i=1,2,...,I;
  //        7. sigma2 (vector, dim=c(Q=2)): the variance of the Laplace error term for Y (q=1) and A (q=2), where sigma2 ~ Inverse-Gamma(a_sigma, b_sigma), a_sigma=b_sigma=1;
  //        8. I (scalar, integer): the total number of individuals;
  //        9. Q (scalar, integer): the number of outcome and treatment, where Q=2;
  //        10. S (scalar, integer): the dimension of continuous covariates X;
  //        11. S_sum (scalar, integer): S_sum=S*2+P, where S and P denote the dimensions of continuous and discrete covariates X and Z, respectively;
  //        12. B (scalar, integer): the degrees of freedom B=B0-1 for the cubic B-spline expansion after imposing the sum-to-zero constraint;
  //        13. P (scalar, integer): the dimension of the discrete covariates Z;
  //        14. Y (vector, dim=c(I)): the outcome vector, where Y[i] denotes the value for individual i, i=1,2,...,I;
  //        15. A (vector, dim=c(I)): the treatment vector, where A[i] denotes the value for individual i, i=1,2,...,I;
  //        16. A_bs (matrix, dim=c(I, B)): the cubic B-spline expansion for the treatment A, where A_bs[i,] denotes the vector value for individual i, i=1,2,...,I;
  //        17. X_bs (cube, dim=c(I, S_sum, B)): the cubic B-spline expansion for the continuous covariates and interactions, where X_bs[i,s,] denotes the vector value for individual i and covariate s, i=1,2,...,I, s=1,2,...,S_sum;
  //        18. Z (matrix, dim=c(I, P)): the discrete covariates Z, where Z[i,p] denotes the value for individual i and discrete covariate p, i=1,2,...,I, p=1,2,...,P.
  // Output: U_update (vector, dim=c(I)): the updated value of U for the current MCMC iteration.
  
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
    
  // Update the parameter lambda in the MCMC algorithm.
  // Input: 1. gamma (vector, dim=c(B)): the estimated treatment effect of A on Y using cubic B-spline expansion, where the degrees of freedom B=B0-1 after imposing the sum-to-zero constraint;
  //        2. alpha (cube, dim=c(Q=2, S_sum, B)): the estimated effects of the continuous covariates X, and the interactions of both the continuous and discrete covariates X and Z with the continuous treatment A (i.e., A*X and A*Z) on the outcome Y (q=1) and the treatment A (q=2) using cubic B-spline expansions, where S_sum=S*2+P, S and P denote the dimensions of X and Z, respectively;
  //        3. beta (matrix, dim=c(Q=2, P)): the estimated effects of the discrete covariates Z on Y (q=1) and A (q=2);
  //        4. U (vector, dim=c(I)): the unmeasured confounder vector, where U[i] denotes the value for individual i, i=1,2,...,I, and I denotes the total number of individuals;
  //        5. mu (vector, dim=c(Q=2)): the global intercept terms for the outcome Y (q=1) and the treatment A (q=2);
  //        6. tau (matrix, dim=c(I, Q=2)): the auxiliary variable associated with the Laplace error term for Y (q=1) and A (q=2), where tau[i,1] and tau[i,2] ~ Inverse-Gamma(1, 1/8) for individual i, i=1,2,...,I;
  //        7. sigma2 (vector, dim=c(Q=2)): the variance of the Laplace error term for Y (q=1) and A (q=2), where sigma2 ~ Inverse-Gamma(a_sigma, b_sigma), a_sigma=b_sigma=1;
  //        8. I (scalar, integer): the total number of individuals;
  //        9. Q (scalar, integer): the number of outcome and treatment, where Q=2;
  //        10. S (scalar, integer): the dimension of continuous covariates X;
  //        11. S_sum (scalar, integer): S_sum=S*2+P, where S and P denote the dimensions of continuous and discrete covariates X and Z, respectively;
  //        12. B (scalar, integer): the degrees of freedom B=B0-1 for the cubic B-spline expansion after imposing the sum-to-zero constraint;
  //        13. P (scalar, integer): the dimension of the discrete covariates Z;
  //        14. Y (vector, dim=c(I)): the outcome vector, where Y[i] denotes the value for individual i, i=1,2,...,I;
  //        15. A (vector, dim=c(I)): the treatment vector, where A[i] denotes the value for individual i, i=1,2,...,I;
  //        16. A_bs (matrix, dim=c(I, B)): the cubic B-spline expansion for the treatment A, where A_bs[i,] denotes the vector value for individual i, i=1,2,...,I;
  //        17. X_bs (cube, dim=c(I, S_sum, B)): the cubic B-spline expansion for the continuous covariates and interactions, where X_bs[i,s,] denotes the vector value for individual i and covariate s, i=1,2,...,I, s=1,2,...,S_sum;
  //        18. Z (matrix, dim=c(I, P)): the discrete covariates Z, where Z[i,p] denotes the value for individual i and discrete covariate p, i=1,2,...,I, p=1,2,...,P.
  // Output: lambda_update (vector, dim=c(I)): the updated value of lambda for the current MCMC iteration.
  
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
    
  // Update the parameter tau in the MCMC algorithm.
  // Input: 1. gamma (vector, dim=c(B)): the estimated treatment effect of A on Y using cubic B-spline expansion, where the degrees of freedom B=B0-1 after imposing the sum-to-zero constraint;
  //        2. alpha (cube, dim=c(Q=2, S_sum, B)): the estimated effects of the continuous covariates X, and the interactions of both the continuous and discrete covariates X and Z with the continuous treatment A (i.e., A*X and A*Z) on the outcome Y (q=1) and the treatment A (q=2) using cubic B-spline expansions, where S_sum=S*2+P, S and P denote the dimensions of X and Z, respectively;
  //        3. beta (matrix, dim=c(Q=2, P)): the estimated effects of the discrete covariates Z on Y (q=1) and A (q=2);
  //        4. U (vector, dim=c(I)): the unmeasured confounder vector, where U[i] denotes the value for individual i, i=1,2,...,I, and I denotes the total number of individuals;
  //        5. lambda (vector, dim=c(Q=2)): the effects of unmeasured confounder U on the outcome Y (q=1) and the treatment A (q=2);
  //        6. mu (vector, dim=c(Q=2)): the global intercept terms for the outcome Y (q=1) and the treatment A (q=2);
  //        7. sigma2 (vector, dim=c(Q=2)): the variance of the Laplace error term for Y (q=1) and A (q=2), where sigma2 ~ Inverse-Gamma(a_sigma, b_sigma), a_sigma=b_sigma=1;
  //        8. I (scalar, integer): the total number of individuals;
  //        9. Q (scalar, integer): the number of outcome and treatment, where Q=2;
  //        10. S (scalar, integer): the dimension of continuous covariates X;
  //        11. S_sum (scalar, integer): S_sum=S*2+P, where S and P denote the dimensions of continuous and discrete covariates X and Z, respectively;
  //        12. B (scalar, integer): the degrees of freedom B=B0-1 for the cubic B-spline expansion after imposing the sum-to-zero constraint;
  //        13. P (scalar, integer): the dimension of the discrete covariates Z;
  //        14. Y (vector, dim=c(I)): the outcome vector, where Y[i] denotes the value for individual i, i=1,2,...,I;
  //        15. A (vector, dim=c(I)): the treatment vector, where A[i] denotes the value for individual i, i=1,2,...,I;
  //        16. A_bs (matrix, dim=c(I, B)): the cubic B-spline expansion for the treatment A, where A_bs[i,] denotes the vector value for individual i, i=1,2,...,I;
  //        17. X_bs (cube, dim=c(I, S_sum, B)): the cubic B-spline expansion for the continuous covariates and interactions, where X_bs[i,s,] denotes the vector value for individual i and covariate s, i=1,2,...,I, s=1,2,...,S_sum;
  //        18. Z (matrix, dim=c(I, P)): the discrete covariates Z, where Z[i,p] denotes the value for individual i and discrete covariate p, i=1,2,...,I, p=1,2,...,P.
  // Output: tau_update (matrix, dim=c(I, Q=2)): the updated value of tau for the current MCMC iteration.
    
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
