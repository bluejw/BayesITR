####################################################### MCMC ########################################################

Nit <- 10000
burn.in <- 5000

mcmc <- NULL

mcmc$gamma <- array(NA, dim=c(Nit, B))
mcmc$eta_gamma <- array(NA, dim=c(Nit))
mcmc$zeta_gamma <- array(NA, dim=c(Nit, B))
mcmc$xi_gamma <- array(NA, dim=c(Nit))
mcmc$nu_gamma <- array(NA, dim=c(Nit))
mcmc$rho_gamma <- array(NA, dim=c(Nit))
mcmc$m_gamma <- array(NA, dim=c(Nit, B))

mcmc$alpha <- array(NA, dim=c(Nit, Q, S_sum, B))
mcmc$eta_alpha <- array(NA, dim=c(Nit, Q, S_sum))
mcmc$zeta_alpha <- array(NA, dim=c(Nit, Q, S_sum, B))
mcmc$xi_alpha <- array(NA, dim=c(Nit, Q, S_sum))
mcmc$nu_alpha <- array(NA, dim=c(Nit, Q, S_sum))
mcmc$rho_alpha <- array(NA, dim=c(Nit))
mcmc$m_alpha <- array(NA, dim=c(Nit, Q, S_sum, B))

mcmc$beta <- array(NA, dim=c(Nit, Q, P))
mcmc$xi_beta <- array(NA, dim=c(Nit, Q, P))
mcmc$nu_beta <- array(NA, dim=c(Nit, Q, P))
mcmc$rho_beta <- array(NA, dim=c(Nit))

mcmc$U <- array(NA, dim=c(Nit, I))
mcmc$lambda <- array(NA, dim=c(Nit, Q))
mcmc$mu <- array(NA, dim=c(Nit, Q))
mcmc$tau <- array(NA, dim=c(Nit, I, Q))
mcmc$sigma2 <- array(NA, dim=c(Nit, Q))

### Initialize
initial <- init()

mcmc$gamma[1,] <- initial$gamma
mcmc$eta_gamma[1] <- initial$eta_gamma
mcmc$zeta_gamma[1,] <- initial$zeta_gamma
mcmc$xi_gamma[1] <- initial$xi_gamma
mcmc$nu_gamma[1] <- initial$nu_gamma
mcmc$rho_gamma[1] <- initial$rho_gamma
mcmc$m_gamma[1,] <- initial$m_gamma

mcmc$alpha[1,,,] <- initial$alpha
mcmc$eta_alpha[1,,] <- initial$eta_alpha
mcmc$zeta_alpha[1,,,] <- initial$zeta_alpha
mcmc$xi_alpha[1,,] <- initial$xi_alpha
mcmc$nu_alpha[1,,] <- initial$nu_alpha
mcmc$rho_alpha[1] <- initial$rho_alpha
mcmc$m_alpha[1,,,] <- initial$m_alpha

mcmc$beta[1,,] <- initial$beta
mcmc$xi_beta[1,,] <- initial$xi_beta
mcmc$nu_beta[1,,] <- initial$nu_beta
mcmc$rho_beta[1] <- initial$rho_beta

mcmc$U[1,] <- initial$U
mcmc$lambda[1,] <- initial$lambda
mcmc$mu[1,] <- initial$mu
mcmc$tau[1,,] <- initial$tau
mcmc$sigma2[1,] <- initial$sigma2


### Start of the chain
start.time = proc.time()

for (nit in 2:Nit){
  
  print(paste(seed_index,nit,sep=","))
  
  # update gamma
  mcmc$eta_gamma[nit] <- update_eta_gamma_rcpp(mcmc$zeta_gamma[nit-1,], mcmc$xi_gamma[nit-1],
                                               mcmc$nu_gamma[nit-1], mcmc$alpha[nit-1,,,], mcmc$beta[nit-1,,], mcmc$U[nit-1,], 
                                               mcmc$lambda[nit-1,], mcmc$mu[nit-1,], mcmc$tau[nit-1,,], mcmc$sigma2[nit-1,],
                                               I, S_sum, B, P, Y, A_bs, X_bs, Z)
  mcmc$zeta_gamma[nit,] <- update_zeta_gamma_rcpp(mcmc$eta_gamma[nit], mcmc$m_gamma[nit-1,], 
                                                  mcmc$alpha[nit-1,,,], mcmc$beta[nit-1,,], mcmc$U[nit-1,], 
                                                  mcmc$lambda[nit-1,], mcmc$mu[nit-1,], mcmc$tau[nit-1,,], mcmc$sigma2[nit-1,],
                                                  I, S_sum, B, P, Y, A_bs, X_bs, Z)
  gamma_update_list <- update_gamma(mcmc$eta_gamma[nit], mcmc$zeta_gamma[nit,])
  mcmc$eta_gamma[nit] <- gamma_update_list$eta_gamma
  mcmc$zeta_gamma[nit,] <- gamma_update_list$zeta_gamma
  mcmc$gamma[nit,] <- gamma_update_list$gamma
  mcmc$xi_gamma[nit] <- update_xi_gamma(mcmc$eta_gamma[nit], mcmc$nu_gamma[nit-1], mcmc$rho_gamma[nit-1])
  mcmc$nu_gamma[nit] <- update_nu_gamma(mcmc$eta_gamma[nit], mcmc$xi_gamma[nit])
  mcmc$rho_gamma[nit] <- update_rho_gamma(mcmc$xi_gamma[nit])
  mcmc$m_gamma[nit,] <- update_m_gamma(mcmc$zeta_gamma[nit,])
  
  # update alpha
  mcmc$eta_alpha[nit,,] <- update_eta_alpha_rcpp(mcmc$gamma[nit,], mcmc$eta_alpha[nit-1,,], mcmc$zeta_alpha[nit-1,,,], 
                                                 mcmc$xi_alpha[nit-1,,], mcmc$nu_alpha[nit-1,,], mcmc$beta[nit-1,,], mcmc$U[nit-1,],
                                                 mcmc$lambda[nit-1,], mcmc$mu[nit-1,], mcmc$tau[nit-1,,], mcmc$sigma2[nit-1,],
                                                 I, Q, S, S_sum, B, P, Y, A, A_bs, X_bs, Z)
  mcmc$zeta_alpha[nit,,,] <- update_zeta_alpha_rcpp(mcmc$gamma[nit,], mcmc$eta_alpha[nit,,], mcmc$zeta_alpha[nit-1,,,], 
                                                    mcmc$m_alpha[nit-1,,,], mcmc$beta[nit-1,,], mcmc$U[nit-1,], mcmc$lambda[nit-1,],
                                                    mcmc$mu[nit-1,], mcmc$tau[nit-1,,], mcmc$sigma2[nit-1,],
                                                    I, Q, S, S_sum, B, P, Y, A, A_bs, X_bs, Z)
  alpha_update_list <- update_alpha(mcmc$eta_alpha[nit,,], mcmc$zeta_alpha[nit,,,])
  mcmc$eta_alpha[nit,,] <- alpha_update_list$eta_alpha
  mcmc$zeta_alpha[nit,,,] <- alpha_update_list$zeta_alpha
  mcmc$alpha[nit,,,] <- alpha_update_list$alpha
  
  mcmc$alpha[nit,,(S*2+1):S_sum,] <- 0
  mcmc$xi_alpha[nit,,] <- update_xi_alpha(mcmc$eta_alpha[nit,,], mcmc$nu_alpha[nit-1,,], mcmc$rho_alpha[nit-1])
  mcmc$nu_alpha[nit,,] <- update_nu_alpha(mcmc$eta_alpha[nit,,], mcmc$xi_alpha[nit,,])
  mcmc$rho_alpha[nit] <- update_rho_alpha(mcmc$xi_alpha[nit,,])
  mcmc$m_alpha[nit,,,] <- update_m_alpha(mcmc$zeta_alpha[nit,,,])
  
  # update beta
  # mcmc$beta[nit,,] <- update_beta_rcpp(mcmc$xi_beta[nit-1,,], mcmc$nu_beta[nit-1,,], mcmc$gamma[nit,], mcmc$alpha[nit,,,], 
  #                                      mcmc$U[nit-1,], mcmc$lambda[nit-1,], mcmc$mu[nit-1,], mcmc$tau[nit-1,,], mcmc$sigma2[nit-1,],
  #                                      I, Q, S, S_sum, B, P, Y, A, A_bs, X_bs, Z)
  # mcmc$xi_beta[nit,,] <- update_xi_beta(mcmc$beta[nit,,], mcmc$nu_beta[nit-1,,], mcmc$rho_beta[nit-1])
  # mcmc$nu_beta[nit,,] <- update_nu_beta(mcmc$beta[nit,,], mcmc$xi_beta[nit,,])
  # mcmc$rho_beta[nit] <- update_rho_beta(mcmc$xi_beta[nit,,])
  mcmc$beta[nit,,] <- 0
  
  # update U and lambda
  if (noUC){
    mcmc$U[nit,] <- mcmc$lambda[nit,] <- 0
  }else{
    mcmc$U[nit,] <- update_U_rcpp(mcmc$gamma[nit,], mcmc$alpha[nit,,,], mcmc$beta[nit,,], mcmc$lambda[nit-1,], mcmc$mu[nit-1,],
                                  mcmc$tau[nit-1,,], mcmc$sigma2[nit-1,], I, Q, S, S_sum, B, P, Y, A, A_bs, X_bs, Z)
    mcmc$lambda[nit,] <- update_lambda_rcpp(mcmc$gamma[nit,], mcmc$alpha[nit,,,], mcmc$beta[nit,,], mcmc$U[nit,],
                                            mcmc$mu[nit-1,], mcmc$tau[nit-1,,], mcmc$sigma2[nit-1,],
                                            I, Q, S, S_sum, B, P, Y, A, A_bs, X_bs, Z)
    mcmc$lambda[nit,] <- abs(mcmc$lambda[nit,])
  }
  
  # update mu, tau, and sigma2
  mcmc$mu[nit,] <- update_mu(mcmc$gamma[nit,], mcmc$alpha[nit,,,], mcmc$beta[nit,,], mcmc$U[nit,], 
                             mcmc$lambda[nit,], mcmc$tau[nit-1,,], mcmc$sigma2[nit-1,])
  mcmc$tau[nit,,] <- update_tau_rcpp(mcmc$gamma[nit,], mcmc$alpha[nit,,,], mcmc$beta[nit,,], mcmc$U[nit,], mcmc$lambda[nit,], 
                                     mcmc$mu[nit,], mcmc$sigma2[nit-1,], I, Q, S, S_sum, B, P, Y, A, A_bs, X_bs, Z)
  mcmc$sigma2[nit,] <- update_sigma2(mcmc$gamma[nit,], mcmc$alpha[nit,,,], mcmc$beta[nit,,], mcmc$U[nit,], 
                                     mcmc$lambda[nit,], mcmc$mu[nit,], mcmc$tau[nit,,])
}

duration = proc.time()-start.time
### End of the chain


### posterior samples
thin.fac <- 10 # thinning factor 
post_index <- seq(burn.in+1, Nit, by=thin.fac) # index of posterior samples

post <- NULL

post$gamma <- mcmc$gamma[post_index,]
post$eta_gamma <- mcmc$eta_gamma[post_index]
post$zeta_gamma <- mcmc$zeta_gamma[post_index,]
post$xi_gamma <- mcmc$xi_gamma[post_index]
post$nu_gamma <- mcmc$nu_gamma[post_index]
post$rho_gamma <- mcmc$rho_gamma[post_index]
post$m_gamma <- mcmc$m_gamma[post_index,]

post$alpha <- mcmc$alpha[post_index,,,]
post$eta_alpha <- mcmc$eta_alpha[post_index,,]
post$zeta_alpha <- mcmc$zeta_alpha[post_index,,,]
post$xi_alpha <- mcmc$xi_alpha[post_index,,]
post$nu_alpha <- mcmc$nu_alpha[post_index,,]
post$rho_alpha <- mcmc$rho_alpha[post_index]
post$m_alpha <- mcmc$m_alpha[post_index,,,]

post$beta <- mcmc$beta[post_index,,]
post$xi_beta <- mcmc$xi_beta[post_index,,]
post$nu_beta <- mcmc$nu_beta[post_index,,]
post$rho_beta <- mcmc$rho_beta[post_index]

post$U <- mcmc$U[post_index,]
post$lambda <- mcmc$lambda[post_index,]
post$mu <- mcmc$mu[post_index,]
post$tau <- mcmc$tau[post_index,,]
post$sigma2 <- mcmc$sigma2[post_index,]