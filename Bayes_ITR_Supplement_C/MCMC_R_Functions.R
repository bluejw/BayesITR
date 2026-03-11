################################################## Simulation Functions ###########################################


init <- function(){
  
  # MCMC initialization 
  # initialize parameters need to be estimated in MCMC
  
  gamma_init <- array(NA, dim=c(B))
  eta_gamma_init <- NA
  zeta_gamma_init <- array(NA, dim=c(B))
  xi_gamma_init <- NA
  nu_gamma_init <- NA
  rho_gamma_init <- 0.5
  m_gamma_init <- array(NA, dim=c(B))
  
  alpha_init <- array(NA, dim=c(Q, S_sum, B))
  eta_alpha_init <- array(NA, dim=c(Q, S_sum))
  zeta_alpha_init <- array(NA, dim=c(Q, S_sum, B))
  xi_alpha_init <- array(NA, dim=c(Q, S_sum))
  nu_alpha_init <- array(NA, dim=c(Q, S_sum))
  rho_alpha_init <- 0.5
  m_alpha_init <- array(NA, dim=c(Q, S_sum, B))
  
  beta_init <- array(NA, dim=c(Q, P))
  xi_beta_init <- array(NA, dim=c(Q, P))
  nu_beta_init <- array(NA, dim=c(Q, P))
  rho_beta_init <- 0.5
  
  U_init <- rnorm(I, 0, 1)
  lambda_init <- rep(0, Q)
  mu_init <- rep(0, Q)
  tau_init <- matrix(1, nrow=I, ncol=Q) 
  sigma2_init <- rep(1, Q)
  
  # initialize gamma
  xi_gamma_init <- sample(c(nu_0,1), 1, prob=c(0.5,0.5)); nu_gamma_init <- 0.01
  eta_gamma_init <- rnorm(1, 0, sqrt(xi_gamma_init*nu_gamma_init))
  for (b in 1:B){
    m_gamma_init[b] <- sample(c(1,-1), 1, prob=c(0.5, 0.5))
    zeta_gamma_init[b] <- rnorm(1, m_gamma_init[b], 1) 
  }
  gamma_init <- eta_gamma_init*zeta_gamma_init
  
  # initialize alpha
  for (q in 1:Q){
    for (s in 1:S_sum){
      xi_alpha_init[q,s] <- sample(c(nu_0,1), 1, prob=c(0.5,0.5)); nu_alpha_init[q,s] <- 0.01
      eta_alpha_init[q,s] <- rnorm(1, 0, sqrt(xi_alpha_init[q,s]*nu_alpha_init[q,s]))
      for (b in 1:B){
        m_alpha_init[q,s,b] <- sample(c(1,-1), 1, prob=c(0.5, 0.5))
        zeta_alpha_init[q,s,b] <- rnorm(1, m_alpha_init[q,s,b], 1) 
      }
      alpha_init[q,s,] <- eta_alpha_init[q,s]*zeta_alpha_init[q,s,]
    }
  }
  
  # initialize beta
  for (q in 1:Q){
    for (p in 1:P){
      xi_beta_init[q,p] <- sample(c(nu_0,1), 1, prob=c(0.5,0.5)); nu_beta_init[q,p] <- 0.01
      beta_init[q,p] <- rnorm(1, 0, sqrt(xi_beta_init[q,p]*nu_beta_init[q,p]))
    }
  }
  
  init_list <- list(gamma=gamma_init, eta_gamma=eta_gamma_init, zeta_gamma=zeta_gamma_init, 
                    xi_gamma=xi_gamma_init, nu_gamma=nu_gamma_init, rho_gamma=rho_gamma_init, m_gamma=m_gamma_init, 
                    alpha=alpha_init, eta_alpha=eta_alpha_init, zeta_alpha=zeta_alpha_init, 
                    xi_alpha=xi_alpha_init, nu_alpha=nu_alpha_init, rho_alpha=rho_alpha_init, m_alpha=m_alpha_init, 
                    beta=beta_init, xi_beta=xi_beta_init, nu_beta=nu_beta_init, rho_beta=rho_beta_init, 
                    U=U_init, lambda=lambda_init, mu=mu_init, tau=tau_init, sigma2=sigma2_init)
  return(init_list)            
}


update_eta_gamma <- function(zeta_gamma, xi_gamma, nu_gamma, alpha, beta, U, lambda, mu, tau, sigma2){
  
  # calculate A_sum and Ay_sum
  A_sum <- t(A_bs%*%zeta_gamma)%*%diag(tau[,1])%*%(A_bs%*%zeta_gamma)
  Y_tilde <- Y - mu[1] - lambda[1]*U - Z%*%beta[1,]
  for (s in 1:S_sum){ Y_tilde <- Y_tilde - X_bs[,s,]%*%alpha[1,s,] }
  Ay_sum <- t(A_bs%*%zeta_gamma)%*%diag(tau[,1])%*%Y_tilde
  
  # Gaussian posterior distribution
  V_n <- 1/(A_sum/sigma2[1] + 1/(xi_gamma*nu_gamma))
  mu_n <- V_n*(Ay_sum/sigma2[1])
  eta_gamma_update <- rnorm(1, mu_n, sqrt(V_n))
  
  return(eta_gamma_update)
}


update_zeta_gamma <- function(eta_gamma, m_gamma, alpha, beta, U, lambda, mu, tau, sigma2){
  
  # calculate A_sum and Ay_sum
  A_sum <- t(A_bs*eta_gamma)%*%diag(tau[,1])%*%(A_bs*eta_gamma)
  Y_tilde <- Y - mu[1] - lambda[1]*U - Z%*%beta[1,]
  for (s in 1:S_sum){ Y_tilde <- Y_tilde - X_bs[,s,]%*%alpha[1,s,] }
  Ay_sum <- t(A_bs*eta_gamma)%*%diag(tau[,1])%*%Y_tilde
  
  # Gaussian posterior distribution
  V_n <- chol2inv(chol(A_sum/sigma2[1] + diag(1,B)))
  mu_n <- V_n%*%(Ay_sum/sigma2[1] + m_gamma)
  zeta_gamma_update <- as.vector(rmvn_rcpp(1, mu_n, V_n))
  
  return(zeta_gamma_update)
}


update_xi_gamma <- function(eta_gamma, nu_gamma, rho_gamma){
  
  # prob = pr(xi_gamma=1) / pr(xi_gamma=nu_0)
  prob <- (sqrt(nu_0)*rho_gamma)/(1-rho_gamma)*exp((1-nu_0)*eta_gamma^2/(2*nu_0*nu_gamma))
  if (prob == Inf){ xi_gamma_update <- 1
  }else{ xi_gamma_update <- sample(c(1,nu_0), 1, prob=c(prob, 1)) }
  
  return(xi_gamma_update)
}


update_nu_gamma <- function(eta_gamma, xi_gamma){
  
  a_nu_star <- a_nu + 1/2
  b_nu_star <- b_nu + eta_gamma^2/(2*xi_gamma)
  nu_gamma_update <- rinvgamma(1, a_nu_star, b_nu_star)
  
  return(nu_gamma_update)
}


update_rho_gamma <- function(xi_gamma){
  
  a_rho_star <- a_rho + sum(xi_gamma == 1)
  b_rho_star <- b_rho + sum(xi_gamma == nu_0)
  rho_gamma_update <- rbeta(1, a_rho_star, b_rho_star)
  
  return(rho_gamma_update)
}


update_m_gamma <- function(zeta_gamma){
  
  m_gamma_update <- rep(NA, B)
  for (b in 1:B){
    # prob = pr(m_gamma_b = 1)
    prob <- 1/(1+exp(-2*zeta_gamma[b]))
    m_gamma_update[b] <- sample(c(1,-1), 1, prob=c(prob, 1-prob)) 
  }
  
  return(m_gamma_update)
}


update_gamma <- function(eta_gamma, zeta_gamma){
  
  eta_gamma_update <- eta_gamma
  zeta_gamma_update <- zeta_gamma
  gamma_update <- rep(NA, B)
  
  # normalization for identifiability
  zeta_gamma_bar <- mean(abs(zeta_gamma))
  zeta_gamma_update <- zeta_gamma_update/zeta_gamma_bar
  eta_gamma_update <- eta_gamma_update*zeta_gamma_bar
  gamma_update <- eta_gamma_update*zeta_gamma_update
  
  returnlist <- list(eta_gamma=eta_gamma_update, 
                     zeta_gamma=zeta_gamma_update, 
                     gamma=gamma_update)
  return(returnlist)
}


update_eta_alpha <- function(gamma, eta_alpha, zeta_alpha, xi_alpha, nu_alpha, beta, U, lambda, mu, tau, sigma2){
  
  eta_alpha_update <- eta_alpha
  
  for (s in 1:S_sum){
    # calculate X_sum and Xy_sum
    X_sum <- t(X_bs[,s,]%*%zeta_alpha[1,s,])%*%diag(tau[,1])%*%(X_bs[,s,]%*%zeta_alpha[1,s,])
    Y_tilde <- Y - mu[1] - lambda[1]*U - A_bs%*%gamma - Z%*%beta[1,] + X_bs[,s,]%*%(eta_alpha_update[1,s]*zeta_alpha[1,s,])
    for (s2 in 1:S_sum){ Y_tilde <- Y_tilde - X_bs[,s2,]%*%(eta_alpha_update[1,s2]*zeta_alpha[1,s2,]) }
    Xy_sum <- t(X_bs[,s,]%*%zeta_alpha[1,s,])%*%diag(tau[,1])%*%Y_tilde
    
    # Gaussian posterior distribution
    V_n <- 1/(X_sum/sigma2[1] + 1/(xi_alpha[1,s]*nu_alpha[1,s]))
    mu_n <- V_n*(Xy_sum/sigma2[1])
    eta_alpha_update[1,s] <- rnorm(1, mu_n, sqrt(V_n))
  }
  
  for (s in 1:S){
    # calculate X_sum and Xa_sum
    X_sum <- t(X_bs[,s,]%*%zeta_alpha[2,s,])%*%diag(tau[,2])%*%(X_bs[,s,]%*%zeta_alpha[2,s,])
    A_tilde <- A - mu[2] - lambda[2]*U - Z%*%beta[2,] + X_bs[,s,]%*%(eta_alpha_update[2,s]*zeta_alpha[2,s,])
    for (s2 in 1:S){ A_tilde <- A_tilde - X_bs[,s2,]%*%(eta_alpha_update[2,s2]*zeta_alpha[2,s2,]) }
    Xa_sum <- t(X_bs[,s,]%*%zeta_alpha[2,s,])%*%diag(tau[,2])%*%A_tilde
    
    # Gaussian posterior distribution
    V_n <- 1/(X_sum/sigma2[2] + 1/(xi_alpha[2,s]*nu_alpha[2,s]))
    mu_n <- V_n*(Xa_sum/sigma2[2])
    eta_alpha_update[2,s] <- rnorm(1, mu_n, sqrt(V_n))
  }
  eta_alpha_update[2,(S+1):S_sum] <- 0
  
  return(eta_alpha_update)
}


update_zeta_alpha <- function(gamma, eta_alpha, zeta_alpha, m_alpha, beta, U, lambda, mu, tau, sigma2){
  
  zeta_alpha_update <- zeta_alpha
  
  for (s in 1:S_sum){
    # calculate X_sum and Xy_sum
    X_sum <- t(X_bs[,s,]*eta_alpha[1,s])%*%diag(tau[,1])%*%(X_bs[,s,]*eta_alpha[1,s])
    Y_tilde <- Y - mu[1] - lambda[1]*U - A_bs%*%gamma - Z%*%beta[1,] + X_bs[,s,]%*%(eta_alpha[1,s]*zeta_alpha_update[1,s,])
    for (s2 in 1:S_sum){ Y_tilde <- Y_tilde - X_bs[,s2,]%*%(eta_alpha[1,s2]*zeta_alpha_update[1,s2,]) }
    Xy_sum <- t(X_bs[,s,]*eta_alpha[1,s])%*%diag(tau[,1])%*%Y_tilde
    
    # Gaussian posterior distribution
    V_n <- chol2inv(chol(X_sum/sigma2[1] + diag(1,B)))
    mu_n <- V_n%*%(Xy_sum/sigma2[1] + m_alpha[1,s,])
    zeta_alpha_update[1,s,] <- as.vector(rmvn_rcpp(1, mu_n, V_n))
  }
  
  for (s in 1:S){
    # calculate X_sum and Xa_sum
    X_sum <- t(X_bs[,s,]*eta_alpha[2,s])%*%diag(tau[,2])%*%(X_bs[,s,]*eta_alpha[2,s])
    A_tilde <- A - mu[2] - lambda[2]*U - Z%*%beta[2,] + X_bs[,s,]%*%(eta_alpha[2,s]*zeta_alpha_update[2,s,])
    for (s2 in 1:S){ A_tilde <- A_tilde - X_bs[,s2,]%*%(eta_alpha[2,s2]*zeta_alpha_update[2,s2,]) }
    Xa_sum <- t(X_bs[,s,]*eta_alpha[2,s])%*%diag(tau[,2])%*%A_tilde
    
    # Gaussian posterior distribution
    V_n <- chol2inv(chol(X_sum/sigma2[2] + diag(1,B)))
    mu_n <- V_n%*%(Xa_sum/sigma2[2] + m_alpha[2,s,])
    zeta_alpha_update[2,s,] <- as.vector(rmvn_rcpp(1, mu_n, V_n))
  }
  zeta_alpha_update[2,(S+1):S_sum,] <- 0
  
  return(zeta_alpha_update)
}


update_xi_alpha <- function(eta_alpha, nu_alpha, rho_alpha){
  
  xi_alpha_update <- matrix(NA, nrow=Q, ncol=S_sum) 
  
  for (q in 1:Q){
    for (s in 1:S_sum){
      # prob = pr(xi_alpha_qs=1) / pr(xi_alpha_qs=nu_0)
      prob <- (sqrt(nu_0)*rho_alpha)/(1-rho_alpha)*exp((1-nu_0)*eta_alpha[q,s]^2/(2*nu_0*nu_alpha[q,s]))
      if (prob == Inf){ xi_alpha_update[q,s] <- 1
      }else{ xi_alpha_update[q,s] <- sample(c(1,nu_0), 1, prob=c(prob, 1)) }
    }
  }
  
  return(xi_alpha_update)
}


update_nu_alpha <- function(eta_alpha, xi_alpha){
  
  nu_alpha_update <- matrix(NA, nrow=Q, ncol=S_sum) 
  
  for (q in 1:Q){
    for (s in 1:S_sum){
      a_nu_star <- a_nu + 1/2
      b_nu_star <- b_nu + eta_alpha[q,s]^2/(2*xi_alpha[q,s])
      nu_alpha_update[q,s] <- rinvgamma(1, a_nu_star, b_nu_star)
    }
  }
  
  return(nu_alpha_update)
}


update_rho_alpha <- function(xi_alpha){
  
  a_rho_star <- a_rho + sum(xi_alpha == 1)
  b_rho_star <- b_rho + sum(xi_alpha == nu_0)
  rho_alpha_update <- rbeta(1, a_rho_star, b_rho_star)
  
  return(rho_alpha_update)
}


update_m_alpha <- function(zeta_alpha){
  
  m_alpha_update <-  array(NA, dim=c(Q, S_sum, B))
  
  for (q in 1:Q){
    for (s in 1:S_sum){
      for (b in 1:B){
        # prob = pr(m_alpha_qsb = 1)
        prob <- 1/(1+exp(-2*zeta_alpha[q,s,b]))
        m_alpha_update[q,s,b] <- sample(c(1,-1), 1, prob=c(prob, 1-prob)) 
      }
    }
  }
  
  return(m_alpha_update)
}


update_alpha <- function(eta_alpha, zeta_alpha){
  
  eta_alpha_update <- eta_alpha
  zeta_alpha_update <- zeta_alpha
  alpha_update <- array(NA, dim=c(Q, S_sum, B))
  
  for (s in 1:S_sum){
    # normalization for identifiability
    zeta_alpha_1s_bar <- mean(abs(zeta_alpha[1,s,]))
    zeta_alpha_update[1,s,] <- zeta_alpha_update[1,s,]/zeta_alpha_1s_bar
    eta_alpha_update[1,s] <- eta_alpha_update[1,s]*zeta_alpha_1s_bar
    alpha_update[1,s,] <- eta_alpha_update[1,s]*zeta_alpha_update[1,s,]
  }
  
  for (s in 1:S){
    # normalization for identifiability
    zeta_alpha_2s_bar <- mean(abs(zeta_alpha[2,s,]))
    zeta_alpha_update[2,s,] <- zeta_alpha_update[2,s,]/zeta_alpha_2s_bar
    eta_alpha_update[2,s] <- eta_alpha_update[2,s]*zeta_alpha_2s_bar
    alpha_update[2,s,] <- eta_alpha_update[2,s]*zeta_alpha_update[2,s,]
  }
  alpha_update[2,(S+1):S_sum,] <- 0
  
  returnlist <- list(eta_alpha=eta_alpha_update, 
                     zeta_alpha=zeta_alpha_update, 
                     alpha=alpha_update)
  return(returnlist)
}


update_beta <- function(xi_beta, nu_beta, gamma, alpha, U, lambda, mu, tau, sigma2){
  
  beta_update <- matrix(NA, nrow=Q, ncol=P) 
  
  # calculate mu_m and V_m
  Sigma_inv <- diag(tau[,1]/sigma2[1]) 
  Sigma_beta_inv <- diag(1/(xi_beta[1,]*nu_beta[1,]))
  V_m <- chol2inv(chol(t(Z)%*%Sigma_inv%*%Z + Sigma_beta_inv))
  Y_tilde <- Y - mu[1] - lambda[1]*U - A_bs%*%gamma
  for (s in 1:S_sum){ Y_tilde <- Y_tilde - X_bs[,s,]%*%alpha[1,s,] }
  mu_m <- V_m%*%(t(Z)%*%Sigma_inv%*%Y_tilde)
  
  # Gaussian posterior distribution
  beta_update[1,] <- rmvn_rcpp(1, mu_m, V_m)
  
  # calculate mu_n and V_n
  Sigma_inv <- diag(tau[,2]/sigma2[2]) 
  Sigma_beta_inv <- diag(1/(xi_beta[2,]*nu_beta[2,]))
  V_n <- chol2inv(chol(t(Z)%*%Sigma_inv%*%Z + Sigma_beta_inv))
  A_tilde <- A - mu[2] - lambda[2]*U 
  for (s in 1:S){ A_tilde <- A_tilde - X_bs[,s,]%*%alpha[2,s,] }
  mu_n <- V_n%*%(t(Z)%*%Sigma_inv%*%A_tilde)
  
  # Gaussian posterior distribution
  beta_update[2,] <- rmvn_rcpp(1, mu_n, V_n)
  
  return(beta_update)
}


update_xi_beta <- function(beta, nu_beta, rho_beta){
  
  xi_beta_update <- matrix(NA, nrow=Q, ncol=P) 
  
  for (q in 1:Q){
    for (p in 1:P){
      # prob = pr(xi_beta_qp=1) / pr(xi_beta_qp=nu_0)
      prob <- (sqrt(nu_0)*rho_beta)/(1-rho_beta)*exp((1-nu_0)*beta[q,p]^2/(2*nu_0*nu_beta[q,p]))
      if (prob == Inf){ xi_beta_update[q,p] <- 1
      }else{ xi_beta_update[q,p] <- sample(c(1,nu_0), 1, prob=c(prob, 1)) }
    }
  }
  
  return(xi_beta_update)
}


update_nu_beta <- function(beta, xi_beta){
  
  nu_beta_update <- matrix(NA, nrow=Q, ncol=P) 
  
  for (q in 1:Q){
    for (p in 1:P){
      a_nu_star <- a_nu + 1/2
      b_nu_star <- b_nu + beta[q,p]^2/(2*xi_beta[q,p])
      nu_beta_update[q,p] <- rinvgamma(1, a_nu_star, b_nu_star)
    }
  }
  
  return(nu_beta_update)
}


update_rho_beta <- function(xi_beta){
  
  a_rho_star <- a_rho + sum(xi_beta == 1)
  b_rho_star <- b_rho + sum(xi_beta == nu_0)
  rho_beta_update <- rbeta(1, a_rho_star, b_rho_star)
  
  return(rho_beta_update)
}


update_U <- function(gamma, alpha, beta, lambda, mu, tau, sigma2){
  
  U_update <- rep(NA, I)
  
  # caculate YA_tilde 
  YA_tilde <- matrix(NA, nrow=I, ncol=Q)
  YA_tilde[,1] <- Y - mu[1] - A_bs%*%gamma - Z%*%beta[1,]
  for (s in 1:S_sum){ YA_tilde[,1] <- YA_tilde[,1] - X_bs[,s,]%*%alpha[1,s,] }
  YA_tilde[,2] <- A - mu[2] - Z%*%beta[2,]
  for (s in 1:S){ YA_tilde[,2] <- YA_tilde[,2] - X_bs[,s,]%*%alpha[2,s,] }
  
  # update U_i
  for (i in 1:I){
    Sigma_i_inv <- diag(tau[i,]/sigma2) 
    # calculate mu_u and V_u
    V_u <- 1/(t(lambda)%*%Sigma_i_inv%*%lambda + 1)
    mu_u <- V_u*(t(lambda)%*%Sigma_i_inv%*%YA_tilde[i,])
    # Gaussian posterior distribution
    U_update[i] <- rnorm(1, mu_u, sqrt(V_u))
  }
  
  return(U_update)
}


update_lambda <- function(gamma, alpha, beta, U, mu, tau, sigma2){
  
  lambda_update <- rep(NA, Q)
  
  # calculate mu_m and V_m
  Sigma_inv <- diag(tau[,1]/sigma2[1]) 
  V_m <- 1/(1/100 + t(U)%*%Sigma_inv%*%U)
  Y_tilde <- Y - mu[1] - A_bs%*%gamma - Z%*%beta[1,]
  for (s in 1:S_sum){ Y_tilde <- Y_tilde - X_bs[,s,]%*%alpha[1,s,] }
  mu_m <- V_m%*%(t(U)%*%Sigma_inv%*%Y_tilde)
  
  # Gaussian posterior distribution
  lambda_update[1] <- rnorm(1, mu_m, sqrt(V_m))
  
  # calculate mu_n and V_n
  Sigma_inv <- diag(tau[,2]/sigma2[2]) 
  V_n <- 1/(1/100 + t(U)%*%Sigma_inv%*%U)
  A_tilde <- A - mu[2] - Z%*%beta[2,]
  for (s in 1:S){ A_tilde <- A_tilde - X_bs[,s,]%*%alpha[2,s,] }
  mu_n <- V_n%*%(t(U)%*%Sigma_inv%*%A_tilde)
  
  # Gaussian posterior distribution
  lambda_update[2] <- rnorm(1, mu_n, sqrt(V_n))
  
  # Gaussian posterior distribution
  # mu_mn <- c(mu_m, mu_n); V_mn <- diag(c(V_m, V_n))
  # if (trunc_norm){ lambda_update <- rtmvnorm(n=1, mu=mu_mn, sigma=V_mn, lb=rep(0,Q), ub=rep(Inf,Q))
  # }else{ lambda_update <- rmvn_rcpp(1, mu_mn, V_mn) }
  
  return(lambda_update)
}


update_mu <- function(gamma, alpha, beta, U, lambda, tau, sigma2){
  
  mu_update <- rep(NA, Q)
  
  Y_tilde <- Y - lambda[1]*U - A_bs%*%gamma - Z%*%beta[1,]
  for (s in 1:S_sum){ Y_tilde <- Y_tilde - X_bs[,s,]%*%alpha[1,s,] }
  Y_tilde_sum <- sum(Y_tilde*tau[,1])
  
  # Gaussian posterior distribution
  V_n <- 1/(1/100 + sum(tau[,1])/sigma2[1])
  mu_n <- V_n*Y_tilde_sum/sigma2[1]
  mu_update[1] <- rnorm(1, mu_n, sqrt(V_n))
  
  A_tilde <- A - lambda[2]*U - Z%*%beta[2,]
  for (s in 1:S){ A_tilde <- A_tilde - X_bs[,s,]%*%alpha[2,s,] }
  A_tilde_sum <- sum(A_tilde*tau[,2])
  
  # Gaussian posterior distribution
  V_n <- 1/(1/100 + sum(tau[,2])/sigma2[2]);
  mu_n <- V_n*A_tilde_sum/sigma2[2];
  mu_update[2] <- rnorm(1, mu_n, sqrt(V_n))
  
  return(mu_update)
}


update_tau <- function(gamma, alpha, beta, U, lambda, mu, sigma2){
  
  tau_update <- matrix(NA, nrow=I, ncol=Q) 
  eps <- 1e-3 # avoid numerical issues
  
  Y_tilde <- Y - mu[1] - lambda[1]*U - A_bs%*%gamma - Z%*%beta[1,]
  for (s in 1:S_sum){ Y_tilde <- Y_tilde - X_bs[,s,]%*%alpha[1,s,] }
  for (i in 1:I){
    # Inverse-Gaussian posterior distribution
    tau_update[i,1] <- max(rinvgaussian(1, mu=sqrt(sigma2[1])/(2*abs(Y_tilde[i])), lambda=1/4), eps) 
  }
  
  A_tilde <- A - mu[2] - lambda[2]*U - Z%*%beta[2,]
  for (s in 1:S){ A_tilde <- A_tilde - X_bs[,s,]%*%alpha[2,s,] }
  for (i in 1:I){
    # Inverse-Gaussian posterior distribution
    tau_update[i,2] <- max(rinvgaussian(1, mu=sqrt(sigma2[2])/(2*abs(A_tilde[i])), lambda=1/4), eps) 
  }
  
  return(tau_update)
}


update_sigma2 <- function(gamma, alpha, beta, U, lambda, mu, tau){
  
  sigma2_update <- rep(NA, Q)
  
  Y_tilde <- Y - mu[1] - lambda[1]*U - A_bs%*%gamma - Z%*%beta[1,]
  for (s in 1:S_sum){ Y_tilde <- Y_tilde - X_bs[,s,]%*%alpha[1,s,] }
  Y_tilde_sum <- sum(Y_tilde^2*tau[,1])
  
  # Inverse-Gamma posterior distribution
  a_star <- a_sigma + I/2
  b_star <- b_sigma + Y_tilde_sum/2
  sigma2_update[1] <- rinvgamma(1, a_star, b_star)
  
  A_tilde <- A - mu[2] - lambda[2]*U - Z%*%beta[2,]
  for (s in 1:S){ A_tilde <- A_tilde - X_bs[,s,]%*%alpha[2,s,] }
  A_tilde_sum <- sum(A_tilde^2*tau[,2])
  
  # Inverse-Gamma posterior distribution
  a_star <- a_sigma + I/2
  b_star <- b_sigma + A_tilde_sum/2
  sigma2_update[2] <- rinvgamma(1, a_star, b_star)
  
  return(sigma2_update)
}


log_likelihood <- function(gamma, alpha, beta, lambda, mu, sigma2){
  
  niter <- 100 # number of samples for latent variables
  ll <- matrix(1, nrow=niter, ncol=I)
  
  for (iter in 1:niter){
    
    U <- rnorm(I, 0, 1)
    
    Y_mean <- mu[1] + lambda[1]*U + A_bs%*%gamma + Z%*%beta[1,]
    for (s in 1:S_sum){ Y_mean <- Y_mean + X_bs[,s,]%*%alpha[1,s,] }
    ll[iter,] <- ll[iter,] * dlaplace(Y, location=Y_mean, scale=2*sqrt(sigma2[1]))
    
    A_mean <- mu[2] + lambda[2]*U + Z%*%beta[2,]
    for (s in 1:S){ A_mean <- A_mean + X_bs[,s,]%*%alpha[2,s,] }
    ll[iter,] <- ll[iter,] * dlaplace(A, location=A_mean, scale=2*sqrt(sigma2[2]))
  }
  
  logll <- log(colMeans(ll))
  return(logll)
}