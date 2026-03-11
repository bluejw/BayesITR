##########################################################################################
#                            Simulation Study: Scenario II                               #
##########################################################################################

library(Rcpp)
library(RcppArmadillo)
library(Matrix)
library(MCMCpack)
library(LaplacesDemon)
library(splines2)
library(truncnorm)
library(ggplot2)
library(twangContinuous)

source("MCMC_R_Functions.R") # R functions for MCMC 
sourceCpp("MCMC_Rcpp_Functions.cpp") # Rcpp functions for MCMC 

##########################################################################################
#                          Run MCMC and SGD on Simulated Data                            #
##########################################################################################

# Simulation setups
noUC <- F # with unmeasured confounders (proposed method)
# noUC <- T # no unmeasured confounders (alternative method)

# Probability of selected data outside the region
prob <- 1 # Scenario II(a), balanced data
# prob <- 0.5 # Scenario II(b), mild unbalanced data
# prob <- 0.1 Scenario II(c), severe unbalanced data

# Average regret
num_seed <- 100 # number of seeds/experiments
regret_conservative <- rep(NA, num_seed)
regret_expected <- rep(NA, num_seed)
regret_ql <- rep(NA, num_seed)
regret_iptw <- rep(NA, num_seed)
regret_aiptw <- rep(NA, num_seed)

# Generate testing data 
set.seed(102); I_test <- 500
x_test <- rtruncnorm(I_test,a=-2,b=2,0,1)

# Start of running repeated experiments: Proposed, NonCons, NoUC, NonCons+NoUC
for (seed_index in 1:num_seed){
  
  # MCMC
  set.seed(seed_index) # set the seed to generate different simulated datasets
  source("Simu_Data_Generate.R") # generate simulation truths: Scenario II
  source("Select_Partial_Data.R") # select partial data to determine sub-scenarios
  source("MCMC_Main.R") # MCMC main function: less than 2min to run 10,000 iterations 
  
  # construct credible region
  num_samples <- 500; logll <- rep(NA, num_samples)
  for (nit in 1:num_samples){
    logll[nit] <- log_likelihood(post$gamma[nit,], post$alpha[nit,,,], post$beta[nit,,],
                                 post$lambda[nit,], post$mu[nit,], post$sigma2[nit,])
  }
  phi <- cbind(post$mu[,1],post$gamma,post$alpha[,1,1,],post$alpha[,1,3,]) # posterior samples
  index_credible_region <- order(logll,decreasing=T)[1:(num_samples*0.95)] # credible_region
  phi_cr <- phi[index_credible_region,] # posterior samples in credible_region
  
  # SGD
  source("SGD_Functions.R")
  Niter <- 1000 # number of sgd iterations
  x <- X[,1] # covariates/predictors
  std_sgd <- 0.2 # standard deviation of the policy function
  step_size <- 1e-4 # step size for sgd
  # data and transformation
  a_spline <- bSpline(A, df=B0); a_trans <- Za
  x_spline <- bSpline(X[,1], df=B0); x_trans <- Zx[1,,]
  xa_spline <- bSpline(X[,S+1], df=B0); xa_trans <- Zx[S+1,,]
  a_lower <- min(A); a_upper <- max(A)
  xa_lower <- min(X[,S+1]); xa_upper <- max(X[,S+1])
  
  # Simulation truth
  a_truth_est <- rep(NA, I_test); y_truth_est <- rep(NA, I_test)
  for (i in 1:I_test){
    opt_truth <- optim(par=c(0),fn=q_truth,lower=a_lower,upper=a_upper,
                       x=x_test[i],max_to_min=T,method='L-BFGS-B')
    y_truth_est[i] <- -opt_truth$value # from min to max
  }
  
  # Conservative reward (Proposed and NoUC)
  source("SGD_Main_Conservative.R") # SGD main function: about 10s to run 1,000 iterations 
  theta_star_conservative <- colMeans(sgd$theta[seq((Niter/2+1),Niter,by=1),])
  a_est_conservative <- policy_function(x_test, I_test, theta_star_conservative, 0) # optimal action
  y_est_conservative_truth <- q_truth(a_est_conservative, x_test, max_to_min=F) # true value of the optimal action 
  regret_conservative[seed_index] <- mean(y_truth_est-y_est_conservative_truth)
  
  # Expected reward (NonCons and NonCons+NoUC)
  source("SGD_Main_Expected.R") # SGD main function: about 5s to run 1,000 iterations 
  theta_star_expected <- colMeans(sgd$theta[seq((Niter/2+1),Niter,by=1),])
  a_est_expected <- policy_function(x_test, I_test, theta_star_expected, 0) # optimal action
  y_est_expected_truth <- q_truth(a_est_expected, x_test, max_to_min=F) # true value of the optimal action 
  regret_expected[seed_index] <- mean(y_truth_est-y_est_expected_truth)
}

# Start of running repeated experiments: Q-learning, IPTW, and AIPTW
for (seed_index in 1:num_seed){
  
  # Q-learning
  set.seed(seed_index) # set the seed to generate different simulated datasets
  source("Simu_Data_Generate.R") # generate simulation truths: Scenario II
  source("Select_Partial_Data.R") # select partial data to determine sub-scenarios
  phi_all <- lm(Y~A_bs+X_bs[,1,]+X_bs[,2,]+X_bs[,S+1,]+X_bs[,S+2,])$coef
  phi <- phi_all[c(1,2:5,6:9,14:17)]
  source('SGD_Functions.R')
  
  q_function_expected <- function(a, x, phi, sample_phi=F){
    # transfrom data to spline basis
    a_bs <- predict(a_spline, a) %*% a_trans
    x_bs <- predict(x_spline, x) %*% x_trans
    xa_bs <- predict(xa_spline, x*a) %*% xa_trans
    # calculate q-value
    y <- cbind(rep(1,length(x)),a_bs,x_bs,xa_bs) %*% phi 
    # check whether out of the boundary
    out_index <- unique(c(which(a<a_lower|a>a_upper),
                          which((x*a)<xa_lower|(x*a)>xa_upper)))
    y[out_index] <- 0
    return(y)
  }
  
  Niter <- 1000 # number of sgd iterations
  x <- X[,1] # covariates/predictors
  std_sgd <- 0.2 # standard deviation of the policy function
  step_size <- 1e-4 # step size for sgd
  # data and transformation
  a_spline <- bSpline(A, df=B0); a_trans <- Za
  x_spline <- bSpline(X[,1], df=B0); x_trans <- Zx[1,,]
  xa_spline <- bSpline(X[,S+1], df=B0); xa_trans <- Zx[S+1,,]
  a_lower <- min(A); a_upper <- max(A)
  xa_lower <- min(X[,S+1]); xa_upper <- max(X[,S+1])
  
  source("SGD_Main_Expected.R") # SGD main function: about 1s to run 1,000 iterations 
  theta_star_ql <- colMeans(sgd$theta[seq((Niter/2+1),Niter,by=1),])
  a_est_ql <- policy_function(x_test, I_test, theta_star_ql, 0) # optimal action
  y_est_ql_truth <- q_truth(a_est_ql, x_test, max_to_min=F) # true value of the optimal action 
  regret_ql[seed_index] <- mean(y_truth_est-y_est_ql_truth)
  
  # IPTW
  set.seed(seed_index) # set the seed to generate different simulated datasets
  source("Simu_Data_Generate.R") # generate simulation truths: Scenario II
  source("Select_Partial_Data.R") # select partial data to determine sub-scenarios
  df <- data.frame(A=A, X=X[,1:S])
  ps.model <- ps.cont(A~X[,1]+X[,2], data=df)
  weights <- get.weights(ps.model)
  outcome.model <- lm(Y~A_bs+X_bs[,S+1,]+X_bs[,S+2,], weights=weights)
  phi_all <- outcome.model$coef
  phi <- phi_all[c(1,2:5,6:9)]
  source('SGD_Functions.R')
  
  q_function_expected <- function(a, x, phi, sample_phi=F){
    # transfrom data to spline basis
    a_bs <- predict(a_spline, a) %*% a_trans
    xa_bs <- predict(xa_spline, x*a) %*% xa_trans
    # calculate q-value
    y <- cbind(rep(1,length(x)),a_bs,xa_bs) %*% phi 
    # check whether out of the boundary
    out_index <- unique(c(which(a<a_lower|a>a_upper),
                          which((x*a)<xa_lower|(x*a)>xa_upper)))
    y[out_index] <- 0
    return(y)
  }
  
  Niter <- 1000 # number of sgd iterations
  x <- X[,1] # covariates/predictors
  std_sgd <- 0.2 # standard deviation of the policy function
  step_size <- 1e-4 # step size for sgd
  # data and transformation
  a_spline <- bSpline(A, df=B0); a_trans <- Za
  x_spline <- bSpline(X[,1], df=B0); x_trans <- Zx[1,,]
  xa_spline <- bSpline(X[,S+1], df=B0); xa_trans <- Zx[S+1,,]
  a_lower <- min(A); a_upper <- max(A)
  xa_lower <- min(X[,S+1]); xa_upper <- max(X[,S+1])
  
  source("SGD_Main_Expected.R") # SGD main function: about 1s to run 1,000 iterations 
  theta_star_iptw <- colMeans(sgd$theta[seq((Niter/2+1),Niter,by=1),])
  a_est_iptw <- policy_function(x_test, I_test, theta_star_iptw, 0) # optimal action
  y_est_iptw_truth <- q_truth(a_est_iptw, x_test, max_to_min=F) # true value of the optimal action 
  regret_iptw[seed_index] <- mean(y_truth_est-y_est_iptw_truth)
  
  # AIPTW
  set.seed(seed_index) # set the seed to generate different simulated datasets
  source("Simu_Data_Generate.R") # generate simulation truths: Scenario II
  source("Select_Partial_Data.R") # select partial data to determine sub-scenarios
  df <- data.frame(A=A, X=X[,1:S])
  ps.model <- ps.cont(A~X[,1]+X[,2], data=df)
  weights <- get.weights(ps.model)
  outcome.model <- lm(Y~A_bs+X_bs[,1,]+X_bs[,2,]+X_bs[,S+1,]+X_bs[,S+2,], weights=weights)
  phi_all <- outcome.model$coef
  phi <- phi_all[c(1,2:5,6:9,14:17)]
  source('SGD_Functions.R')
  
  q_function_expected <- function(a, x, phi, sample_phi=F){
    # transfrom data to spline basis
    a_bs <- predict(a_spline, a) %*% a_trans
    x_bs <- predict(x_spline, x) %*% x_trans
    xa_bs <- predict(xa_spline, x*a) %*% xa_trans
    # calculate q-value
    y <- cbind(rep(1,length(x)),a_bs,x_bs,xa_bs) %*% phi 
    # check whether out of the boundary
    out_index <- unique(c(which(a<a_lower|a>a_upper),
                          which((x*a)<xa_lower|(x*a)>xa_upper)))
    y[out_index] <- 0
    return(y)
  }
  
  Niter <- 1000 # number of sgd iterations
  x <- X[,1] # covariates/predictors
  std_sgd <- 0.2 # standard deviation of the policy function
  step_size <- 1e-4 # step size for sgd
  # data and transformation
  a_spline <- bSpline(A, df=B0); a_trans <- Za
  x_spline <- bSpline(X[,1], df=B0); x_trans <- Zx[1,,]
  xa_spline <- bSpline(X[,S+1], df=B0); xa_trans <- Zx[S+1,,]
  a_lower <- min(A); a_upper <- max(A)
  xa_lower <- min(X[,S+1]); xa_upper <- max(X[,S+1])
  
  source("SGD_Main_Expected.R") # SGD main function: about 1s to run 1,000 iterations 
  theta_star_aiptw <- colMeans(sgd$theta[seq((Niter/2+1),Niter,by=1),])
  a_est_aiptw <- policy_function(x_test, I_test, theta_star_aiptw, 0) # optimal action
  y_est_aiptw_truth <- q_truth(a_est_aiptw, x_test, max_to_min=F) # true value of the optimal action 
  regret_aiptw[seed_index] <- mean(y_truth_est-y_est_aiptw_truth)
}

##########################################################################################
#                                   Figure 4                                             #
##########################################################################################

load("SGD_Results_II(a).RData") # simulation results: Scenario II(a)
load("SGD_Results_II(b).RData") # simulation results: Scenario II(b)
load("SGD_Results_II(c).RData") # simulation results: Scenario II(c)

# Figure 4
df <- data.frame(x=c(rep("Scenario II(a)",num_seed*7),
                     rep("Scenario II(b)",num_seed*7),
                     rep("Scenario II(c)",num_seed*7)),
                 y=c(res10$regret_conservative, res10$regret_expected, 
                     res10$regret_conservative_noUC, res10$regret_expected_noUC, 
                     res10$regret_ql, res10$regret_iptw, res10$regret_aiptw,
                     res5$regret_conservative, res5$regret_expected, 
                     res5$regret_conservative_noUC, res5$regret_expected_noUC, 
                     res5$regret_ql, res5$regret_iptw, res5$regret_aiptw,
                     res1$regret_conservative, res1$regret_expected, 
                     res1$regret_conservative_noUC, res1$regret_expected_noUC, 
                     res1$regret_ql, res1$regret_iptw, res1$regret_aiptw),
                 z=factor(rep(c(rep("Proposed",num_seed), rep("NonCons",num_seed),
                                rep("NoUC",num_seed), rep("NonCons+NoUC",num_seed),
                                rep("Q-Learning",num_seed), rep("IPTW",num_seed), rep("AIPTW",num_seed)),3),
                          levels=c("Proposed","NonCons","NoUC","NonCons+NoUC","Q-Learning","IPTW","AIPTW")))
p <- ggplot(data=df, aes(x=x, y=y, fill=z)) + 
  stat_boxplot(geom ='errorbar', width=0.35, position=position_dodge(0.75)) + 
  geom_boxplot() + xlab(" ") + ylab("Average Regret") + theme_light() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) + 
  theme(legend.title=element_text(size=25), legend.text=element_text(size=25)) + 
  scale_y_continuous(limits=c(0,0.3)) + 
  scale_fill_manual(name="Method", values=c("white","lightblue1","skyblue3","steelblue4","lightgrey","darkgrey","slategray"), 
                    labels=c("Proposed","NonCons","NoUC","NonCons+NoUC","Q-Learning","IPTW","AIPTW"))
name <- paste("~/Desktop/simu2_average_regret.pdf")
pdf(name, width = 16, height = 8, onefile = TRUE)
p
dev.off()

##########################################################################################
#                                   Figure S5                                            #
##########################################################################################

q_truth <- function(a, x, max_to_min){
  # calculate q-value
  y <- cos(a) + sin(x*a) 
  if (max_to_min){ return(-y) # from max to min
  }else { return(y) }
}

set.seed(102)
n <- 2000; x <- sort(rtruncnorm(n,a=-3,b=3,0,1))
a_truth <- rep(NA, n); y_truth <- rep(NA, n)
for (i in 1:n){
  opt_truth <- optim(par=c(0),fn=q_truth,lower=-3,upper=3,
                     x=x[i],max_to_min=T,method='L-BFGS-B')
  a_truth[i] <- opt_truth$par
  y_truth[i] <- -opt_truth$value # from min to max
}

# Figure S5 - left
df <- data.frame(x=x,y=a_truth)
p <- ggplot(data=df, aes(x,y)) + geom_line(linewidth=2) +
  xlab(expression(X[1])) + ylab("Treatment") + theme_light() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25))
name <- paste("~/Desktop/simu2_optimal_treatment.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S5 - right
df <- data.frame(x=x,y=y_truth)
p <- ggplot(data=df, aes(x,y)) + geom_line(linewidth=2) +
  xlab(expression(X[1])) + ylab("Outcome") + theme_light() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25))
name <- paste("~/Desktop/simu2_optimal_outcome.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()