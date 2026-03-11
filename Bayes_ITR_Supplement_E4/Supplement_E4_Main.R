##########################################################################################
#                            Simulation Study: Supplement E4                             #
##########################################################################################

library(Matrix)
library(MCMCpack)
library(LaplacesDemon)
library(splines2)
library(truncnorm)

source("MCMC_R_Functions.R") # R functions for MCMC 
sourceCpp("MCMC_Rcpp_Functions.cpp") # Rcpp functions for MCMC 

##########################################################################################
#                          Run MCMC and SGD on Simulated Data                            #
##########################################################################################

# Simulation setups
noUC <- F # with unmeasured confounders

# Probability of selected data outside the region
prob <- 1 # Scenario II(a), balanced data
# prob <- 0.5 # Scenario II(b), mild unbalanced data
# prob <- 0.1 Scenarios II(c) and II(d), severe unbalanced data

# average regret
num_seed <- 100 # number of seeds/experiments
regret_conservative <- rep(NA, num_seed) # proposed method
regret_conservative_mean <- rep(NA, num_seed) # alternative method

# generate testing data
set.seed(102); I_test <- 500
# x_test <- rtruncnorm(I_test,a=-2,b=2,0,1) # Scenario II(a-c)
x_test <- rnorm(I_test,1,1) # Scenario II(d)

# Start of running repeated experiments: Proposed and Alternative
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
  
  # Conservative reward (Proposed)
  source("SGD_Main_Conservative.R") # SGD main function: about 10s to run 1,000 iterations 
  theta_star_conservative <- colMeans(sgd$theta[seq((Niter/2+1),Niter,by=1),])
  a_est_conservative <- policy_function(x_test, I_test, theta_star_conservative, 0) # optimal action
  y_est_conservative_truth <- q_truth(a_est_conservative, x_test, max_to_min=F) # true value of the optimal action 
  regret_conservative[seed_index] <- mean(y_truth_est-y_est_conservative_truth)
  
  # Conservative mean reward (Alternative)
  source("SGD_Main_Conservative_Mean.R") # SGD main function: about 10s to run 1,000 iterations 
  theta_star_conservative_mean <- colMeans(sgd$theta[seq((Niter/2+1),Niter,by=1),])
  a_est_conservative_mean <- policy_function(x_test, I_test, theta_star_conservative_mean, 0) # optimal action
  y_est_conservative_mean_truth <- q_truth(a_est_conservative_mean, x_test, max_to_min=F) # true value of the optimal action 
  regret_conservative_mean[seed_index] <- mean(y_truth_est-y_est_conservative_mean_truth)
}

##########################################################################################
#                                   Figure E5                                            #
##########################################################################################

load("SGD_Results_II(a).RData") # simulation results: Scenario II(a)
load("SGD_Results_II(b).RData") # simulation results: Scenario II(b)
load("SGD_Results_II(c).RData") # simulation results: Scenario II(c)
load("SGD_Results_II(d).RData") # simulation results: Scenario II(d)

df <- data.frame(x=c(rep("Scenario II(a)",num_seed*2),
                     rep("Scenario II(b)",num_seed*2),
                     rep("Scenario II(c)",num_seed*2),
                     rep("Scenario II(d)",num_seed*2)),
                 y=c(res10$regret_conservative, res10$regret_conservative_mean,
                     res5$regret_conservative, res5$regret_conservative_mean,
                     res1$regret_conservative, res1$regret_conservative_mean,
                     res$regret_conservative, res$regret_conservative_mean),
                 z=factor(rep(c(rep("Proposed",num_seed), rep("Alternative",num_seed)),4),
                          levels=c("Proposed","Alternative")))
p <- ggplot(data=df, aes(x=x, y=y, fill=z)) +
  stat_boxplot(geom ='errorbar', width=0.35, position=position_dodge(0.75)) +
  geom_boxplot() + xlab(" ") + ylab("Average Regret") + theme_light() +
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) +
  theme(legend.title=element_text(size=25), legend.text=element_text(size=25)) +
  scale_y_continuous(limits=c(0,0.22)) +
  scale_fill_manual(name="Method", values=c("white","steelblue"),
                    labels=c("Proposed","Alternative"))
name <- paste("~/Desktop/simu2_method_compare_average_regret.pdf")
pdf(name, width = 16, height = 8, onefile = TRUE)
p
dev.off()
