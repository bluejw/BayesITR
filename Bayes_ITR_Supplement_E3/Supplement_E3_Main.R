##########################################################################################
#                            Simulation Study: Supplement E3                             #
##########################################################################################

library(Matrix)
library(MCMCpack)
library(LaplacesDemon)
library(splines2)
library(ggplot2)

source("MCMC_R_Functions.R") # R functions for MCMC 
sourceCpp("MCMC_Rcpp_Functions.cpp") # Rcpp functions for MCMC 

##########################################################################################
#                          Run MCMC and SGD on Simulated Data                            #
##########################################################################################

# Simulation setups
noUC <- F # with unmeasured confounders (proposed method)
# Probability of selected data outside the region
prob <- 1 # Scenario II(a), balanced data
# Average regret
num_seed <- 100 # number of seeds/experiments
regret_conservative <- rep(NA, num_seed)
# Generate testing data 
set.seed(102); I_test <- 500
x_test <- rtruncnorm(I_test,a=-2,b=2,0,1)

for (seed_index in 1:num_seed){
  
  # MCMC
  set.seed(seed_index) # set the seed to generate different simulated datasets
  source("Simu_Data_Generate.R") # generate simulation truths: Scenario II
  source("Select_Partial_Data.R") # select partial data to determine sub-scenarios
  source("MCMC_Main.R") # MCMC main function: less than 2min to run 10,000 iterations 
  
  # construct credible region
  # alpha <- 0.9 # confidence level 
  alpha <- 0.95 # confidence level 
  # alpha <- 0.99 # confidence level 
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
  # std_sgd <- 0.1 # standard deviation of the policy function
  std_sgd <- 0.2 # standard deviation of the policy function
  # std_sgd <- 0.3 # standard deviation of the policy function
  # step_size <- 7.5e-5 # step size for sgd
  step_size <- 1e-4 # step size for sgd
  # step_size <- 2.5e-4 # step size for sgd
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
  
  # Conservative reward 
  source("SGD_Main_Conservative.R") # SGD main function
  theta_star_conservative <- colMeans(sgd$theta[seq((Niter/2+1),Niter,by=1),])
  a_est_conservative <- policy_function(x_test, I_test, theta_star_conservative, 0) # optimal action
  y_est_conservative_truth <- q_truth(a_est_conservative, x_test, max_to_min=F) # true value of the optimal action 
  regret_conservative[seed_index] <- mean(y_truth_est-y_est_conservative_truth)
}

##########################################################################################
#                                  Figure E4                                             #
##########################################################################################

# Figure E4(a)
load("SGD_Results_Alpha.RData")
df <- data.frame(x=c(1,2,3),
                 y=c(res$regret_alpha99,res$regret_alpha95,res$regret_alpha90),
                 z=factor(c("alpha=0.01","alpha=0.05","alpha=0.10"),
                          levels=c("alpha=0.01","alpha=0.05","alpha=0.10")))
p <- ggplot(data=df, aes(x=x, y=y, fill=z)) + 
  stat_boxplot(geom ='errorbar', width=0.35, position=position_dodge(0.75)) + 
  geom_boxplot() + xlab(" ") + ylab("Average Regret") + theme_light() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_blank()) + theme(axis.text.y = element_text(size=25)) + 
  theme(legend.title=element_text(size=25), legend.text=element_text(size=25)) +
  scale_fill_manual(name=expression(alpha), values=c("white","lightblue1","skyblue3"), 
                    labels=c("0.01","0.05","0.1"))
name <- paste("~/Desktop/simu2_average_regret_alpha.pdf")
pdf(name, width = 16, height = 6, onefile = TRUE)
p
dev.off()

# Figure E4(b)
load("SGD_Results_Kappa.RData")
df <- data.frame(x=c(1,2,3),
                 y=c(res$regret_step1,res$regret_step2,res$regret_step3),
                 z=factor(c("step=7.5e-5","step=1e-4","step=2.5e-4"),
                          levels=c("step=7.5e-5","step=1e-4","step=2.5e-4")))
p <- ggplot(data=df, aes(x=x, y=y, fill=z)) + 
  stat_boxplot(geom ='errorbar', width=0.35, position=position_dodge(0.75)) + 
  geom_boxplot() + xlab(" ") + ylab("Average Regret") + theme_light() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_blank()) + theme(axis.text.y = element_text(size=25)) + 
  theme(legend.title=element_text(size=25), legend.text=element_text(size=25)) + 
  scale_fill_manual(name=expression(kappa[t]), values=c("white","lightblue1","skyblue3"), 
                    labels=c("7.5e-5","1e-4","2.5e-4"))
name <- paste("~/Desktop/simu2_average_regret_step.pdf")
pdf(name, width = 16, height = 6, onefile = TRUE)
p
dev.off()

# Figure E4(c)
load("SGD_Results_Sigma.RData")
df <- data.frame(x=c(1,2,3),
                 y=c(res$regret_std_sgd1,res$regret_std_sgd2,res$regret_std_sgd3),
                 z=factor(c("std_sg=0.1","std_sg=0.2","std_sg=0.3"),
                          levels=c("std_sg=0.1","std_sg=0.2","std_sg=0.3")))
p <- ggplot(data=df, aes(x=x, y=y, fill=z)) + 
  stat_boxplot(geom ='errorbar', width=0.35, position=position_dodge(0.75)) + 
  geom_boxplot() + xlab(" ") + ylab("Average Regret") + theme_light() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_blank()) + theme(axis.text.y = element_text(size=25)) + 
  theme(legend.title=element_text(size=25), legend.text=element_text(size=25)) + 
  scale_fill_manual(name=expression(sigma[pi]), values=c("white","lightblue1","skyblue3"), 
                    labels=c("0.1","0.2","0.3"))
name <- paste("~/Desktop/simu2_average_regret_std_sgd.pdf")
pdf(name, width = 16, height = 6, onefile = TRUE)
p
dev.off()