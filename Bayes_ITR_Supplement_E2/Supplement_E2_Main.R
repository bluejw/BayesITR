##########################################################################################
#                            Simulation Study: Supplement E2                             #
##########################################################################################

library(Matrix)
library(MCMCpack)
library(LaplacesDemon)
library(splines2)

source("MCMC_R_Functions.R") # R functions for MCMC 
sourceCpp("MCMC_Rcpp_Functions.cpp") # Rcpp functions for MCMC 

##########################################################################################
#                           Run MCMC on Simulated Data                                   #
##########################################################################################

# Generate simulation truths
source("Simu_Data_Generate_I(a).R") 

# Setup hyperparameters for the spike-and-slab priors
# nu_0 <- 1e-4; a_nu <- 5; b_nu <- 50
nu_0 <- 2.5e-4; a_nu <- 5; b_nu <- 50
# nu_0 <- 5e-4; a_nu <- 5; b_nu <- 50

# nu_0 <- 2.5e-4; a_nu <- 5; b_nu <- 25
nu_0 <- 2.5e-4; a_nu <- 5; b_nu <- 50
# nu_0 <- 2.5e-4; a_nu <- 5; b_nu <- 100

# MCMC main function
source("MCMC_Main.R")

##########################################################################################
#                                   Table E3                                             #
##########################################################################################

# nu_0 <- 1e-4; a_nu <- 5; b_nu <- 50
load("MCMC_Results_S1.RData") # posterior samples
# gamma, delta_1, delta_2, beta_1^Y, beta_2^Y, beta_1^A, beta_2^A
mean(post$xi_gamma==1) # 1
mean(post$xi_alpha[,1,3]==1) # 1
mean(post$xi_alpha[,1,4]==1) # 0.11
mean(post$xi_alpha[,1,1]==1) # 1
mean(post$xi_alpha[,1,2]==1) # 0.17
mean(post$xi_alpha[,2,1]==1) # 1
mean(post$xi_alpha[,2,2]==1) # 0.05

# nu_0 <- 2.5e-4; a_nu <- 5; b_nu <- 50
load("MCMC_Results_S2.RData") # posterior samples
# gamma, delta_1, delta_2, beta_1^Y, beta_2^Y, beta_1^A, beta_2^A
mean(post$xi_gamma==1) # 1
mean(post$xi_alpha[,1,3]==1) # 1
mean(post$xi_alpha[,1,4]==1) # 0.09
mean(post$xi_alpha[,1,1]==1) # 1
mean(post$xi_alpha[,1,2]==1) # 0.14
mean(post$xi_alpha[,2,1]==1) # 1
mean(post$xi_alpha[,2,2]==1) # 0.05

# nu_0 <- 5e-4; a_nu <- 5; b_nu <- 50
load("MCMC_Results_S3.RData") # posterior samples
# gamma, delta_1, delta_2, beta_1^Y, beta_2^Y, beta_1^A, beta_2^A
mean(post$xi_gamma==1) # 1
mean(post$xi_alpha[,1,3]==1) # 1
mean(post$xi_alpha[,1,4]==1) # 0.09
mean(post$xi_alpha[,1,1]==1) # 1
mean(post$xi_alpha[,1,2]==1) # 0.09
mean(post$xi_alpha[,2,1]==1) # 1
mean(post$xi_alpha[,2,2]==1) # 0.03

# nu_0 <- 2.5e-4; a_nu <- 5; b_nu <- 25
load("MCMC_Results_S4.RData") # posterior samples
# gamma, delta_1, delta_2, beta_1^Y, beta_2^Y, beta_1^A, beta_2^A
mean(post$xi_gamma==1) # 1
mean(post$xi_alpha[,1,3]==1) # 1
mean(post$xi_alpha[,1,4]==1) # 0.10
mean(post$xi_alpha[,1,1]==1) # 1
mean(post$xi_alpha[,1,2]==1) # 0.18
mean(post$xi_alpha[,2,1]==1) # 1
mean(post$xi_alpha[,2,2]==1) # 0.08

# nu_0 <- 2.5e-4; a_nu <- 5; b_nu <- 50
load("MCMC_Results_S5.RData") # posterior samples
# gamma, delta_1, delta_2, beta_1^Y, beta_2^Y, beta_1^A, beta_2^A
mean(post$xi_gamma==1) # 1
mean(post$xi_alpha[,1,3]==1) # 1
mean(post$xi_alpha[,1,4]==1) # 0.09
mean(post$xi_alpha[,1,1]==1) # 1
mean(post$xi_alpha[,1,2]==1) # 0.14
mean(post$xi_alpha[,2,1]==1) # 1
mean(post$xi_alpha[,2,2]==1) # 0.05

# nu_0 <- 2.5e-4; a_nu <- 5; b_nu <- 100
load("MCMC_Results_S6.RData") # posterior samples
# gamma, delta_1, delta_2, beta_1^Y, beta_2^Y, beta_1^A, beta_2^A
mean(post$xi_gamma==1) # 1
mean(post$xi_alpha[,1,3]==1) # 1
mean(post$xi_alpha[,1,4]==1) # 0.06
mean(post$xi_alpha[,1,1]==1) # 1
mean(post$xi_alpha[,1,2]==1) # 0.09
mean(post$xi_alpha[,2,1]==1) # 1
mean(post$xi_alpha[,2,2]==1) # 0.04
