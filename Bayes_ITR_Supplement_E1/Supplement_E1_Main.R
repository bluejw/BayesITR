##########################################################################################
#                            Simulation Study: Supplement E1                             #
##########################################################################################

library(Matrix)
library(MCMCpack)
library(LaplacesDemon)
library(splines2)
library(ggplot2)

source("MCMC_R_Functions.R") # R functions for MCMC 
sourceCpp("MCMC_Rcpp_Functions.cpp") # Rcpp functions for MCMC 

##########################################################################################
#                           Run MCMC on Simulated Data                                   #
##########################################################################################

# Scenario I(a) + Non-Gaussian Unmeasured Confounder
source("Simu_Data_Generate_I(a).R") # generate simulation truths: Scenario I(a) 
source("MCMC_Main.R") # MCMC main function: less than 2min to run 10,000 iterations 

# Scenario I(b) + Non-Gaussian Unmeasured Confounder
source("Simu_Data_Generate_I(b).R") # generate simulation truths: Scenario I(b)
source("MCMC_Main.R") # MCMC main function: less than 2min to run 10,000 iterations 

##########################################################################################
#                                  Figure E1                                             #
##########################################################################################

num_samples <- 500
source("Simu_Data_Generate_I(a).R") # generate simulation truths: Scenario I(a)
load("MCMC_Results_I(a)_Mode1.RData") # posterior samples: Scenario I(a), Mode 1

# Figure E1(a) - left
index <- !(1:I %in% c(order(A)[1:5],which.max(A)))
x <- A[index]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-6)
y_pred <- post$gamma %*% t(A_bs[index,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) +
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) + 
  xlab(" ") + ylab(" ") + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu1.1_positive_main_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure E1(b) - left
index <- !(1:I %in% c(order(X[,S+1])[1:2],which.max(X[,S+1])))
x <- X[index,S+1]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-3)
y_pred <- post$alpha[,1,S+1,] %*% t(X_bs[index,S+1,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) + 
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  xlab(" ") + ylab(" ") + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu1.1_positive_interaction_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure E1(c) - left
index <- !(1:I %in% c(order(X[,1])[1:2],which.max(X[,1])))
x <- X[index,1]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-3)
y_pred <- post$alpha[,1,1,] %*% t(X_bs[index,1,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) +
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  xlab(" ") + ylab(" ") + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25))
name <- paste("~/Desktop/simu1.1_positive_covariate_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

source("Simu_Data_Generate_I(a).R") # generate simulation truths: Scenario I(a)
load("MCMC_Results_I(a)_Mode2.RData") # posterior samples: Scenario I(a), Mode 2

# Figure E1(a) - right
index <- !(1:I %in% c(order(A)[1:5],which.max(A)))
x <- A[index]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-6)
y_pred <- post$gamma %*% t(A_bs[index,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) +
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) + 
  xlab(" ") + ylab(" ") + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu1.1_negative_main_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure E1(b) - right
index <- !(1:I %in% c(order(X[,S+1])[1:2],which.max(X[,S+1])))
x <- X[index,S+1]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-3)
y_pred <- post$alpha[,1,S+1,] %*% t(X_bs[index,S+1,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) + 
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  xlab(" ") + ylab(" ") + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu1.1_negative_interaction_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure E1(c) - right
index <- !(1:I %in% c(order(X[,1])[1:2],which.max(X[,1])))
x <- X[index,1]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-3)
y_pred <- post$alpha[,1,1,] %*% t(X_bs[index,1,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) +
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  xlab(" ") + ylab(" ") + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25))
name <- paste("~/Desktop/simu1.1_negative_covariate_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

##########################################################################################
#                                  Figure E2                                             #
##########################################################################################

num_samples <- 500
source("Simu_Data_Generate_I(b).R") # generate simulation truths: Scenario I(b)
load("MCMC_Results_I(b)_Mode1.RData") # posterior samples: Scenario I(b), Mode 1

# Figure E2(a) - left
index <- !(1:I %in% c(order(A)[1:5],which.max(A)))
x <- A[index]; y <- 0.1*exp(x)
y_pred <- matrix(NA, nrow=num_samples, ncol=I-6)
y_pred <- post$gamma %*% t(A_bs[index,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) +
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) + 
  xlab(" ") + ylab(" ") + theme_classic() + 
  scale_x_continuous(breaks=c(-2,0,2,4)) + scale_y_continuous(breaks=c(0,3,6,9,12)) + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu1.2_positive_main_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure E2(b) - left
index <- !(1:I %in% c(order(X[,S+1])[1:5],which.max(X[,S+1])))
x <- X[index,S+1]; y <- 0.2*(x^2)
y_pred <- matrix(NA, nrow=num_samples, ncol=I-6)
y_pred <- post$alpha[,1,S+1,] %*% t(X_bs[index,S+1,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) + 
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  xlab(" ") + ylab(" ") + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu1.2_positive_interaction_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure E2(c) - left
index <- !(1:I %in% c(order(X[,1])[1:5],which.max(X[,1])))
x <- X[index,1]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-6)
y_pred <- post$alpha[,1,1,] %*% t(X_bs[index,1,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) + 
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  xlab(" ") + ylab(" ") + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25))
name <- paste("~/Desktop/simu1.2_positive_covariate_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

num_samples <- 500
source("Simu_Data_Generate_I(b).R") # generate simulation truths: Scenario I(b)
load("MCMC_Results_I(b)_Mode2.RData") # posterior samples: Scenario I(b), Mode 2

# Figure E2(a) - right
index <- !(1:I %in% c(order(A)[1:5],which.max(A)))
x <- A[index]; y <- 0.1*exp(x)
y_pred <- matrix(NA, nrow=num_samples, ncol=I-6)
y_pred <- post$gamma %*% t(A_bs[index,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) +
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) + 
  xlab(" ") + ylab(" ") + theme_classic() + 
  scale_x_continuous(breaks=c(-2,0,2,4)) + scale_y_continuous(breaks=c(0,3,6,9,12)) + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu1.2_negative_main_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure E2(b) - right
index <- !(1:I %in% c(order(X[,S+1])[1:5],which.max(X[,S+1])))
x <- X[index,S+1]; y <- 0.2*(x^2)
y_pred <- matrix(NA, nrow=num_samples, ncol=I-6)
y_pred <- post$alpha[,1,S+1,] %*% t(X_bs[index,S+1,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) + 
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  xlab(" ") + ylab(" ") + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu1.2_negative_interaction_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure E2(c) - right
index <- !(1:I %in% c(order(X[,1])[1:5],which.max(X[,1])))
x <- X[index,1]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-6)
y_pred <- post$alpha[,1,1,] %*% t(X_bs[index,1,]) + mean(y)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y=y,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y)) + geom_line(linewidth=1.5) + 
  geom_line(aes(x=x,y=y_mean),col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  xlab(" ") + ylab(" ") + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25))
name <- paste("~/Desktop/simu1.2_negative_covariate_nonnormal_U.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()