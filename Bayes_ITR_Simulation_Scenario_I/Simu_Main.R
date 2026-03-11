##########################################################################################
#                             Simulation Study: Scenario I                               #
##########################################################################################

library(Rcpp)
library(RcppArmadillo)
library(Matrix)
library(MCMCpack)
library(LaplacesDemon)
library(splines2)
library(ggplot2)
library(gridExtra)

source("MCMC_R_Functions.R") # R functions for MCMC 
sourceCpp("MCMC_Rcpp_Functions.cpp") # Rcpp functions for MCMC 

##########################################################################################
#                           Run MCMC on Simulated Data                                   #
##########################################################################################

# Scenario I(a)
source("Simu_Data_Generate_I(a).R") # generate simulation truths: Scenario I(a)
source("MCMC_Main.R") # MCMC main function: less than 2min to run 10,000 iterations 

# Scenario I(b)
source("Simu_Data_Generate_I(b).R") # generate simulation truths: Scenario I(b)
source("MCMC_Main.R") # MCMC main function: less than 2min to run 10,000 iterations 

##########################################################################################
#                               Figure 3 and Figure S3                                   #
##########################################################################################

num_samples <- 500
source("Simu_Data_Generate_I(a).R") # generate simulation truths: Scenario I(a)
load("MCMC_Results_I(a)_Mode1.RData") # posterior samples: Scenario I(a), Mode 1

# Figure 3 - left
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
name <- paste("~/Desktop/simu1.1_positive_main.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S3(a) - left
index <- !(1:I %in% c(which.min(X[,S+1]),which.max(X[,S+1])))
x <- X[index,S+1]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-2)
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
name <- paste("~/Desktop/simu1.1_positive_interaction.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S3(b) - left
index <- !(1:I %in% c(which.min(X[,1]),which.max(X[,1])))
x <- X[index,1]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-2)
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
name <- paste("~/Desktop/simu1.1_positive_covariate.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

source("Simu_Data_Generate_I(a).R") # generate simulation truths: Scenario I(a)
load("MCMC_Results_I(a)_Mode2.RData") # posterior samples: Scenario I(a), Mode 2

# Figure 3 - right
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
name <- paste("~/Desktop/simu1.1_negative_main.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S3(a) - right
index <- !(1:I %in% c(which.min(X[,S+1]),which.max(X[,S+1])))
x <- X[index,S+1]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-2)
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
name <- paste("~/Desktop/simu1.1_negative_interaction.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S3(b) - right
index <- !(1:I %in% c(which.min(X[,1]),which.max(X[,1])))
x <- X[index,1]; y <- 0.5*x
y_pred <- matrix(NA, nrow=num_samples, ncol=I-2)
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
name <- paste("~/Desktop/simu1.1_negative_covariate.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

##########################################################################################
#                                   Figure S4                                            #
##########################################################################################

source("Simu_Data_Generate_I(b).R") # generate simulation truths: Scenario I(b)
load("MCMC_Results_I(b)_Mode1.RData") # posterior samples: Scenario I(b), Mode 1

# Figure S4(a) - left
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
name <- paste("~/Desktop/simu1.2_positive_main.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S4(b) - left
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
name <- paste("~/Desktop/simu1.2_positive_interaction.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S4(c) - left
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
name <- paste("~/Desktop/simu1.2_positive_covariate.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

source("Simu_Data_Generate_I(b).R") # generate simulation truths: Scenario I(b)
load("MCMC_Results_I(b)_Mode2.RData") # posterior samples: Scenario I(b), Mode 2

# Figure S4(a) - right
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
name <- paste("~/Desktop/simu1.2_negative_main.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S4(b) - right
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
name <- paste("~/Desktop/simu1.2_negative_interaction.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S4(c) - right
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
name <- paste("~/Desktop/simu1.2_negative_covariate.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

##########################################################################################
#                                   Figure S1                                            #
##########################################################################################

source("Simu_Data_Generate_I(a).R") # generate simulation truths: Scenario I(a)
load("MCMC_Results_I(a)_Mode1.RData") # posterior samples: Scenario I(a), Mode 1

# Figure S1(a)
theme_update(plot.title = element_text(hjust = 0.5, size = 25))
df11 <- data.frame(iteration = seq(1, 500, by=1), value = post$gamma[,1])
p11 <- ggplot(df11, aes(iteration, value)) + geom_line() +  
  labs(title=expression(tilde(gamma)['1'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df12 <- data.frame(iteration = seq(1, 500, by=1), value = post$gamma[,2])
p12 <- ggplot(df12, aes(iteration, value)) + geom_line() +  
  labs(title=expression(tilde(gamma)['2'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df13 <- data.frame(iteration = seq(1, 500, by=1), value = post$alpha[,1,3,2])
p13 <- ggplot(df13, aes(iteration, value)) + geom_line() +  
  labs(title=expression(tilde(delta)['2,1'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df21 <- data.frame(iteration = seq(1, 500, by=1), value = post$alpha[,1,3,3])
p21 <- ggplot(df21, aes(iteration, value)) + geom_line() +  
  labs(title=expression(tilde(delta)['3,1'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df22 <- data.frame(iteration = seq(1, 500, by=1), value = abs(post$lambda[,2]))
p22 <- ggplot(df22, aes(iteration, value)) + geom_line() +  
  labs(title=expression(lambda['A'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df23 <- data.frame(iteration = seq(1, 500, by=1), value = abs(post$sigma2[,1]))
p23 <- ggplot(df23, aes(iteration, value)) + geom_line() +  
  labs(title=expression(tilde(sigma)['Y']^{2})) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
name <- paste("~/Desktop/trace_plot_simu1.1.pdf")
pdf(name, width = 24, height = 12, onefile = TRUE)
grid.arrange(p11, p12, p13, p21, p22, p23, widths = c(1, 1, 1),
             layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6)))
dev.off()

source("Simu_Data_Generate_I(b).R") # generate simulation truths: Scenario I(b)
load("MCMC_Results_I(b)_Mode1.RData") # posterior samples: Scenario I(b), Mode 1

# Figure S1(b)
theme_update(plot.title = element_text(hjust = 0.5, size = 25))
df11 <- data.frame(iteration = seq(1, 500, by=1), value = post$gamma[,1])
p11 <- ggplot(df11, aes(iteration, value)) + geom_line() +  
  labs(title=expression(tilde(gamma)['1'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df12 <- data.frame(iteration = seq(1, 500, by=1), value = post$gamma[,2])
p12 <- ggplot(df12, aes(iteration, value)) + geom_line() +  
  labs(title=expression(tilde(gamma)['2'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df13 <- data.frame(iteration = seq(1, 500, by=1), value = post$alpha[,1,3,2])
p13 <- ggplot(df13, aes(iteration, value)) + geom_line() +  
  labs(title=expression(tilde(delta)['2,1'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df21 <- data.frame(iteration = seq(1, 500, by=1), value = post$alpha[,1,3,3])
p21 <- ggplot(df21, aes(iteration, value)) + geom_line() +  
  labs(title=expression(tilde(delta)['3,1'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df22 <- data.frame(iteration = seq(1, 500, by=1), value = abs(post$lambda[,2]))
p22 <- ggplot(df22, aes(iteration, value)) + geom_line() +  
  labs(title=expression(lambda['A'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df23 <- data.frame(iteration = seq(1, 500, by=1), value = abs(post$sigma2[,1]))
p23 <- ggplot(df23, aes(iteration, value)) + geom_line() +  
  labs(title=expression(tilde(sigma)['Y']^{2})) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
name <- paste("~/Desktop/trace_plot_simu1.2.pdf")
pdf(name, width = 24, height = 12, onefile = TRUE)
grid.arrange(p11, p12, p13, p21, p22, p23, widths = c(1, 1, 1),
             layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6)))
dev.off()
