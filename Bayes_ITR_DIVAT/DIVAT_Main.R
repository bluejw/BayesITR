##########################################################################################
#                                  DIVAT Data Analysis                                   #
##########################################################################################

library(Matrix)
library(MCMCpack)
library(LaplacesDemon)
library(splines2)
library(truncnorm)
library(ggplot2)
library(twangContinuous)

source("MCMC_R_Functions.R") # R functions for MCMC 
sourceCpp("MCMC_Rcpp_Functions.cpp") # Rcpp functions for MCMC 
load("DIVAT_Data.RData") # load DIVAT data
source("Data_Process.R") # R functions for processing DIVAT data

##########################################################################################
#                           Run MCMC and Optim on DIVAT Data                             #
##########################################################################################

source("MCMC_Main.R") # MCMC main function: about 1min to run 10,000 iterations 
load("MCMC_Results.RData") # load posterior samples

# construct credible region
set.seed(123); num_samples <- 500; logll <- rep(NA, num_samples)
for (nit in 1:num_samples){
  logll[nit] <- log_likelihood(post$gamma[nit,], post$alpha[nit,,,], post$beta[nit,,],
                               post$lambda[nit,], post$mu[nit,], post$sigma2[nit,])
}
phi <- cbind(post$gamma,post$alpha[,1,S*2+1,]) # posterior samples
index_credible_region <- order(logll,decreasing=T)[1:(num_samples*0.95)] # credible_region
phi_cr <- phi[index_credible_region,] # posterior samples in credible_region

# data and transformation
x <- Z[,1] # covariates/predictors: DGF
a_spline <- bSpline(A, df=B0); a_trans <- Za
xa_spline <- bSpline(X[,S*2+1], df=B0); xa_trans <- Zx[S*2+1,,]
a_lower <- min(A); a_upper <- max(A)
xa_lower <- min(X[,S*2+1]); xa_upper <- max(X[,S*2+1])

# conservative q-function
q_function_conservative <- function(a, x){
  # transfrom data to spline basis
  a_bs <- predict(a_spline, a) %*% a_trans
  xa_bs <- predict(xa_spline, x*a) %*% xa_trans
  # calculate q-value
  # from minimize creatinine to maximize (negative) creatinine
  y <- -(cbind(a_bs,xa_bs) %*% t(phi_cr)) 
  y_min <- apply(y, 1, min) # lower bound
  return(-y_min) # from max to min for optim
}

# DGF = 0
opt0 <- optim(par=c(0),fn=q_function_conservative,lower=a_lower,upper=a_upper,x=0,method='L-BFGS-B')
opt0$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 7.957523
# DGF = 1
opt1 <- optim(par=c(0),fn=q_function_conservative,lower=a_lower,upper=a_upper,x=1,method='L-BFGS-B')
opt1$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 10.41719

##########################################################################################
#                               Figure 5 and Figure S6                                   #
##########################################################################################

# Figure 5 - left
index <- 1:I; x <- (A[index]*sd(data$tacro_tl))+mean(data$tacro_tl)
y_pred <- matrix(NA, nrow=num_samples, ncol=I)
y_pred <- (post$gamma %*% t(A_bs[index,]))*sd(data$creat)+mean(data$creat)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y_mean)) + 
  geom_line(col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) + 
  geom_jitter(aes(x=x,y=50),shape=4,size=4) + 
  xlab(expression("Dosage of Tacrolimus"~(ng/ml))) + 
  ylab(expression("Creatinine Level"~(umol/l))) + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/real_main.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure 5 - right
s <- 4; index <- 1:I 
x <- (X[index,S+s]*sd(data$tacro_tl))+mean(data$tacro_tl)
y_pred <- matrix(NA, nrow=num_samples, ncol=I)
y_pred <- (post$alpha[,1,S+s,] %*% t(X_bs[index,S+s,]))*sd(data$creat)+mean(data$creat)
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y_mean)) + 
  geom_line(col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  geom_jitter(aes(x=x,y=-22),shape=4,size=4) + 
  xlab(expression("Dosage of Tacrolimus"~(ng/ml) %*% "DGF")) + 
  ylab(expression("Creatinine Level"~(umol/l))) + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/real_interaction.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S6 - left
s <- 2; index <- 1:I; x <- (X[index,s]*sd(data$bmiR))+mean(data$bmiR)
y_pred <- matrix(NA, nrow=num_samples, ncol=I)
y_pred <- (post$alpha[,1,s,] %*% t(X_bs[index,s,]))*sd(data$creat)+mean(data$creat) 
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y_mean)) + 
  geom_line(col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  geom_jitter(aes(x=x,y=95),shape=4,size=4) + 
  xlab(expression("BMI"~(kg/m^2))) + ylab(expression("Creatinine Level"~(umol/l))) + 
  theme_classic() + scale_y_continuous(limits=c(94,160)) + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25))
name <- paste("~/Desktop/real_covariate_bmiR.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S6 - right
s <- 3; index <- 1:I; x <- (X[index,s]*sd(data$ageD))+mean(data$ageD)
y_pred <- matrix(NA, nrow=num_samples, ncol=I)
y_pred <- (post$alpha[,1,s,] %*% t(X_bs[index,s,]))*sd(data$creat)+mean(data$creat) 
y_mean <- colMeans(y_pred)
y_lower <- apply(y_pred, 2, quantile, probs=0.025)
y_upper <- apply(y_pred, 2, quantile, probs=0.975)
df <- data.frame(x=x,y_mean=y_mean,y_lower=y_lower,y_upper=y_upper)
p <- ggplot(data=df, aes(x=x, y=y_mean)) + 
  geom_line(col='steelblue4',linewidth=1.5,linetype='dashed') +
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper),fill='skyblue',alpha=0.25) +
  geom_jitter(aes(x=x,y=100),shape=4,size=4) + 
  xlab("AgeD (years)") + ylab(expression("Creatinine Level"~(umol/l))) + theme_classic() + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25))
name <- paste("~/Desktop/real_covariate_ageD.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

##########################################################################################
#                                    Table 1                                             #
##########################################################################################

# NonCons 
q_function_expected <- function(a, x){
  # transfrom data to spline basis
  a_bs <- predict(a_spline, a) %*% a_trans
  xa_bs <- predict(xa_spline, x*a) %*% xa_trans
  # calculate q-value
  y <- cbind(a_bs,xa_bs) %*% t(phi); y <- rowMeans(y) 
  return(y)
}
# DGF = 0
opt0 <- optim(par=c(5),fn=q_function_expected,lower=a_lower,upper=a_upper,x=0,method='L-BFGS-B')
opt0$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 11.07976
# DGF = 1
opt1 <- optim(par=c(5),fn=q_function_expected,lower=a_lower,upper=a_upper,x=1,method='L-BFGS-B')
opt1$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 33.2

# NoUC
load("MCMC_Results_NoUC.RData") # load posterior samples: NoUC method
# construct credible region
set.seed(123); num_samples <- 500; logll <- rep(NA, num_samples)
for (nit in 1:num_samples){
  logll[nit] <- log_likelihood(post$gamma[nit,], post$alpha[nit,,,], post$beta[nit,,],
                               post$lambda[nit,], post$mu[nit,], post$sigma2[nit,])
}
phi <- cbind(matrix(0,num_samples,B),post$alpha[,1,S*2+1,]) # posterior samples: noUC
index_credible_region <- order(logll,decreasing=T)[1:(num_samples*0.95)] # credible_region
phi_cr <- phi[index_credible_region,] # posterior samples in credible_region
# DGF = 0 
# No treatment main effect -> NA
# DGF = 1
opt1 <- optim(par=c(0),fn=q_function_conservative,lower=a_lower,upper=a_upper,x=1,method='L-BFGS-B')
opt1$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 9.202631

# NonCons+NoUC
# DGF = 0 
# No treatment main effect -> NA
# DGF = 1
opt1 <- optim(par=c(5),fn=q_function_expected,lower=a_lower,upper=a_upper,x=1,method='L-BFGS-B')
opt1$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 33.2

# Q-Learning
phi_all <- lm(Y~A_bs+X_bs[,1,]+X_bs[,2,]+X_bs[,3,]+X_bs[,S+1,]+X_bs[,S+2,]+X_bs[,S+3,]+ 
                X_bs[,S*2+1,]+X_bs[,S*2+2,]+Z[,1]+Z[,2])$coef
phi <- phi_all[c(1,2:4,23:25)]
q_function_expected <- function(a, x){
  # transfrom data to spline basis
  a_bs <- predict(a_spline, a) %*% a_trans
  xa_bs <- predict(xa_spline, x*a) %*% xa_trans
  # calculate q-value
  y <- cbind(rep(1,length(x)),a_bs,xa_bs) %*% phi 
  return(y)
}
# DGF = 0
opt0 <- optim(par=c(5),fn=q_function_expected,lower=a_lower,upper=a_upper,x=0,method='L-BFGS-B')
opt0$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 33.2
# DGF = 1
opt1 <- optim(par=c(5),fn=q_function_expected,lower=a_lower,upper=a_upper,x=1,method='L-BFGS-B')
opt1$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 33.2

# IPTW
df <- data.frame(A=A, X=X[,1:S], Z=Z)
ps.model <- ps.cont(A~X[,1]+X[,2]+X[,3]+Z[,1]+Z[,2], data=df)
weights <- get.weights(ps.model)
outcome.model <- lm(Y~A_bs+X_bs[,S+1,]+X_bs[,S+2,]+X_bs[,S+3,]+
                      X_bs[,S*2+1,]+X_bs[,S*2+2,], weights=weights)
phi_all <- outcome.model$coef
phi <- phi_all[c(1,2:4,14:16)]
q_function_expected <- function(a, x){
  # transfrom data to spline basis
  a_bs <- predict(a_spline, a) %*% a_trans
  xa_bs <- predict(xa_spline, x*a) %*% xa_trans
  # calculate q-value
  y <- cbind(rep(1,length(x)),a_bs,xa_bs) %*% phi 
  return(y)
}
# DGF = 0
opt0 <- optim(par=c(5),fn=q_function_expected,lower=a_lower,upper=a_upper,x=0,method='L-BFGS-B')
opt0$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 18.70545
# DGF = 1
opt1 <- optim(par=c(5),fn=q_function_expected,lower=a_lower,upper=a_upper,x=1,method='L-BFGS-B')
opt1$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 33.2

# AIPTW
df <- data.frame(A=A, X=X[,1:S], Z=Z)
ps.model <- ps.cont(A~X[,1]+X[,2]+X[,3]+Z[,1]+Z[,2], data=df)
weights <- get.weights(ps.model)
outcome.model <- lm(Y~A_bs+X_bs[,1,]+X_bs[,2,]+X_bs[,3,]+X_bs[,S+1,]+X_bs[,S+2,]+X_bs[,S+3,]+
                      X_bs[,S*2+1,]+X_bs[,S*2+2,]+Z[,1]+Z[,2], weights=weights)
phi_all <- outcome.model$coef
phi <- phi_all[c(1,2:4,23:25)]
q_function_expected <- function(a, x){
  # transfrom data to spline basis
  a_bs <- predict(a_spline, a) %*% a_trans
  xa_bs <- predict(xa_spline, x*a) %*% xa_trans
  # calculate q-value
  y <- cbind(rep(1,length(x)),a_bs,xa_bs) %*% phi 
  return(y)
}
# DGF = 0
opt0 <- optim(par=c(5),fn=q_function_expected,lower=a_lower,upper=a_upper,x=0,method='L-BFGS-B')
opt0$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 33.2
# DGF = 1
opt1 <- optim(par=c(5),fn=q_function_expected,lower=a_lower,upper=a_upper,x=1,method='L-BFGS-B')
opt1$par*sd(data$tacro_tl)+mean(data$tacro_tl) # 33.2