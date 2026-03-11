##########################################################################################
#                            Simulation Study: Supplement C                              #
##########################################################################################

library(Matrix)
library(MCMCpack)
library(LaplacesDemon)
library(splines2)
library(ggplot2)
library(loo)

source("MCMC_R_Functions.R") # R functions for MCMC 
sourceCpp("MCMC_Rcpp_Functions.cpp") # Rcpp functions for MCMC 
source("MCMC_Main.R") # MCMC main function

##########################################################################################
#                                   Figure C1                                            #
##########################################################################################

# calculate elpd values
num_samples <- 500; num_df <- 8; num_data <- 2000
logll <- array(NA, dim=c(num_df, num_samples, num_data))
elpd <- rep(NA, num_df)

for (B0_index in 3:10){
  
  source("Simu_Data_Generate.R")
  file_name <- paste("MCMC_Results_B0=", B0_index, ".RData", sep="")
  load(file=file_name)
  
  B0 <- B0_index # B0 spline expansion
  B <- B0 - 1 # sum-to-zero identifiability constraint
  
  A_bs0 <- bSpline(A, df=B0)  # b-spline expansion for treatment
  A_bs <- array(NA, dim=c(I, B)) # b-spline expansion with constraint
  Za <- array(NA, dim=c(B0, B)) # transformation matrix
  Ca <- t(rep(1,I)) %*% A_bs0
  qr_decomp_a <- qr(t(Ca))
  Za <- qr.Q(qr_decomp_a, complete=T)[,2:B0]
  A_bs <- A_bs0 %*% Za
  
  X_bs0 <- array(NA, dim=c(I, S_sum, B0)) # b-spline expansion for continuous covariates
  for (s in 1:S_sum){ X_bs0[,s,] <- bSpline(X[,s], df=B0) }
  X_bs <- array(NA, dim=c(I, S_sum, B)) # b-spline expansion with constraint
  Zx <- array(NA, dim=c(S_sum, B0, B)) # transformation matrix
  for (s in 1:S_sum){
    Cx <- t(rep(1,I)) %*% X_bs0[,s,]
    qr_decomp_x <- qr(t(Cx))
    Zx[s,,] <- qr.Q(qr_decomp_x, complete=T)[,2:B0]
    X_bs[,s,] <- X_bs0[,s,] %*% Zx[s,,]
  }
  
  # calculate log-likelihood
  for (nit in 1:num_samples){
    print(paste(B0,nit,sep=","))
    logll[B0-2,nit,] <- log_likelihood(post$gamma[nit,], post$alpha[nit,,,], post$beta[nit,,],
                                       post$lambda[nit,], post$mu[nit,], post$sigma2[nit,])
  }

  # calculate elpd
  elpd[B0-2] <- loo(logll[B0-2,,])$elpd
}

df <- data.frame(x=3:10, y=elpd)
p <- ggplot(data=df, aes(x=x,y=y)) + geom_line(linewidth=1.5) + geom_point(size=4) + 
  geom_point(aes(x=which.max(elpd)+2,y=max(elpd)),col='red',size=4) + 
  xlab("Degrees of Freedom") + ylab("elpd") + theme_light() + 
  scale_x_continuous(labels=as.character(3:10), breaks=3:10) + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu0_optimal_df.pdf")
pdf(name, width = 16, height = 8, onefile = TRUE)
p
dev.off()

##########################################################################################
#                                   Figure S2                                            #
##########################################################################################

load("ELPD_Values.RData") # elpd values

# Figure S2(a)
df <- data.frame(x=3:10, y=elpd$Simu1a)
p <- ggplot(data=df, aes(x=x,y=y)) + geom_line(linewidth=1.5) + geom_point(size=4) + 
  geom_point(aes(x=which.max(elpd$Simu1a)+2,y=max(elpd$Simu1a)),col='red',size=4) + 
  xlab("Degrees of Freedom") + ylab("elpd") + theme_light() + 
  scale_x_continuous(labels=as.character(3:10), breaks=3:10) + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu1.1_optimal_df.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S2(b)
df <- data.frame(x=3:10, y=elpd$Simu1b)
p <- ggplot(data=df, aes(x=x,y=y)) + geom_line(linewidth=1.5) + geom_point(size=4) + 
  geom_point(aes(x=which.max(elpd$Simu1b)+2,y=max(elpd$Simu1b)),col='red',size=4) + 
  xlab("Degrees of Freedom") + ylab("elpd") + theme_light() + 
  scale_x_continuous(labels=as.character(3:10), breaks=3:10) + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu1.2_optimal_df.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S2(c)
df <- data.frame(x=3:10, y=elpd$Simu2)
p <- ggplot(data=df, aes(x=x,y=y)) + geom_line(linewidth=1.5) + geom_point(size=4) + 
  geom_point(aes(x=which.max(elpd$Simu2)+2,y=max(elpd$Simu2)),col='red',size=4) + 
  xlab("Degrees of Freedom") + ylab("elpd") + theme_light() + 
  scale_x_continuous(labels=as.character(3:10), breaks=3:10) + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/simu2_optimal_df.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()

# Figure S2(d)
df <- data.frame(x=3:10, y=elpd$DIVAT)
p <- ggplot(data=df, aes(x=x,y=y)) + geom_line(linewidth=1.5) + geom_point(size=4) + 
  geom_point(aes(x=which.max(elpd$DIVAT)+2,y=max(elpd$DIVAT)),col='red',size=4) + 
  xlab("Degrees of Freedom") + ylab("elpd") + theme_light() + 
  scale_x_continuous(labels=as.character(3:10), breaks=3:10) + 
  theme(axis.title.x=element_text(size=25)) + theme(axis.title.y = element_text(size=25)) +
  theme(axis.text.x=element_text(size=25)) + theme(axis.text.y = element_text(size=25)) 
name <- paste("~/Desktop/real_optimal_df.pdf")
pdf(name, width = 10, height = 8, onefile = TRUE)
p
dev.off()