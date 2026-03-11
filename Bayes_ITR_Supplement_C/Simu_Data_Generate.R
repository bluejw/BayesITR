
# start to generate simulated data 
set.seed(102) 

I <- 2000 # number of data points
Q <- 2 # outcome and treatment
S <- 2 # number of continuous covariates
P <- 2 # number of binary covariates
S_sum <- S*2+P # number of continuous covariates and interactions

Y <- rep(NA, I) # outcome
A <- rep(NA, I) # treatment 
U <- rep(NA, I) # latent variable (Gaussian)
E <- matrix(NA, nrow=I, ncol=Q) # outcome and treatment errors (non-Gaussian)
X <- matrix(NA, nrow=I, ncol=S_sum) # continuous covariates (include continuous interactions)
Z <- matrix(NA, nrow=I, ncol=P) # binary covariates

# independent errors
sigma2 <- rep(0.125, Q) # variance of errors
for (q in 1:Q){
  E[,q] <- rlaplace(I, location=0, scale=2*sqrt(sigma2[q])) 
}

# latent variables and covariates
U <- rnorm(I, 0, 1)
X[,1:S] <- rmvn_rcpp(I, rep(0,S), diag(0.5,S))
for (p in 1:P){ Z[,p] <- sample(c(0,1),I,replace=T,prob=c(0.4,0.6)) }

# spline coefficient
B0 <- 5 # dimension of the spline expansion
B <- B0 - 1 # sum-to-zero identifiability constraint
gamma <- rmvn_rcpp(1, rep(1,B), diag(1,B)) * sample(c(1,-1),B,replace=T)
alpha <- array(0, dim=c(Q, S_sum, B))
alpha[1,1,] <- rmvn_rcpp(1, rep(2,B), diag(1,B)) * sample(c(1,-1),B,replace=T)
alpha[1,3,] <- rmvn_rcpp(1, rep(2,B), diag(1,B)) * sample(c(1,-1),B,replace=T)
alpha[2,1,] <- rmvn_rcpp(1, rep(2,B), diag(1,B)) * sample(c(1,-1),B,replace=T)

# spline expansion: X[1:S]
X_bs0 <- array(NA, dim=c(I, S_sum, B0)) # b-spline expansion for continuous covariates 
for (s in 1:S){ X_bs0[,s,] <- bSpline(X[,s], df=B0) }
X_bs <- array(NA, dim=c(I, S_sum, B)) # b-spline expansion with constraint
Zx <- array(NA, dim=c(S_sum, B0, B)) # transformation matrix 
for (s in 1:S){
  Cx <- t(rep(1,I)) %*% X_bs0[,s,]
  qr_decomp_x <- qr(t(Cx))
  Zx[s,,] <- qr.Q(qr_decomp_x, complete=T)[,2:B0]
  X_bs[,s,] <- X_bs0[,s,] %*% Zx[s,,]
}

# treatment
A <- as.vector(X_bs[,1,]%*%alpha[2,1,] + 0.5*U + E[,2])
X[,(S+1):(S*2)] <- A*X[,1:S]; X[,(S*2+1):S_sum] <- A*Z

# spline expansion: A, X[S+1,S_sum]
A_bs0 <- bSpline(A, df=B0)  # b-spline expansion for treatment
A_bs <- array(NA, dim=c(I, B)) # b-spline expansion with constraint
Za <- array(NA, dim=c(B0, B)) # transformation matrix 
Ca <- t(rep(1,I)) %*% A_bs0
qr_decomp_a <- qr(t(Ca))
Za <- qr.Q(qr_decomp_a, complete=T)[,2:B0]
A_bs <- A_bs0 %*% Za

for (s in (S+1):S_sum){ X_bs0[,s,] <- bSpline(X[,s], df=B0) }
for (s in (S+1):S_sum){
  Cx <- t(rep(1,I)) %*% X_bs0[,s,]
  qr_decomp_x <- qr(t(Cx))
  Zx[s,,] <- qr.Q(qr_decomp_x, complete=T)[,2:B0]
  X_bs[,s,] <- X_bs0[,s,] %*% Zx[s,,]
}

# outcome
Y <- A_bs%*%as.vector(gamma) + X_bs[,S+1,]%*%alpha[1,3,] + 
  X_bs[,1,]%*%alpha[1,1,] + 0.5*U + E[,1]

# hyper-parameters
nu_0 <- 2.5e-4
a_nu <- 5; b_nu <- 50
a_rho <- 0.5; b_rho <- 0.5
a_sigma <- 1; b_sigma <- 1