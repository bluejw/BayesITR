
# select partial data
index_out <- which((X[,1] > 0 & A < 0) | (X[,1] < 0 & A > 0)) # data index outside the region
index_in <- which(!((1:I) %in% index_out)) # data index inside the region
index_select <- sample(index_out, floor(prob*length(index_out)), replace=F)
index <- sort(c(index_in,index_select))

I <- length(index)
Y <- Y[index]; A <- A[index]; X <- X[index,]; Z <- Z[index,]
U <- U[index]; E <- E[index,]

# spline expansion
B0 <- 5 # dimension of the spline expansion
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

# hyper-parameters
nu_0 <- 2.5e-4
a_nu <- 5; b_nu <- 50
a_rho <- 0.5; b_rho <- 0.5 
a_sigma <- 1; b_sigma <- 1
