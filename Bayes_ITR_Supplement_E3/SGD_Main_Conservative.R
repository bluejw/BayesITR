####################################################### SGD ########################################################

sgd <- NULL
sgd$theta <- matrix(0, nrow=Niter, ncol=B0)
sgd$grad <- matrix(0, nrow=Niter, ncol=B0)
sgd$y <- matrix(0, nrow=Niter, ncol=I)
sgd$a <- matrix(0, nrow=Niter, ncol=I)
sgd$r <- rep(0, Niter)

### Start of the algorithm
start.time = proc.time()

for (iter in 2:Niter){
  
  print(paste(seed_index,iter,sep=","))
  
  # generate state and action 
  sgd$a[iter,] <- policy_function(x, I, sgd$theta[iter-1,], std_sgd)
  sgd$y[iter,] <- q_function_conservative(sgd$a[iter,], x, phi_cr)
  sgd$r[iter] <- sum(sgd$y[iter,])
  
  # calculate the gradient of expected reward
  sgd$grad[iter,] <- grad_expected_reward(sgd$a[iter,], x, sgd$theta[iter-1,], sgd$y[iter,], std_sgd)
  
  # update the policy parameter
  sgd$theta[iter,] <- sgd$theta[iter-1,] + step_size*sgd$grad[iter,]
}

duration = proc.time()-start.time
print(duration) # print the running time 
### End of the algorithm
