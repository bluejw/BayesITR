
### q-functions
q_truth <- function(a, x, max_to_min){
  # calculate q-value
  y <- cos(a) + sin(x*a) 
  if (max_to_min){ return(-y) # from max to min
  }else { return(y) }
}

q_function_conservative <- function(a, x, phi){
  # transfrom data to spline basis
  a_bs <- predict(a_spline, a) %*% a_trans
  x_bs <- predict(x_spline, x) %*% x_trans
  xa_bs <- predict(xa_spline, x*a) %*% xa_trans
  # calculate q-value
  y <- cbind(rep(1,length(x)),a_bs,x_bs,xa_bs) %*% t(phi)
  y_min <- apply(y, 1, min) # lower bound
  # check whether out of the boundary
  out_index <- unique(c(which(a<a_lower|a>a_upper), 
                      which((x*a)<xa_lower|(x*a)>xa_upper)))
  y_min[out_index] <- 0
  return(y_min)
}

q_function_expected <- function(a, x, phi, sample_phi){
  # transfrom data to spline basis
  a_bs <- predict(a_spline, a) %*% a_trans
  x_bs <- predict(x_spline, x) %*% x_trans
  xa_bs <- predict(xa_spline, x*a) %*% xa_trans
  # calculate q-value
  if (sample_phi){ y <- as.vector(cbind(rep(1,length(x)),a_bs,x_bs,xa_bs) %*% phi[sample(1:num_samples,1),]) 
  }else{ y <- cbind(rep(1,length(x)),a_bs,x_bs,xa_bs) %*% t(phi); y <- rowMeans(y) }
  # check whether out of the boundary
  out_index <- unique(c(which(a<a_lower|a>a_upper),
                        which((x*a)<xa_lower|(x*a)>xa_upper)))
  y[out_index] <- 0
  return(y)
}

### sgd functions
policy_function <- function(x, I, theta, std_sgd){
  # transfrom data to spline basis
  x_bs <- predict(x_spline, x) %*% x_trans
  x_pred <- cbind(rep(1,length(x)), x_bs)
  # generate action
  a <- x_pred %*% theta + rnorm(I,0,std_sgd)
}

grad_expected_reward <- function(a, x, theta, y, std_sgd){
  # transfrom data to spline basis
  x_bs <- predict(x_spline, x) %*% x_trans
  x_pred <- cbind(rep(1,length(x)), x_bs)
  # calculate gradient of log policy
  grad_log_policy <- as.vector((a - x_pred %*% theta)/std_sgd^2) * x_pred
  # calculate the gradient of expected reward
  grad <- colSums(as.vector(y) * grad_log_policy)
  return(grad)
}