### Model functions ###

# Geometric triggering function
geom_trig = function(t, beta){
  beta * (1-beta)^(t-1)
}

# Self-exciting component of conditional intensity function
calc_decay_vals = function(beta, decay_times, decay_times_counts, calc_grad=FALSE){
  
  # Calculations for decay function
  t = 1:max(unlist(decay_times))
  calc_decay = sapply(t, function(t) geom_trig(t, beta))
  decay_vals = lapply(decay_times, function(x) calc_decay[x])
  decay_vals_counts = Map('*',decay_vals,decay_times_counts)
  sum_decay = sapply(1:length(decay_vals_counts), function(x) sum(decay_vals_counts[[x]]))
  
  if (calc_grad==TRUE){
    # for gradients
    grad_vals_i = lapply(decay_times, function(i) (1-beta*i) * (1-beta)^(i-2))
    grad_vals_counts_i = Map('*',grad_vals_i,decay_times_counts)
    sum_grad_i = sapply(1:length(grad_vals_counts_i), function(x) sum(grad_vals_counts_i[[x]]))
    
    return(list(sum_decay, sum_grad_i))
  } else {
    return(list(sum_decay))
  }
  
}

# Conditional intensity function
lambda_cond_int = function(mu, alpha, beta, decay_times, decay_times_counts){
  
  sum_decay = calc_decay_vals(beta, decay_times, decay_times_counts)[[1]]
  lambda = mu + alpha * sum_decay
  
  return(lambda)
}

# Gradients of the log full conditional with respect to alpha and beta
grad_log_post = function(country, mu, alpha, beta, event_count, decay_times, decay_times_counts){
  # Set up
  grad = zeros(2, 1)
  
  # Calculations
  decay_vals = calc_decay_vals(beta, decay_times, decay_times_counts, calc_grad=TRUE)
  sum_decay = decay_vals[[1]]
  lambda = mu + alpha * sum_decay
  d_llambda_beta_num = alpha * decay_vals[[2]]
  
  # Gradient with respect to alpha
  grad[1] =  sum(as.vector(event_count[,2]) * (sum_decay / (lambda)) - sum_decay)
  # Gradient with respect to beta
  grad[2] = sum(as.vector(event_count[,2]) * (d_llambda_beta_num / (lambda)) - d_llambda_beta_num)
  
  
  return(grad)
}

# Log likelihood
log_likelihood = function(country, param, data){
  
  M = length(data$decay_times)
  
  f = function(m){
    lambda = lambda_cond_int(param$mu[m], param$alpha[m], param$beta[m], data$decay_times[[m]], data$decay_times_counts[[m]])
    return(lambda)
  }
  eval_lambda = lapply(1:M, f)
  log_lik= sum(vapply(1:M,function(m){sum(dpois(data$event_count[[m]][,2],lambda=eval_lambda[[m]],log=TRUE))},1))
  return(log_lik)
}
