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


# Log likelihood
log_likelihood = function(country, mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2,event_count_up, decay_times_up, decay_times_counts_up, 
                          event_count_down, decay_times_down, decay_times_counts_down, i_star=NULL, cps=NULL){
  
  if ((!(country %in% c("Brazil", "India"))) & (is.null(i_star) || ((!is.null(i_star)) & i_star > cps[1]))){
    eval_lambda_up = lambda_cond_int(mu_1, alpha_1, beta_1, decay_times_up, decay_times_counts_up)
    eval_lambda_down = lambda_cond_int(mu_2, alpha_2, beta_2, decay_times_down, decay_times_counts_down)
    eval_lambda = c(eval_lambda_up, eval_lambda_down)
    event_counts = c(event_count_up[,2], event_count_down[,2])
  } else {
    eval_lambda_up = lambda_cond_int(mu_1, alpha_1, beta_1, decay_times_up, decay_times_counts_up)
    eval_lambda = eval_lambda_up
    event_counts = c(event_count_up[,2])
  }
  log_lik = sum(dpois(event_counts, lambda=eval_lambda, log=TRUE))
  return(log_lik)
}

# No summation
log_likelihood_i = function(i, country, mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2, event_count_up, decay_times_up, decay_times_counts_up, 
                            event_count_down, decay_times_down, decay_times_counts_down, cps=NULL){
  if ((!(country %in% c("Brazil", "India"))) & (i > cps[1])){
    eval_lambda_up = lambda_cond_int(mu_1, alpha_1, beta_1, decay_times_up, decay_times_counts_up)
    eval_lambda_down = lambda_cond_int(mu_2, alpha_2, beta_2, decay_times_down, decay_times_counts_down)
    eval_lambda = c(eval_lambda_up, eval_lambda_down)
    event_counts = c(event_count_up[,2], event_count_down[,2])
  } else {
    eval_lambda_up = lambda_cond_int(mu_1, alpha_1, beta_1, decay_times_up, decay_times_counts_up)
    eval_lambda = eval_lambda_up
    event_counts = c(event_count_up[,2])
  }
  log_lik = dpois(event_counts[i], lambda=eval_lambda[i], log=TRUE)
  return(log_lik)
}

# Gradients of the log full conditional with respect to alpha and beta - first phase
grad_log_post_1 = function(country, mu_1, alpha_1, beta_1, event_count_up, decay_times_up, decay_times_counts_up){

  # Set up
  grad = zeros(2, 1)
  
  # Calculations
  decay_vals_up = calc_decay_vals(beta_1, decay_times_up, decay_times_counts_up, calc_grad=TRUE)
  sum_decay_up = decay_vals_up[[1]]
  lambda_up = mu_1 + alpha_1 * sum_decay_up
  d_llambda_beta_num_up = alpha_1 * decay_vals_up[[2]]
  
  # Gradient with respect to alpha
  grad[1] =  sum(as.vector(event_count_up[,2]) * (sum_decay_up / (lambda_up)) - sum_decay_up)
  # Gradient with respect to beta
  grad[2] = sum(as.vector(event_count_up[,2]) * (d_llambda_beta_num_up / (lambda_up)) - d_llambda_beta_num_up)

  return(grad)
}

# Gradients of the log full conditional with respect to alpha and beta - second phase
grad_log_post_2 = function(mu_2, alpha_2, beta_2, event_count_down, decay_times_down, decay_times_counts_down){

  # Set up
  grad = zeros(2, 1)
 
  # Calculations
  decay_vals_down = calc_decay_vals(beta_2, decay_times_down, decay_times_counts_down, calc_grad=TRUE)
  sum_decay_down = decay_vals_down[[1]]
  lambda_down = mu_2 + alpha_2 * sum_decay_down
  d_llambda_beta_num_down = alpha_2 * decay_vals_down[[2]]

  # Gradient with respect to alpha
  grad[1] = sum(as.vector(event_count_down[,2]) * (sum_decay_down / (lambda_down)) - sum_decay_down)
  # Gradient with respect to beta
  grad[2] = sum(as.vector(event_count_down[,2]) * (d_llambda_beta_num_down / (lambda_down)) - d_llambda_beta_num_down)

  return(grad)
}
