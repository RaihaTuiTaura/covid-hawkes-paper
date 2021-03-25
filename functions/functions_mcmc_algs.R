
# MALA algorithm 
mala_cp = function(seed_id, country, input, par_algo, desc, prior_mu, m_start=3){

  # Start timer
  ptm=proc.time()
  
  # Set seed
  set.seed(NULL)
  set.seed(seed_id)
  
  param = input$param
  N = par_algo$N
  burnin = par_algo$burnin
  thin = par_algo$thin
  
  # Define training size
  if (!is.null(input$train)){
    i_star = input$train
  } else {
    i_star = NULL
  }
  
  # Prepare data
  if (is.null(i_star)){
    filter_data_cp_peak(country, par_algo$cps)
  } else if (!is.null(i_star)){
    filter_data_cp_peak(country, par_algo$cps, i_star)
  }
  
  M = length(data$decay_times)
  
  # Define random walk standard deviations
  sigma_mu = par_algo$s_mu
  
  # Initialise counters
  counter_mu = rep(0,M)
  counter_alpha_beta = rep(0,M)
  mu_out_of_bounds = rep(0,M)
  alpha_beta_out_of_bounds = rep(0,M)
  
  
  # Create shell for results
  parameters = list()
  parameters$mu = matrix(0., N, M)
  parameters$alpha = matrix(0., N, M)
  parameters$beta = matrix(0., N, M)
  

  # Calculate initial log likelihood and prior values
  log_lik = log_likelihood(country, param, data)

  if(prior_mu=="lognormal"){
    log_prior_mu = vapply(1:M, function(m) dlnorm(param$mu[m], mean_lmu, s_lmu, log = TRUE),1)
  } else if (prior_mu == "gamma"){
    log_prior_mu = vapply(1:M, function(m) dgamma(param$mu[m], shape=a_mu, scale=b_mu, log = TRUE),1) 
  } else if (prior_mu == "uniform"){
    log_prior_mu = rep(0,M)
  } 

  for (i in 1:N){
    
    ####### Update for mu #######
    
    for (m in m_start:M){
      
      sigma_mu_sample = c(0.1, 0.5, 1) * sigma_mu[m]
      S = sample(1:3,1)
      
      # Generate proposal for mu
      param_prop = param
      param_prop$mu[m] = rnorm(1, param$mu[m], sigma_mu_sample[S])
      
      if(param_prop$mu[m] > 0){
        
        # Evaluate posterior at proposed mu
        log_lik_prop = log_likelihood(country, param_prop, data)
        
        if(prior_mu=="lognormal"){
          log_prior_mu_prop = vapply(1:M, function(m) dlnorm(param_prop$mu[m], mean_lmu, s_lmu, log = TRUE),1)
        } else if (prior_mu == "gamma"){
          log_prior_mu_prop = vapply(1:M, function(m) dgamma(param_prop$mu[m], shape=a_mu, scale=b_mu, log = TRUE),1) 
        } else if (prior_mu == "uniform"){
          log_prior_mu_prop = rep(0,M)
        } 
  
        # Accept-reject step
        prob_accept_mu = sum(log_lik_prop - log_lik) + sum(log_prior_mu_prop - log_prior_mu)
        
        if (log(runif(1)) <= prob_accept_mu){
          param = param_prop
          log_lik = log_lik_prop
          log_prior_mu = log_prior_mu_prop
          counter_mu[m] = counter_mu[m] + 1
        }
      }else{
        mu_out_of_bounds[m] = mu_out_of_bounds[m] + 1
      }
    }
    
    ####### Joint update for alpha and beta with MALA step #######
    
    for (m in m_start:M){
      
      param_prop = param
      
      # Calculate gradients at current values of alpha and beta
      grad = grad_log_post(country, param$mu[m], param$alpha[m], param$beta[m], data$event_count[[m]], data$decay_times[[m]], data$decay_times_counts[[m]])
 

      alpha_beta_prop = as.vector(rmvnorm(1, mean = c(param$alpha[m],param$beta[m]) + 0.5 * par_algo$eps[m]^2 * par_algo$g_alpha_beta[[m]] %*% grad, sigma = par_algo$eps[m]^2 * par_algo$g_alpha_beta[[m]]))
      param_prop$alpha[m] = alpha_beta_prop[1]
      param_prop$beta[m] = alpha_beta_prop[2]
      
      if (param_prop$alpha[m] >= 0 & param_prop$beta[m] >= 0 & param_prop$beta[m] <= 1){
        
        # Calculate gradients at proposed alpha and beta
        grad_prop = grad_log_post(country, param_prop$mu[m], param_prop$alpha[m], param_prop$beta[m], data$event_count[[m]], data$decay_times[[m]], data$decay_times_counts[[m]])
        # Evaluate posterior at proposed alpha and beta
        log_lik_prop = log_likelihood(country, param_prop, data)


        # Evaluate candidate density at the current given proposed c(alpha,beta)
        q_ab_ab_c = dmvnorm(c(param$alpha[m],param$beta[m]), mean = alpha_beta_prop + 0.5 * par_algo$eps[m]^2 * par_algo$g_alpha_beta[[m]] %*% grad_prop, sigma = par_algo$eps[m]^2 * par_algo$g_alpha_beta[[m]], log=TRUE)
        # Evaluate candidate density at the proposed given current c(alpha,beta)
        q_ab_c_ab = dmvnorm(alpha_beta_prop, mean = c(param$alpha[m],param$beta[m]) + 0.5 * par_algo$eps[m]^2 * par_algo$g_alpha_beta[[m]] %*% grad, sigma = par_algo$eps[m]^2 * par_algo$g_alpha_beta[[m]], log=TRUE)
        
        # Accept-reject step
        prob_accept_alpha_beta = sum(log_lik_prop - log_lik) + sum(q_ab_ab_c - q_ab_c_ab)  
        
        if (log(runif(1)) <= prob_accept_alpha_beta){
          param = param_prop
          log_lik = log_lik_prop
          counter_alpha_beta[m] = counter_alpha_beta[m] + 1
        }
      }else{
        alpha_beta_out_of_bounds[m] = alpha_beta_out_of_bounds[m] + 1
      }
    }
    
    # Save results
    parameters$mu[i,] = param$mu[1:M]
    parameters$alpha[i,] = param$alpha[1:M]
    parameters$beta[i,] = param$beta[1:M]
  }
  
  if (is.null(input$test)){
    # Print acceptance rates
    for (m in 1:M){
      print(paste0("Acceptance rate for mu_", m))
      print(counter_mu[m]/N)
      print(paste0("Acceptance rate for alpha_beta_", m))
      print(counter_alpha_beta[m]/N)
    }
  }

  
  # Print posterior means
  print(colMeans(parameters$mu[seq(burnin,N,thin),]))
  print(colMeans(parameters$alpha[seq(burnin,N,thin),]))
  print(colMeans(parameters$beta[seq(burnin,N,thin),]))
  acc_rates = list()
  acc_rates$mu = counter_mu/N
  acc_rates$alpha_beta = counter_alpha_beta/N

  # End timer
  print(proc.time() - ptm)
  
  return(list(parameters, acc_rates))
}

