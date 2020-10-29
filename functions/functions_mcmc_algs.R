##### MCMC algorithms #####

# MALA algorithm with single change point
mala_cp = function(seed_id, country, N, burnin, thin, mu1_0, alpha1_0, beta1_0, mu2_0, alpha2_0, beta2_0, 
                   s_mu_1, g_alpha_beta_1, eps_1, s_mu_2, g_alpha_beta_2, eps_2, 
                   desc, cps, prior_mu, cv=NULL, train=NULL, test=NULL, interpol=FALSE){

  # Start timer
  ptm=proc.time()
  
  # Set seed
  set.seed(NULL)
  set.seed(seed_id)
  
  # Initialise parameters
  mu_1 = mu1_0 
  alpha_1 = alpha1_0
  beta_1 = beta1_0
  mu_2 = mu2_0 
  alpha_2 = alpha2_0
  beta_2 = beta2_0
  
  # Define random walk standard deviations
  sigma_mu_1 = s_mu_1
  sigma_mu_2 = s_mu_2

  # Initialise counters
  counter_mu_1 = 0
  counter_alpha_beta_1 = 0
  counter_mu_2 = 0
  counter_alpha_beta_2 = 0

  # Define training size
  if (!is.null(cv)){
    i_star = cv
  } else if (!is.null(train)){
    i_star = train
  } else {
    i_star = NULL
  }
  
  # Prepare data
  if (is.null(i_star)){
    filter_data_cp_peak(country, cps)
  } else if (!is.null(i_star)){
    filter_data_cp_peak(country, cps, i_star)
  }
  
  # Simulate missing data if doing interpolation
  if (interpol == TRUE){
    
    # Empty matrix to store simulated data
    if (country %in% c("Brazil", "India")){
      n_missing = ceil(0.1*length(decay_times_upward))
      missing_data = matrix(0., N+1, (n_missing))
      missing_ind = sample(1:length(decay_times_upward), n_missing, replace=FALSE)
    } else {
      n_missing = ceil(0.1*(length(decay_times_upward)+length(decay_times_downward)))
      missing_data = matrix(0., N+1, (n_missing))    
      missing_ind = sample(1:(length(decay_times_upward)+length(decay_times_downward)), n_missing, replace=FALSE)

    }
    
    # Set 10% of observations to missing
    missing_ind = missing_ind[order(missing_ind)]
    colnames(missing_data) = missing_ind
    event_count_sim = event_count
    event_count_sim[event_count_sim[,1] %in% missing_ind, 2] = NA
    event_count_upward[event_count_upward[,1] %in% missing_ind, 2] = NA
    if (!(country %in% c("Brazil", "India"))){
      event_count_downward[event_count_downward[,1] %in% missing_ind, 2] = NA
    }
    
    # Simulate missing data initially
    
    # Calculate decay function
    if (!(country %in% c("Brazil", "India"))){
      t = 1:max(unlist(get(paste0("decay_times_downward"))))
    } else {
      t = 1:max(unlist(get(paste0("decay_times_upward"))))
    }
    
    for (j in 1:length(missing_ind)){
      i = missing_ind[j]
      
      if ((!(country %in% c("Brazil", "India"))) & (i > cps[1])){
        calc_decay = sapply(t, function(t) as.numeric(geom_trig(t, beta_2)))  
        sim_decay_times = as.numeric(i- event_count_sim[event_count_sim[,1] < i & event_count_sim[,2] >0, 1])
        sim_decay_times_counts = as.numeric(event_count_sim[event_count_sim[,1] < i & event_count_sim[,2] >0, 2])
        sim_decay_vals = sapply(sim_decay_times, function(x) calc_decay[x])
        sim_decay_vals_counts = sim_decay_vals*sim_decay_times_counts
        lambda = as.numeric(mu_2) + as.numeric(alpha_2) * sum(sim_decay_vals_counts)
   
        
      } else {
        calc_decay = sapply(t, function(t) as.numeric(geom_trig(t, beta_1)))  
        sim_decay_times = as.numeric(i- event_count_sim[event_count_sim[,1] < i & event_count_sim[,2] >0, 1])
        sim_decay_times_counts = as.numeric(event_count_sim[event_count_sim[,1] < i & event_count_sim[,2] >0, 2])
        sim_decay_vals = sapply(sim_decay_times, function(x) calc_decay[x])
        sim_decay_vals_counts = sim_decay_vals*sim_decay_times_counts
        
        lambda = as.numeric(mu_1) + as.numeric(alpha_1) * sum(sim_decay_vals_counts)
      } 
      sim = rpois(1, lambda)
      event_count_sim[event_count_sim[,1]==i,2] = missing_data[1,j] = sim
    } 
    
    # Recalculate decay times
    event_count_upward = event_count_sim[(which(event_count_sim[,1] == 1):(which(event_count_sim[,1] == cps[1]))),]    
    decay_times_upward = lapply(1:cps[1], function(t) as.numeric(t- event_count_sim[event_count_sim[,1] <t & event_count_sim[,2] >0 ,1]))
    decay_times_counts_upward = lapply(1:cps[1], function(t) as.numeric(event_count_sim[event_count_sim[,1] <t & event_count_sim[,2] >0 ,2]))
    
    if (!(country %in% c("Brazil", "India"))){
      event_count_downward = event_count_sim[(which(event_count_sim[,1] == (cps[1]+1))):(which(event_count_sim[,1] == max(event_count_sim[,1]))),]
      decay_times_downward = lapply((cps[1]+1):max(event_count_downward[,1]), function(t) as.numeric(t- event_count_sim[event_count_sim[,1] <t & event_count_sim[,2] >0 ,1]))
      decay_times_counts_downward = lapply((cps[1]+1):max(event_count_downward[,1]), function(t) as.numeric(event_count_sim[event_count_sim[,1] <t & event_count_sim[,2] >0 ,2]))
    }
  }

  # Create shell for results
  parameters = matrix(0., N, 6)
  colnames(parameters) = c("mu_1", "alpha_1", "beta_1", "mu_2", "alpha_2", "beta_2")
  
  if (!is.null(cv)){
    if (country %in% c("Brazil", "India")){
      # Number of data points training on
      i_total = max(event_count_upward[,1])
    } else {
      # Number of data points training on
      i_total = max(event_count_downward[,1])
    }
    
    if (i_star > i_total){
      print('Error: i_star too high - validation set out of range')
    }
     
    # Shell for predictive log likelihoods
    log_lik_val = matrix(0., (N-burnin), (i_total-i_star))
  }
  
  # Calculate initial log likelihood and prior values
  if (is.null(i_star)){
    if (!(country %in% c("Brazil", "India"))){
      log_lik = log_likelihood(country, mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2, 
                               event_count_upward, decay_times_upward, decay_times_counts_upward, 
                               event_count_downward, decay_times_downward, decay_times_counts_downward)
    } else {
      log_lik = log_likelihood(country, mu_1, alpha_1, beta_1, event_count_up=event_count_upward, 
                               decay_times_up=decay_times_upward, decay_times_counts_up=decay_times_counts_upward)
    }
  } else if (i_star <= cps[1] || (country %in% c("Brazil", "India"))){
      log_lik = log_likelihood(country, mu_1, alpha_1, beta_1, event_count_up=event_count_upward_val, 
                               decay_times_up=decay_times_upward_val, decay_times_counts_up=decay_times_counts_upward_val,
                               i_star=i_star, cps=cps)
  } else if (i_star > cps[1]){
      log_lik = log_likelihood(country, mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2, 
                               event_count_upward_val, decay_times_upward_val, decay_times_counts_upward_val, 
                               event_count_downward_val, decay_times_downward_val, decay_times_counts_downward_val,
                               i_star=i_star, cps=cps)
  }  
  
  if(prior_mu=="lognormal"){
    log_prior_mu_1 = dlnorm(mu_1, mean_lmu, s_lmu, log = TRUE) 
    log_prior_mu_2 = dlnorm(mu_2, mean_lmu, s_lmu, log = TRUE) 
  } else if (prior_mu == "gamma"){
    log_prior_mu_1 = dgamma(mu_1, shape=a_mu, scale=b_mu, log = TRUE) 
    log_prior_mu_2 = dgamma(mu_2, shape=a_mu, scale=b_mu, log = TRUE) 
  } else if (prior_mu == "uniform"){
    log_prior_mu_1 = 0
    log_prior_mu_2 = 0
  } 

  for (i in 1:N){
    
    ####### Update for mu_1 #######
    sigma_mu_1_sample = c(0.1, 0.5, 1) * sigma_mu_1
    S = sample(1:3,1)

    # Generate proposal for mu_1
    mu_1_prop = rnorm(1, mu_1, sigma_mu_1_sample[S])

    if(mu_1_prop > 0){

      # Evaluate posterior at proposed mu_1
      if (is.null(i_star)){
        if (!(country %in% c("Brazil", "India"))){
          log_lik_prop = log_likelihood(country, mu_1_prop, alpha_1, beta_1, mu_2, alpha_2, beta_2, 
                                        event_count_upward, decay_times_upward, decay_times_counts_upward,
                                        event_count_downward, decay_times_downward, decay_times_counts_downward)
        } else {
          log_lik_prop = log_likelihood(country, mu_1_prop, alpha_1, beta_1, event_count_up=event_count_upward, 
                                        decay_times_up=decay_times_upward, decay_times_counts_up=decay_times_counts_upward)
        }
      } else if (i_star <= cps[1] || (country %in% c("Brazil", "India"))){
        log_lik_prop = log_likelihood(country, mu_1_prop, alpha_1, beta_1, event_count_up=event_count_upward_val, 
                                      decay_times_up=decay_times_upward_val, decay_times_counts_up=decay_times_counts_upward_val,
                                      i_star=i_star, cps=cps)
      } else if (i_star > cps[1]){
        log_lik_prop = log_likelihood(country, mu_1_prop, alpha_1, beta_1, mu_2, alpha_2, beta_2, 
                                      event_count_upward_val, decay_times_upward_val, decay_times_counts_upward_val,
                                      event_count_downward_val, decay_times_downward_val, decay_times_counts_downward_val,
                                      i_star=i_star, cps=cps)
      }

      if(prior_mu=="lognormal"){
        log_prior_mu_1_prop = dlnorm(mu_1_prop, mean_lmu, s_lmu, log = TRUE)
      } else if (prior_mu == "gamma"){
        log_prior_mu_1_prop = dgamma(mu_1_prop, shape=a_mu, scale=b_mu, log = TRUE)
      } else if (prior_mu == "uniform"){
        log_prior_mu_1_prop = 0
      }

      # Accept-reject step
      prob_accept_lmu_1 = sum(log_lik_prop - log_lik) + sum(log_prior_mu_1_prop - log_prior_mu_1)

      if (log(runif(1)) <= prob_accept_lmu_1){
        mu_1 = mu_1_prop
        log_lik = log_lik_prop
        log_prior_mu_1 = log_prior_mu_1_prop
        counter_mu_1 = counter_mu_1 + 1
      }
    }
    
    ####### Joint update for alpha_1 and beta_1 with MALA step #######
    
    # Calculate gradients at current values of alpha_1 and beta_1
    if (is.null(i_star)){
        grad_1 = grad_log_post_1(country, mu_1, alpha_1, beta_1, event_count_upward, decay_times_upward, decay_times_counts_upward)
  
    } else {
        grad_1 = grad_log_post_1(country, mu_1, alpha_1, beta_1, event_count_upward_val, decay_times_upward_val, decay_times_counts_upward_val)
    } 

  alpha_beta_1_prop = as.vector(rmvnorm(1, mean = c(alpha_1,beta_1) + 0.5 * eps_1^2 * g_alpha_beta_1 %*% grad_1, sigma = eps_1^2 * g_alpha_beta_1))
  alpha_1_prop = alpha_beta_1_prop[1]
  beta_1_prop = alpha_beta_1_prop[2]
  
     if (alpha_1_prop >= 0 & beta_1_prop >= 0 & beta_1_prop <= 1){

       if (is.null(i_star)){
  
         if (!(country %in% c("Brazil", "India"))){
           # Calculate gradients at proposed alpha_1 and beta_1
           grad_1_prop = grad_log_post_1(country, mu_1, alpha_1_prop, beta_1_prop, event_count_upward, decay_times_upward, decay_times_counts_upward)
           # Evaluate posterior at proposed alpha_1 and beta_1
           log_lik_prop = log_likelihood(country, mu_1, alpha_1_prop, beta_1_prop, mu_2, alpha_2, beta_2, 
                                        event_count_upward, decay_times_upward, decay_times_counts_upward,
                                        event_count_downward, decay_times_downward, decay_times_counts_downward)
         } else {
           # Calculate gradients at proposed alpha_1 and beta_1
           grad_1_prop = grad_log_post_1(country, mu_1, alpha_1_prop, beta_1_prop, event_count_upward, decay_times_upward, decay_times_counts_upward)
           # Evaluate posterior at proposed alpha_1 and beta_1
           log_lik_prop = log_likelihood(country, mu_1, alpha_1_prop, beta_1_prop, event_count_up=event_count_upward, 
                                         decay_times_up=decay_times_upward, decay_times_counts_up=decay_times_counts_upward)
         }
         
      } else if (i_star <= cps[1] || (country %in% c("Brazil", "India"))){
        # Calculate gradients at proposed alpha_1 and beta_1
        grad_1_prop = grad_log_post_1(country, mu_1, alpha_1_prop, beta_1_prop, event_count_upward_val, decay_times_upward_val, decay_times_counts_upward_val)
        # Evaluate posterior at proposed alpha_1 and beta_1
        log_lik_prop = log_likelihood(country, mu_1, alpha_1_prop, beta_1_prop, event_count_up=event_count_upward_val, 
                                      decay_times_up=decay_times_upward_val, decay_times_counts_up=decay_times_counts_upward_val,
                                      i_star=i_star, cps=cps)
        
      } else if (i_star > cps[1]){
        # Calculate gradients at proposed alpha_1 and beta_1
        grad_1_prop = grad_log_post_1(country, mu_1, alpha_1_prop, beta_1_prop, event_count_upward_val, decay_times_upward_val, decay_times_counts_upward_val)
        # Evaluate posterior at proposed alpha_1 and beta_1
        log_lik_prop = log_likelihood(country, mu_1, alpha_1_prop, beta_1_prop, mu_2, alpha_2, beta_2,  
                                      event_count_upward_val, decay_times_upward_val, decay_times_counts_upward_val,
                                      event_count_downward_val, decay_times_downward_val, decay_times_counts_downward_val,
                                      i_star=i_star, cps=cps)
      } 
      
      # Evaluate candidate density at the current given proposed c(alpha_1,beta_1)
      q_ab_ab_c_1 = dmvnorm(c(alpha_1,beta_1), mean = alpha_beta_1_prop + 0.5 * eps_1^2 * g_alpha_beta_1 %*% grad_1_prop, sigma = eps_1^2 * g_alpha_beta_1, log=TRUE)
      # Evaluate candidate density at the proposed given current c(alpha_1,beta_1)
      q_ab_c_ab_1 = dmvnorm(alpha_beta_1_prop, mean = c(alpha_1,beta_1) + 0.5 * eps_1^2 * g_alpha_beta_1 %*% grad_1, sigma = eps_1^2 * g_alpha_beta_1, log=TRUE)

      # Accept-reject step
      prob_accept_alpha_beta_1 = sum(log_lik_prop - log_lik) + sum(q_ab_ab_c_1 - q_ab_c_ab_1)  

      if (log(runif(1)) <= prob_accept_alpha_beta_1){
        alpha_1 = alpha_1_prop
        beta_1 = beta_1_prop
        log_lik = log_lik_prop
        counter_alpha_beta_1 = counter_alpha_beta_1 + 1
      }
    }
    
    ####### Update for mu_2 #######

    if (!(country %in% c("Brazil", "India"))){
      if (is.null(i_star) || ((!is.null(i_star)) & (i_star > cps[1]))){

        sigma_mu_2_sample = c(0.1, 0.5, 1) * sigma_mu_2
        S = sample(1:3,1)

        # Generate proposal for mu_2
        mu_2_prop = rnorm(1, mu_2, sigma_mu_2_sample[S])

        if(mu_2_prop > 0){

          # Evaluate posterior at proposed mu_2
          if (is.null(i_star) ){
            log_lik_prop = log_likelihood(country, mu_1, alpha_1, beta_1, mu_2_prop, alpha_2, beta_2, 
                                          event_count_upward, decay_times_upward, decay_times_counts_upward,
                                          event_count_downward, decay_times_downward, decay_times_counts_downward)
          } else {
            log_lik_prop = log_likelihood(country, mu_1, alpha_1, beta_1, mu_2_prop, alpha_2, beta_2, 
                                          event_count_upward_val, decay_times_upward_val, decay_times_counts_upward_val,
                                          event_count_downward_val, decay_times_downward_val, decay_times_counts_downward_val,
                                          i_star=i_star, cps=cps)
          }

          if(prior_mu=="lognormal"){
            log_prior_mu_2_prop = dlnorm(mu_2_prop, mean_lmu, s_lmu, log = TRUE)
          } else if (prior_mu == "gamma"){
            log_prior_mu_2_prop = dgamma(mu_2_prop, shape=a_mu, scale=b_mu, log = TRUE)
          } else if (prior_mu == "uniform"){
            log_prior_mu_2_prop = 0
          }

          # Accept-reject step
          prob_accept_lmu_2 = sum(log_lik_prop - log_lik) + sum(log_prior_mu_2_prop - log_prior_mu_2)

          if (log(runif(1)) <= prob_accept_lmu_2){
            mu_2 = mu_2_prop
            log_lik = log_lik_prop
            log_prior_mu_2 = log_prior_mu_2_prop
            counter_mu_2 = counter_mu_2 + 1
          }
        }
      }
    }

    ####### Joint update for alpha_2 and beta_2 with MALA step #######

    if (!(country %in% c("Brazil", "India"))){

      if (is.null(i_star) || ((!is.null(i_star)) & i_star > cps[1])){

        # Calculate gradients at current values of alpha_2 and beta_2
        if (is.null(i_star)){
          grad_2 = grad_log_post_2(mu_2, alpha_2, beta_2, event_count_downward, decay_times_downward, decay_times_counts_downward)
        } else {
          grad_2 = grad_log_post_2(mu_2, alpha_2, beta_2,event_count_downward_val, decay_times_downward_val, decay_times_counts_downward_val)
        }

        alpha_beta_2_prop = as.vector(rmvnorm(1, mean = c(alpha_2,beta_2) + 0.5 * eps_2^2 * g_alpha_beta_2 %*% grad_2, sigma = eps_2^2 * g_alpha_beta_2))
        alpha_2_prop = alpha_beta_2_prop[1]
        beta_2_prop = alpha_beta_2_prop[2]
        
      if (alpha_2_prop >= 0 & beta_2_prop >= 0 & beta_2_prop <= 1){

          if (is.null(i_star)){

            # Calculate gradients at proposed alpha_2 and beta_2
            grad_2_prop = grad_log_post_2(mu_2, alpha_2_prop, beta_2_prop, event_count_downward, decay_times_downward, decay_times_counts_downward)

            # Evaluate posterior at proposed alpha_2 and beta_2
            log_lik_prop = log_likelihood(country, mu_1, alpha_1, beta_1, mu_2, alpha_2_prop, beta_2_prop, 
                                          event_count_upward, decay_times_upward, decay_times_counts_upward,
                                          event_count_downward, decay_times_downward, decay_times_counts_downward)
          } else {
            # Calculate gradients at proposed alpha_2 and beta_2
            grad_2_prop = grad_log_post_2(mu_2, alpha_2_prop, beta_2_prop, 
                                          event_count_downward_val, decay_times_downward_val, decay_times_counts_downward_val)

            # Evaluate posterior at proposed alpha_2 and beta_2
            log_lik_prop = log_likelihood(country, mu_1, alpha_1, beta_1, mu_2, alpha_2_prop, beta_2_prop,  
                                          event_count_upward_val, decay_times_upward_val, decay_times_counts_upward_val,
                                          event_count_downward_val, decay_times_downward_val, decay_times_counts_downward_val,
                                          i_star=i_star, cps=cps)
          }

          # Evaluate candidate density at the current given proposed c(alpha_2,beta_2)
          q_ab_ab_c_2 = dmvnorm(c(alpha_2,beta_2), mean = alpha_beta_2_prop + 0.5 * eps_2^2 * g_alpha_beta_2 %*% grad_2_prop, sigma = eps_2^2 * g_alpha_beta_2, log=TRUE)
          # Evaluate candidate density at the proposed given current c(alpha_2,beta_2)
          q_ab_c_ab_2 = dmvnorm(alpha_beta_2_prop, mean = c(alpha_2,beta_2) + 0.5 * eps_2^2 * g_alpha_beta_2 %*% grad_2, sigma = eps_2^2 * g_alpha_beta_2, log=TRUE)

          # Accept-reject step
          prob_accept_alpha_beta_2 = sum(log_lik_prop - log_lik) + sum(q_ab_ab_c_2 - q_ab_c_ab_2)

          if (log(runif(1)) <= prob_accept_alpha_beta_2){
            alpha_2 = alpha_2_prop
            beta_2 = beta_2_prop
            log_lik = log_lik_prop
            counter_alpha_beta_2 = counter_alpha_beta_2 + 1
          }
        }
      }
    }
  
    # Save results
    parameters[i,] = c(mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2)
    
    if (!is.null(cv)){
      # Calculate predictive log_likelihood
      if (i>burnin){
        if (country %in% c("Brazil", "India")){
          log_lik_val[(i-burnin),]=sapply((i_star+1):(i_total),function(x) {log_likelihood_i(x, country, mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2,
                                                                                             event_count_up=event_count_upward, decay_times_up=decay_times_upward, decay_times_counts_up=decay_times_counts_upward,
                                                                                             cps=cps)})
        } else if (i_star < cps[1]){
          log_lik_val[(i-burnin),1:(cps[1]-i_star)]=sapply((i_star+1):cps[1],function(x) {log_likelihood_i(x, country, mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2, 
                                                                                                           event_count_up=event_count_upward, decay_times_up=decay_times_upward, decay_times_counts_up=decay_times_counts_upward,
                                                                                          cps=cps)})
          log_lik_val[(i-burnin),(cps[1]-i_star+1):(i_total-i_star)]=sapply((cps[1]+1):(i_total),function(x) {log_likelihood_i(x, country, mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2, 
                                                                                                                               event_count_upward, decay_times_upward, decay_times_counts_upward, 
                                                                                                                               event_count_downward, decay_times_downward, decay_times_counts_downward, 
                                                                                                                               cps)})
        } else if (i_star >= cps[1]){
          log_lik_val[(i-burnin),]=sapply((i_star+1):(i_total),function(x) {log_likelihood_i(x, country, mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2,  
                                                                                             event_count_upward, decay_times_upward, decay_times_counts_upward, 
                                                                                             event_count_downward, decay_times_downward, decay_times_counts_downward, 
                                                                                             cps)})
        } 
      }
    }
    
    # Simulate missing data 
  
    if (interpol==TRUE){
      
      # Simulate missing data using updated parameters
      event_count_sim[event_count_sim[,1] %in% missing_ind, 2] = NA
      event_count_upward[event_count_upward[,1] %in% missing_ind, 2] = NA
      if (!(country %in% c("Brazil", "India"))){
        event_count_downward[event_count_downward[,1] %in% missing_ind, 2] = NA
      }
      # Calculate decay function
      if (!(country %in% c("Brazil", "India"))){
        t = 1:max(unlist(get(paste0("decay_times_downward"))))
      } else {
        t = 1:max(unlist(get(paste0("decay_times_upward"))))
      }
      
      for (j in 1:length(missing_ind)){
        m = missing_ind[j]
        
        if ((!(country %in% c("Brazil", "India"))) & (m > cps[1])){
          calc_decay = sapply(t, function(t) as.numeric(geom_trig(t, beta_2)))  
          sim_decay_times = as.numeric(m- event_count_sim[event_count_sim[,1] < m  & event_count_sim[,2] >0, 1])
          sim_decay_times_counts = as.numeric(event_count_sim[event_count_sim[,1] < m  & event_count_sim[,2] >0, 2])
          sim_decay_vals = sapply(sim_decay_times, function(x) calc_decay[x])
          sim_decay_vals_counts = sim_decay_vals*sim_decay_times_counts
          lambda = as.numeric(mu_2) + as.numeric(alpha_2) * sum(sim_decay_vals_counts)
 
      } else {
          calc_decay = sapply(t, function(t) as.numeric(geom_trig(t, beta_1)))  
          sim_decay_times = as.numeric(m- event_count_sim[event_count_sim[,1] < m & event_count_sim[,2] >0, 1])
          sim_decay_times_counts = as.numeric(event_count_sim[event_count_sim[,1] < m & event_count_sim[,2] >0, 2])
          sim_decay_vals = sapply(sim_decay_times, function(x) calc_decay[x])
          sim_decay_vals_counts = sim_decay_vals*sim_decay_times_counts
          
          lambda = as.numeric(mu_1) + as.numeric(alpha_1) * sum(sim_decay_vals_counts)
        } 
        sim = rpois(1, lambda)
        event_count_sim[event_count_sim[,1]==m,2] = missing_data[(i+1),j] = sim
      } 
     
      # Recalculate decay times
      event_count_upward = event_count_sim[(which(event_count_sim[,1] == 1):(which(event_count_sim[,1] == cps[1]))),]    
      decay_times_upward = lapply(1:cps[1], function(t) as.numeric(t- event_count_sim[event_count_sim[,1] <t & event_count_sim[,2] >0 ,1]))
      decay_times_counts_upward = lapply(1:cps[1], function(t) as.numeric(event_count_sim[event_count_sim[,1] <t & event_count_sim[,2] >0 ,2]))
      
      if (!(country %in% c("Brazil", "India"))){
        event_count_downward = event_count_sim[(which(event_count_sim[,1] == (cps[1]+1))):(which(event_count_sim[,1] == max(event_count_sim[,1]))),]
        decay_times_downward = lapply((cps[1]+1):max(event_count_downward[,1]), function(t) as.numeric(t- event_count_sim[event_count_sim[,1] <t & event_count_sim[,2] >0 ,1]))
        decay_times_counts_downward = lapply((cps[1]+1):max(event_count_downward[,1]), function(t) as.numeric(event_count_sim[event_count_sim[,1] <t & event_count_sim[,2] >0 ,2]))
      }

      # Update log likelihood
      if (!(country %in% c("Brazil", "India"))){
        log_lik = log_likelihood(country, mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2, 
                                 event_count_upward, decay_times_upward, decay_times_counts_upward, 
                                 event_count_downward, decay_times_downward, decay_times_counts_downward)
      } else {
        log_lik = log_likelihood(country, mu_1, alpha_1, beta_1, event_count_up=event_count_upward, 
                                 decay_times_up=decay_times_upward, decay_times_counts_up=decay_times_counts_upward)
      }
    }
  }

  if (!is.null(test)){
    oos_prediction_intervals = post_pred_oos(country, burnin, N, parameters, cps, 1000, train, test)
  }
  
  # Calculate acceptance rates
  acc_rates = c(as.numeric(counter_mu_1/N),as.numeric(counter_alpha_beta_1/N),
                as.numeric(counter_mu_2/N),as.numeric(counter_alpha_beta_2/N))

  # End timer
  print(proc.time() - ptm)
  
  if (is.null(cv) & interpol==FALSE & is.null(test)){
    return(list(parameters, acc_rates))
  } else if (!is.null(cv)){
    return(list(parameters, acc_rates, log_lik_val))
  } else if (interpol==TRUE){
    return(list(parameters, acc_rates, missing_data))
  } else if (!is.null(test)){
    return(list(parameters, acc_rates, oos_prediction_intervals))
  }
}

# Out-of-sample prediction
sim_data_pred_oos = function(country, par, cps, train, test){

  # Empty data from starting point
  event_count_pred = event_count
  event_count_pred[event_count_pred[,1]>train,2] = NA
  
  # Calculate decay function
  if (!(country %in% c("Brazil", "India"))){
    t = 1:max(unlist(get(paste0("decay_times_downward"))))
  } else {
    t = 1:max(unlist(get(paste0("decay_times_upward"))))
  }
  
  for (j in 1:test){
    i= train+j

    if ((!(country %in% c("Brazil", "India"))) & (train > cps[1])){
      calc_decay = sapply(t, function(t) as.numeric(geom_trig(t, par[6])))  
      pred_decay_times = as.numeric(i- event_count_pred[event_count_pred[,1] < i & event_count_pred[,1] > cps[1] & event_count_pred[,2] >0, 1])
      pred_decay_times_counts = as.numeric(event_count_pred[event_count_pred[,1] < i & event_count_pred[,1] > cps[1] & event_count_pred[,2] >0, 2])
      pred_decay_vals = sapply(pred_decay_times, function(x) calc_decay[x])
      pred_decay_vals_counts = pred_decay_vals*pred_decay_times_counts
 
      lambda = as.numeric(par[4]) + as.numeric(par[5]) * sum(pred_decay_vals_counts) 
      
    } else {
      calc_decay = sapply(t, function(t) as.numeric(geom_trig(t, par[3])))  
      pred_decay_times = as.numeric(i- event_count_pred[event_count_pred[,1] < i & event_count_pred[,2] >0, 1])
      pred_decay_times_counts = as.numeric(event_count_pred[event_count_pred[,1] < i & event_count_pred[,2] >0, 2])
      pred_decay_vals = sapply(pred_decay_times, function(x) calc_decay[x])
      pred_decay_vals_counts = pred_decay_vals*pred_decay_times_counts
         
      lambda = as.numeric(par[1]) + as.numeric(par[2]) * sum(pred_decay_vals_counts)
    } 
    event_count_pred[event_count_pred[,1]==i,2] = rpois(1, lambda)
  } 

  
  # Calculate all decay times
  event_count_pred = as.matrix(event_count_pred)
  event_count_pred[event_count_pred[,1]<=train,2] = NA
  return(event_count_pred[,2])
}

post_pred_oos = function(country, burnin, N, results, cps, n_samples, train, test){
  sample_ind = sample(burnin:N, n_samples, replace=TRUE)
  sim_pred = t(sapply(sample_ind, function(i) {sim_data_pred_oos(country, results[i,], cps, train, test)}))
  quantiles = t(apply(sim_pred, 2, function(x) quantile(x, c(0.025,0.5,0.975), na.rm=TRUE)))
  colnames(quantiles) = c(paste0("pi_0.025_", train), paste0("pi_0.5_", train), paste0("pi_0.975_", train))
  return(quantiles)
}
