
### Code to perform PSIS-LFO-CV (adapted from https://github.com/paul-buerkner/LFO-CV-paper) ###

#some helper functions we'll use throughout
SW <- suppressWarnings
SM <- suppressMessages

# more stable than log(sum(exp(x))) 
log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

# compute log of raw importance ratios
# sums over observations *not* over posterior samples
sum_log_ratios <- function(ll, ids = NULL, mode) {
  if (!is.null(ids)) {
    ll <- ll[, ids, drop = FALSE]
  }
  out <- rowSums(ll)
  if (mode == "backward") {
    out <- -out
  }
  out
}


scale_unit_interval <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

summarize_elpds <- function(elpds) {
  elpds <- na.omit(elpds)
  c(Estimate = sum(elpds), SE = sqrt(length(elpds) * var(elpds)))
}

plot_ks <- function(ks, ids, k_thres = 0.7) {
  dat_ks <- data.frame(ks = ks, ids = ids)
  ggplot(dat_ks, aes(x = ids, y = ks)) + 
    geom_point(aes(color = ks > k_thres), shape = 3, show.legend = FALSE) + 
    geom_hline(yintercept = k_thres, linetype = 2, color = "red2") + 
    scale_color_manual(values = c("cornflowerblue", "darkblue")) + 
    labs(x = "Data point", y = "Pareto k") + 
    ylim(-0.5, 1.5)
}


# Perform approximate leave-future-out cross-validation
# @param M steps to predict into the future
# @param L minimal number of observations required for fitting the model
# @return A vector of pointwise ELPD values

approx_lfo <- function(fit, data, M, L, k_thres = 0.7, 
                       mode = c("forward", "backward", "combined"),
                       criterion = c("elpd", "rmse"), factorize = FALSE,
                       seed_id, country, N_iters, burnin, thin, mu1_0, alpha1_0, beta1_0, 
                       mu2_0, alpha2_0, beta2_0, s_mu_1, g_alpha_beta_1, eps_1,
                       s_mu_2, g_alpha_beta_2, eps_2, cps, prior_mu, ...) {
  require(loo)
  mode <- match.arg(mode)
  criterion <- match.arg(criterion)
  data = as.data.frame(data)
  N <- nrow(data)
  
  # compute approximate LFO likelihoods
  out <- ks <- reffs <- rep(NA, N)
  refits <- numeric(0)
  if (mode %in% c("forward", "combined")) {
    # move forward in time
    # need to refit the model with the first L observations
    reffs[L] <- 1
    
    tmp <- lfo_step_pred(
      i = L, fit_star = fit, data = data, M = M, 
      criterion = criterion, factorize = factorize,
      seed_id=seed_id, country=country, N_iters=N_iters, burnin=burnin, thin=thin, mu1_0=mu1_0, alpha1_0=alpha1_0, beta1_0=beta1_0, 
      mu2_0=mu2_0, alpha2_0=alpha2_0, beta2_0=beta2_0, s_mu_1=s_mu_1, g_alpha_beta_1=g_alpha_beta_1, eps_1=eps_1,
      s_mu_2=s_mu_2, g_alpha_beta_2=g_alpha_beta_2, eps_2=eps_2, 
      cps=cps, prior_mu=prior_mu
    )
    fit_star <- tmp$fit_star
    i_star <- tmp$i_star
    ks[L] <- tmp$k
    refits[L] <- tmp$refit
    out[L] <- tmp$crit 
    # start from L + 1 as we already handled L above
    ids <- (L + 1):(N - M)
  } else {
    # move backward in time
    fit_star <- fit
    i_star <- N
    ids <- (N - M):L
  }
  
  for (i in ids) {
    if (mode == "combined") {
      i_star_old <- i_star
    }
    psis_i <- lfo_step_psis(
      i = i, fit_star = fit_star, i_star = i_star, 
      data = data, mode = mode, factorize = factorize
    )
    reffs[i] <- attr(psis_i, "r_eff_lr")
    
    tmp <- lfo_step_pred(
      i = i, fit_star = fit_star, i_star = i_star, 
      data = data, M = M, psis = psis_i, k_thres = k_thres,
      criterion = criterion, factorize = factorize,
      seed_id=seed_id, country=country, N_iters=N_iters, burnin=burnin, thin=thin, mu1_0=mu1_0, alpha1_0=alpha1_0, beta1_0=beta1_0, 
      mu2_0=mu2_0, alpha2_0=alpha2_0, beta2_0=beta2_0, s_mu_1=s_mu_1, g_alpha_beta_1=g_alpha_beta_1, eps_1=eps_1,
      s_mu_2=s_mu_2, g_alpha_beta_2=g_alpha_beta_2, eps_2=eps_2, 
      cps=cps, prior_mu=prior_mu
    )
    fit_star <- tmp$fit_star
    i_star <- tmp$i_star
    ks[i] <- tmp$k
    refits[i] <- tmp$refit
    out[i] <- tmp$crit
    
    if (i_star == i && mode == "combined") {
      # model was just refitted now also move backwards in time
      ids_bw <- (i_star - 1):(i_star_old + 1) 
      for (j in ids_bw) {
        psis_j <- lfo_step_psis(
          i = j, fit_star = fit_star, i_star = i_star, 
          data = data, mode = "backward", factorize = factorize
        )
        k_j <- pareto_k_values(psis_j)
        if (k_j > k_thres) {
          # stop backward in forward if a refit is necessary
          break
        } else {
          tmp <- lfo_step_pred(
            i = j, fit_star = fit_star, i_star = i_star, 
            data = data, M = M, psis = psis_j, k_thres = k_thres,
            criterion = criterion, factorize = factorize
          )
          out[j] <- (out[j] + tmp$crit) / 2 
        }
      }
    }
  }
  attr(out, "ks") <- ks
  attr(out, "reffs") <- reffs
  attr(out, "refits") <- refits
  out
}

lfo_step_psis <- function(i, fit_star, data, i_star, mode,
                          factorize = FALSE) {
  # compute psis over observations in J_i 
  if (mode %in% c("forward", "combined")) {
    J_i <- (i_star + 1):i 
  } else {
    J_i <- (i + 1):i_star
  }
  if (factorize) {
    data_psis <- data[J_i, , drop = FALSE]
    J_i <- seq_len(NROW(data_psis))
  } else {
    data_psis <- data[1:max(J_i), , drop = FALSE]
  }
  
  ll_psis = as.matrix(fit_star[[3]][,J_i])
  logratio <- sum_log_ratios(ll_psis, ids = J_i, mode = mode)
  
  # compute relative efficiency of logratios
  chains=1
  chain_id <- rep(seq_len(chains), each = NROW(logratio) / chains)
  r_eff_lr <- loo::relative_eff(logratio, chain_id = chain_id)
  
  r_eff <- loo::relative_eff(1 / exp(logratio), chain_id = chain_id)
  psis <- SW(psis(logratio, r_eff = r_eff))
  attr(psis, "r_eff_lr") <- r_eff_lr
  psis
}


lfo_step_pred <- function(i, fit_star, data, M, criterion, factorize = FALSE,
                          k_thres = Inf, i_star = i, psis = NULL, 
                          seed_id, country, N_iters, burnin, thin, mu1_0, alpha1_0, beta1_0, 
                          mu2_0, alpha2_0, beta2_0, s_mu_1, g_alpha_beta_1, eps_1,
                          s_mu_2, g_alpha_beta_2, eps_2, cps, prior_mu, ...){
  
  # observations to predict
  oos <- (i + 1):(i + M)
  data_pred <- data[1:(i + M), , drop = FALSE]
  
  if (!is.null(psis)) {
    k <- pareto_k_values(psis) 
  } else {
    k <- 0
  }
  if (is.null(psis) || k > k_thres) {
    # refit the model based on the first i observations
    print(i)
    refit <- TRUE
    i_star <- i
    fit_star <- mala_cp(seed_id, country, N_iters, burnin, thin, mu1_0, alpha1_0, beta1_0, mu2_0, alpha2_0, beta2_0, 
                              s_mu_1, g_alpha_beta_1, eps_1, s_mu_2, g_alpha_beta_2, eps_2, 
                              "_validation_", cps, prior_mu, cv=i_star)

    # perform exact LFO-CV for the ith observation
    crit <- comp_criterion(
      fit_star, data = data_pred, oos = oos, 
      criterion = criterion, factorize = factorize,
      i_star =i_star
    )
  } else {
    # PSIS approximate LFO-CV is possible
    refit <- FALSE
    crit <- comp_criterion(
      fit_star, data = data_pred, oos = oos,
      criterion = criterion, psis = psis,
      factorize = factorize,
      i_star =i_star
    )
  }
  list(
    crit = crit,
    i_star = i_star,
    fit_star = fit_star,
    k = k, 
    refit = refit
  )
}

comp_criterion <- function(fit, data, oos, criterion = c("elpd", "rmse"),
                           factorize = FALSE, i_star, psis = NULL, ...) {

  data <- as.data.frame(data)
  stopifnot(all(oos %in% seq_len(NROW(data))))
  criterion <- match.arg(criterion)
  stopifnot(is.null(psis) || loo:::is.psis(psis))

  data <- as.data.frame(data[oos,])
  oos <- oos - i_star

  if (criterion == "elpd") {
    ll <- fit[[3]][,oos]
    oos_ind = seq_len(NROW(data))
    sum_ll <- rowSums(ll[, oos_ind, drop = FALSE])
    if (is.null(psis)) {
      out <- log_mean_exp(sum_ll)
    } else {
      lw <- weights(psis, normalize = TRUE)[, 1]
      out <- log_sum_exp(lw + sum_ll)
    }
  } 

  out
}

compute_lfo <- function(fit, type = c("exact", "approx"), file = NULL, 
                        seed_id, country, N_iters, burnin, thin, mu1_0, alpha1_0, beta1_0, 
                        mu2_0, alpha2_0, beta2_0, s_mu_1, g_alpha_beta_1, eps_1,
                        s_mu_2, g_alpha_beta_2, eps_2, cps, prior_mu,...){

  type <- match.arg(type)
  if (!is.null(file) && file.exists(file)) {
    return(readRDS(file))
  }
  filter_data_cp_peak(country, cps)
  if (!(country %in% c("Brazil","India"))){
    data = rbind(event_count_upward, event_count_downward)
  } else {
    data = event_count_upward
  }
  lfo_fun <- get(paste0(type, "_lfo"), mode = "function")
  
  out <- lfo_fun(fit,data,
                 seed_id=seed_id, country=country, N_iters=N_iters, burnin=burnin, thin=thin, mu1_0=mu1_0, alpha1_0=alpha1_0, beta1_0=beta1_0, 
                 mu2_0=mu2_0, alpha2_0=alpha2_0, beta2_0=beta2_0, s_mu_1=s_mu_1, g_alpha_beta_1=g_alpha_beta_1, eps_1=eps_1,
                 s_mu_2=s_mu_2, g_alpha_beta_2=g_alpha_beta_2, eps_2=eps_2, 
                 cps=cps, prior_mu=prior_mu,  ...)
  if (!is.null(file)) {
    saveRDS(out, file)
  }
  out
}

