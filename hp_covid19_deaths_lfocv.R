
######## Perform PSIS-LFO-CV (https://github.com/paul-buerkner/LFO-CV-paper) ######## 

# Re-initialise seed
set.seed(NULL)

# Prepare data and load functions
source('functions/prepare_data_deaths.R')
source('functions/functions_data_prep.R')
source('functions/functions_mcmc_algs.R')
source('functions/functions_model.R')
source('functions/functions_cross_validation.R')

load("cov_mat_deaths.RData")
load("mcmc_pars.RData")


parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

library(tidyverse)
library(brms)

# Set some initial parameters
k_thres = 0.7
burnin = 20000
N_iters = 60000
thin = 10


# Lognormal(mu=1,s=1) prior
s_lmu = 1
mean_lmu = log(1) - (s_lmu^2)/2

lfo_cv_lognormal1_1 = function(seed_id, country, L){
  assign("mcmc_pars", get(paste0(country, "_pars")))
  assign("cov_mat_1", get("cov_mat_deaths_1")[[country]])
  assign("cov_mat_2", get("cov_mat_deaths_2")[[country]])
  elpds <- compute_lfo(
    type = "approx", M = 4, L = L, 
    k_thres = k_thres, mode = "forward", factorize = TRUE,
    file = paste0(country, "_elpd_4sap_approx_lognormal1_1.rds"), 
    seed_id=seed_id, country=country, N_iters=N_iters, burnin=burnin, thin=thin, 
    mu1_0=1, alpha1_0=1, beta1_0=0.5,mu2_0=1, alpha2_0=1, beta2_0=0.5,
    s_mu_1=mcmc_pars$s_mu_1, g_alpha_beta_1=cov_mat_1, eps_1=mcmc_pars$eps_1, 
    s_mu_2=mcmc_pars$s_mu_2, g_alpha_beta_2=cov_mat_2, eps_2=mcmc_pars$eps_2, 
    cps=mcmc_pars$cps, prior_mu="lognormal"
  )
  return(elpds)
}
France_elpd_4sap_approx_lognormal1_1 = lfo_cv_lognormal1_1(seed_id=300, "France", 20)
Italy_elpd_4sap_approx_lognormal1_1 = lfo_cv_lognormal1_1(seed_id=301, "Italy", 20)
Spain_elpd_4sap_approx_lognormal1_1 = lfo_cv_lognormal1_1(seed_id=302, "Spain", 20)
Germany_elpd_4sap_approx_lognormal1_1 = lfo_cv_lognormal1_1(seed_id=303, "Germany", 20)
US_elpd_4sap_approx_lognormal1_1 = lfo_cv_lognormal1_1(seed_id=304, "US", 20)
Brazil_elpd_4sap_approx_lognormal1_1 = lfo_cv_lognormal1_1(seed_id=305, "Brazil", 20)
India_elpd_4sap_approx_lognormal1_1 = lfo_cv_lognormal1_1(seed_id=306, "India", 20)
Sweden_elpd_4sap_approx_lognormal1_1 = lfo_cv_lognormal1_1(seed_id=307, "Sweden", 20)
UK_elpd_4sap_approx_lognormal1_1 = lfo_cv_lognormal1_1(seed_id=308, "UK", 20)
China_elpd_4sap_approx_lognormal1_1 = lfo_cv_lognormal1_1(seed_id=309, "China", 20)


# Lognormal(mu=5,s=1.5) prior
s_lmu = 1.5
mean_lmu = log(5) - (s_lmu^2)/2

lfo_cv_lognormal5_1.5 = function(seed_id, country, L){
  assign("mcmc_pars", get(paste0(country, "_pars")))
  assign("cov_mat_1", get("cov_mat_deaths_1")[[country]])
  assign("cov_mat_2", get("cov_mat_deaths_2")[[country]])
  elpds <- compute_lfo(
    type = "approx", M = 4, L = L, 
    k_thres = k_thres, mode = "forward", factorize = TRUE,
    file = paste0(country, "_elpd_4sap_approx_lognormal5_1.5.rds"), 
    seed_id=seed_id, country=country, N_iters=N_iters, burnin=burnin, thin=thin, 
    mu1_0=1, alpha1_0=1, beta1_0=0.5,mu2_0=1, alpha2_0=1, beta2_0=0.5,
    s_mu_1=mcmc_pars$s_mu_1, g_alpha_beta_1=cov_mat_1, eps_1=mcmc_pars$eps_1, 
    s_mu_2=mcmc_pars$s_mu_2, g_alpha_beta_2=cov_mat_2, eps_2=mcmc_pars$eps_2, 
    cps=mcmc_pars$cps, prior_mu="lognormal"
  )
  return(elpds)
}
France_elpd_4sap_approx_lognormal5_1.5 = lfo_cv_lognormal5_1.5(seed_id=300, "France", 20)
Italy_elpd_4sap_approx_lognormal5_1.5 = lfo_cv_lognormal5_1.5(seed_id=301, "Italy", 20)
Spain_elpd_4sap_approx_lognormal5_1.5 = lfo_cv_lognormal5_1.5(seed_id=302, "Spain", 20)
Germany_elpd_4sap_approx_lognormal5_1.5 = lfo_cv_lognormal5_1.5(seed_id=303, "Germany", 20)
US_elpd_4sap_approx_lognormal5_1.5 = lfo_cv_lognormal5_1.5(seed_id=304, "US", 20)
Brazil_elpd_4sap_approx_lognormal5_1.5 = lfo_cv_lognormal5_1.5(seed_id=305, "Brazil", 20)
India_elpd_4sap_approx_lognormal5_1.5 = lfo_cv_lognormal5_1.5(seed_id=306, "India", 20)
Sweden_elpd_4sap_approx_lognormal5_1.5 = lfo_cv_lognormal5_1.5(seed_id=307, "Sweden", 20)
UK_elpd_4sap_approx_lognormal5_1.5 = lfo_cv_lognormal5_1.5(seed_id=308, "UK", 20)
China_elpd_4sap_approx_lognormal5_1.5 = lfo_cv_lognormal5_1.5(seed_id=309, "China", 20)


# Gamma(2,2) prior
a_mu=2
b_mu=2

lfo_cv_gamma2_2 = function(seed_id, country, L){
  assign("mcmc_pars", get(paste0(country, "_pars")))
  assign("cov_mat_1", get("cov_mat_deaths_1")[[country]])
  assign("cov_mat_2", get("cov_mat_deaths_2")[[country]])
  elpds <- compute_lfo(
    type = "approx", M = 4, L = L, 
    k_thres = k_thres, mode = "forward", factorize = TRUE,
    file = paste0(country, "_elpd_4sap_approx_gamma2_2.rds"), 
    seed_id=seed_id, country=country, N_iters=N_iters, burnin=burnin, thin=thin, 
    mu1_0=1, alpha1_0=1, beta1_0=0.5,mu2_0=1, alpha2_0=1, beta2_0=0.5,
    s_mu_1=mcmc_pars$s_mu_1, g_alpha_beta_1=cov_mat_1, eps_1=mcmc_pars$eps_1, 
    s_mu_2=mcmc_pars$s_mu_2, g_alpha_beta_2=cov_mat_2, eps_2=mcmc_pars$eps_2, 
    cps=mcmc_pars$cps, prior_mu="gamma"
  )
  return(elpds)
}
France_elpd_4sap_approx_gamma2_2 = lfo_cv_gamma2_2(seed_id=300, "France", 20)
Italy_elpd_4sap_approx_gamma2_2 = lfo_cv_gamma2_2(seed_id=301, "Italy", 20)
Spain_elpd_4sap_approx_gamma2_2 = lfo_cv_gamma2_2(seed_id=302, "Spain", 20)
Germany_elpd_4sap_approx_gamma2_2 = lfo_cv_gamma2_2(seed_id=303, "Germany", 20)
US_elpd_4sap_approx_gamma2_2 = lfo_cv_gamma2_2(seed_id=304, "US", 20)
Brazil_elpd_4sap_approx_gamma2_2 = lfo_cv_gamma2_2(seed_id=305, "Brazil", 20)
India_elpd_4sap_approx_gamma2_2 = lfo_cv_gamma2_2(seed_id=306, "India", 20)
Sweden_elpd_4sap_approx_gamma2_2 = lfo_cv_gamma2_2(seed_id=307, "Sweden", 20)
UK_elpd_4sap_approx_gamma2_2 = lfo_cv_gamma2_2(seed_id=308, "UK", 20)
China_elpd_4sap_approx_gamma2_2 = lfo_cv_gamma2_2(seed_id=309, "China", 20)


# Gamma (5,1) prior
a_mu=5
b_mu=1

lfo_cv_gamma5_1 = function(seed_id, country, L){
  assign("mcmc_pars", get(paste0(country, "_pars")))
  assign("cov_mat_1", get("cov_mat_deaths_1")[[country]])
  assign("cov_mat_2", get("cov_mat_deaths_2")[[country]])
  elpds <- compute_lfo(
    type = "approx", M = 4, L = L, 
    k_thres = k_thres, mode = "forward", factorize = TRUE,
    file = paste0(country, "_elpd_4sap_approx_gamma5_1.rds"), 
    seed_id=seed_id, country=country, N_iters=N_iters, burnin=burnin, thin=thin, 
    mu1_0=1, alpha1_0=1, beta1_0=0.5,mu2_0=1, alpha2_0=1, beta2_0=0.5,
    s_mu_1=mcmc_pars$s_mu_1, g_alpha_beta_1=cov_mat_1, eps_1=mcmc_pars$eps_1, 
    s_mu_2=mcmc_pars$s_mu_2, g_alpha_beta_2=cov_mat_2, eps_2=mcmc_pars$eps_2,
    cps=mcmc_pars$cps, prior_mu="gamma"
  )
  return(elpds)
}
France_elpd_4sap_approx_gamma5_1 = lfo_cv_gamma5_1(seed_id=300, "France", 20)
Italy_elpd_4sap_approx_gamma5_1 = lfo_cv_gamma5_1(seed_id=301, "Italy", 20)
Spain_elpd_4sap_approx_gamma5_1 = lfo_cv_gamma5_1(seed_id=302, "Spain", 20)
Germany_elpd_4sap_approx_gamma5_1 = lfo_cv_gamma5_1(seed_id=303, "Germany", 20)
US_elpd_4sap_approx_gamma5_1 = lfo_cv_gamma5_1(seed_id=304, "US", 20)
Brazil_elpd_4sap_approx_gamma5_1 = lfo_cv_gamma5_1(seed_id=305, "Brazil", 20)
India_elpd_4sap_approx_gamma5_1 = lfo_cv_gamma5_1(seed_id=306, "India", 20)
Sweden_elpd_4sap_approx_gamma5_1 = lfo_cv_gamma5_1(seed_id=307, "Sweden", 20)
UK_elpd_4sap_approx_gamma5_1 = lfo_cv_gamma5_1(seed_id=308, "UK", 20)
China_elpd_4sap_approx_gamma5_1 = lfo_cv_gamma5_1(seed_id=309, "China", 20)


# Uniform prior
lfo_cv_uniform = function(seed_id, country, L){
  assign("mcmc_pars", get(paste0(country, "_pars")))
  assign("cov_mat_1", get("cov_mat_deaths_1")[[country]])
  assign("cov_mat_2", get("cov_mat_deaths_2")[[country]])
  elpds <- compute_lfo(
    type = "approx", M = 4, L = L, 
    k_thres = k_thres, mode = "forward", factorize = TRUE,
    file = paste0(country, "_elpd_4sap_approx_uniform.rds"), 
    seed_id=seed_id, country=country, N_iters=N_iters, burnin=burnin, thin=thin, 
    mu1_0=1, alpha1_0=1, beta1_0=0.5,mu2_0=1, alpha2_0=1, beta2_0=0.5, 
    s_mu_1=mcmc_pars$s_mu_1, g_alpha_beta_1=cov_mat_1, eps_1=mcmc_pars$eps_1, 
    s_mu_2=mcmc_pars$s_mu_2, g_alpha_beta_2=cov_mat_2, eps_2=mcmc_pars$eps_2,
    cps=mcmc_pars$cps, prior_mu="uniform"
  )
  return(elpds)
}
France_elpd_4sap_approx_uniform = lfo_cv_uniform(seed_id=300, "France", 20)
Italy_elpd_4sap_approx_uniform = lfo_cv_uniform(seed_id=301, "Italy", 20)
Spain_elpd_4sap_approx_uniform = lfo_cv_uniform(seed_id=302, "Spain", 20)
Germany_elpd_4sap_approx_uniform = lfo_cv_uniform(seed_id=303, "Germany", 20)
US_elpd_4sap_approx_uniform = lfo_cv_uniform(seed_id=304, "US", 20)
Brazil_elpd_4sap_approx_uniform = lfo_cv_uniform(seed_id=305, "Brazil", 20)
India_elpd_4sap_approx_uniform = lfo_cv_uniform(seed_id=306, "India", 20)
Sweden_elpd_4sap_approx_uniform = lfo_cv_uniform(seed_id=307, "Sweden", 20)
UK_elpd_4sap_approx_uniform = lfo_cv_uniform(seed_id=308, "UK", 20)
China_elpd_4sap_approx_uniform = lfo_cv_uniform(seed_id=309, "China", 20)


# Summarise results
France_gamma2_2_sum = round(summarize_elpds(France_elpd_4sap_approx_gamma2_2), digits=2)
France_gamma5_1_sum = round(summarize_elpds(France_elpd_4sap_approx_gamma5_1), digits=2)
France_lognormal1_1_sum = round(summarize_elpds(France_elpd_4sap_approx_lognormal1_1), digits=2)
France_lognormal5_1.5_sum = round(summarize_elpds(France_elpd_4sap_approx_lognormal5_1.5), digits=2)
France_uniform_sum = round(summarize_elpds(France_elpd_4sap_approx_uniform), digits=2)

Italy_gamma2_2_sum = round(summarize_elpds(Italy_elpd_4sap_approx_gamma2_2), digits=2)
Italy_gamma5_1_sum = round(summarize_elpds(Italy_elpd_4sap_approx_gamma5_1), digits=2)
Italy_lognormal1_1_sum = round(summarize_elpds(Italy_elpd_4sap_approx_lognormal1_1), digits=2)
Italy_lognormal5_1.5_sum = round(summarize_elpds(Italy_elpd_4sap_approx_lognormal5_1.5), digits=2)
Italy_uniform_sum = round(summarize_elpds(Italy_elpd_4sap_approx_uniform), digits=2)

Spain_gamma2_2_sum = round(summarize_elpds(Spain_elpd_4sap_approx_gamma2_2), digits=2)
Spain_gamma5_1_sum = round(summarize_elpds(Spain_elpd_4sap_approx_gamma5_1), digits=2)
Spain_lognormal1_1_sum = round(summarize_elpds(Spain_elpd_4sap_approx_lognormal1_1), digits=2)
Spain_lognormal5_1.5_sum = round(summarize_elpds(Spain_elpd_4sap_approx_lognormal5_1.5), digits=2)
Spain_uniform_sum = round(summarize_elpds(Spain_elpd_4sap_approx_uniform), digits=2)

Germany_gamma2_2_sum = round(summarize_elpds(Germany_elpd_4sap_approx_gamma2_2), digits=2)
Germany_gamma5_1_sum = round(summarize_elpds(Germany_elpd_4sap_approx_gamma5_1), digits=2)
Germany_lognormal1_1_sum = round(summarize_elpds(Germany_elpd_4sap_approx_lognormal1_1), digits=2)
Germany_lognormal5_1.5_sum = round(summarize_elpds(Germany_elpd_4sap_approx_lognormal5_1.5), digits=2)
Germany_uniform_sum = round(summarize_elpds(Germany_elpd_4sap_approx_uniform), digits=2)

US_gamma2_2_sum = round(summarize_elpds(US_elpd_4sap_approx_gamma2_2), digits=2)
US_gamma5_1_sum = round(summarize_elpds(US_elpd_4sap_approx_gamma5_1), digits=2)
US_lognormal1_1_sum = round(summarize_elpds(US_elpd_4sap_approx_lognormal1_1), digits=2)
US_lognormal5_1.5_sum = round(summarize_elpds(US_elpd_4sap_approx_lognormal5_1.5), digits=2)
US_uniform_sum = round(summarize_elpds(US_elpd_4sap_approx_uniform), digits=2)

Brazil_gamma2_2_sum = round(summarize_elpds(Brazil_elpd_4sap_approx_gamma2_2), digits=2)
Brazil_gamma5_1_sum = round(summarize_elpds(Brazil_elpd_4sap_approx_gamma5_1), digits=2)
Brazil_lognormal1_1_sum = round(summarize_elpds(Brazil_elpd_4sap_approx_lognormal1_1), digits=2)
Brazil_lognormal5_1.5_sum = round(summarize_elpds(Brazil_elpd_4sap_approx_lognormal5_1.5), digits=2)
Brazil_uniform_sum = round(summarize_elpds(Brazil_elpd_4sap_approx_uniform), digits=2)

India_gamma2_2_sum = round(summarize_elpds(India_elpd_4sap_approx_gamma2_2), digits=2)
India_gamma5_1_sum = round(summarize_elpds(India_elpd_4sap_approx_gamma5_1), digits=2)
India_lognormal1_1_sum = round(summarize_elpds(India_elpd_4sap_approx_lognormal1_1), digits=2)
India_lognormal5_1.5_sum = round(summarize_elpds(India_elpd_4sap_approx_lognormal5_1.5), digits=2)
India_uniform_sum = round(summarize_elpds(India_elpd_4sap_approx_uniform), digits=2)

Sweden_gamma2_2_sum = round(summarize_elpds(Sweden_elpd_4sap_approx_gamma2_2), digits=2)
Sweden_gamma5_1_sum = round(summarize_elpds(Sweden_elpd_4sap_approx_gamma5_1), digits=2)
Sweden_lognormal1_1_sum = round(summarize_elpds(Sweden_elpd_4sap_approx_lognormal1_1), digits=2)
Sweden_lognormal5_1.5_sum = round(summarize_elpds(Sweden_elpd_4sap_approx_lognormal5_1.5), digits=2)
Sweden_uniform_sum = round(summarize_elpds(Sweden_elpd_4sap_approx_uniform), digits=2)

UK_gamma2_2_sum = round(summarize_elpds(UK_elpd_4sap_approx_gamma2_2), digits=2)
UK_gamma5_1_sum = round(summarize_elpds(UK_elpd_4sap_approx_gamma5_1), digits=2)
UK_lognormal1_1_sum = round(summarize_elpds(UK_elpd_4sap_approx_lognormal1_1), digits=2)
UK_lognormal5_1.5_sum = round(summarize_elpds(UK_elpd_4sap_approx_lognormal5_1.5), digits=2)
UK_uniform_sum = round(summarize_elpds(UK_elpd_4sap_approx_uniform), digits=2)

China_gamma2_2_sum = round(summarize_elpds(China_elpd_4sap_approx_gamma2_2), digits=2)
China_gamma5_1_sum = round(summarize_elpds(China_elpd_4sap_approx_gamma5_1), digits=2)
China_lognormal1_1_sum = round(summarize_elpds(China_elpd_4sap_approx_lognormal1_1), digits=2)
China_lognormal5_1.5_sum = round(summarize_elpds(China_elpd_4sap_approx_lognormal5_1.5), digits=2)
China_uniform_sum = round(summarize_elpds(China_elpd_4sap_approx_uniform), digits=2)

