
######## Perform missing data interpolation ######## 

# Re-initialise seed
set.seed(NULL)

# Prepare data and load functions
source('functions/prepare_data_deaths.R')
source('functions/functions_data_prep.R')
source('functions/functions_mcmc_algs.R')
source('functions/functions_model.R')

load("cov_mat_deaths.RData")
load("mcmc_pars.RData")

# Set some initial parameters
burnin = 20000
N = 60000
thin = 10


# Hyper parameters for mu (using lognormal(1,1) prior)
prior_mu = "lognormal"
s_lmu = 1
mean_lmu = log(1) - (s_lmu^2)/2

run_lognormal1_1_interpol = function(seed_id, country){
  for (i in 1:5){
    assign("pars", get(paste0(country, "_pars")))
    results = mala_cp(seed_id=seed_id+(10*(i-1)), country, N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, 
                            pars$s_mu_1, cov_mat_deaths_1[[country]], pars$eps_1, 
                            pars$s_mu_2, cov_mat_deaths_2[[country]], pars$eps_2, 
                            paste0("mala_deaths_", prior_mu, "1_1_interpol", i), pars$cps, 
                            prior_mu, interpol=TRUE)
    assign(paste0("mala_", country, "_deaths_lognormal1_1_interpol", i), results, envir=.GlobalEnv)
  }
}
run_lognormal1_1_interpol(210, "France")
run_lognormal1_1_interpol(211, "Germany")
run_lognormal1_1_interpol(212, "Italy")
run_lognormal1_1_interpol(213, "Spain")
run_lognormal1_1_interpol(214, "US")
run_lognormal1_1_interpol(215, "Sweden")
run_lognormal1_1_interpol(216, "UK")
run_lognormal1_1_interpol(217, "Brazil")
run_lognormal1_1_interpol(218, "India")
run_lognormal1_1_interpol(219, "China")


# Hyper parameters for mu (using lognormal(5,1.5) prior)
prior_mu = "lognormal"
s_lmu = 1.5
mean_lmu = log(5) - (s_lmu^2)/2

run_lognormal5_1.5_interpol = function(seed_id, country){
        for (i in 1:5){
                assign("pars", get(paste0(country, "_pars")))
                results = mala_cp(seed_id=seed_id+(10*(i-1)), country, N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, 
                                        pars$s_mu_1, cov_mat_deaths_1[[country]], pars$eps_1, 
                                        pars$s_mu_2, cov_mat_deaths_2[[country]], pars$eps_2, 
                                        paste0("mala_deaths_", prior_mu, "5_1.5_interpol", i), pars$cps, 
                                        prior_mu, interpol=TRUE)
                assign(paste0("mala_", country, "_deaths_lognormal5_1.5_interpol", i), results, envir=.GlobalEnv)
        }
}
run_lognormal5_1.5_interpol(210, "France")
run_lognormal5_1.5_interpol(211, "Germany")
run_lognormal5_1.5_interpol(212, "Italy")
run_lognormal5_1.5_interpol(213, "Spain")
run_lognormal5_1.5_interpol(214, "US")
run_lognormal5_1.5_interpol(215, "Sweden")
run_lognormal5_1.5_interpol(216, "UK")
run_lognormal5_1.5_interpol(217, "Brazil")
run_lognormal5_1.5_interpol(218, "India")
run_lognormal5_1.5_interpol(219, "China")


# Hyper parameters for mu (using gamma(2,2) prior)
prior_mu = "gamma"
a_mu = 2
b_mu = 2

run_gamma2_2_interpol = function(seed_id, country){
        for (i in 1:5){
                assign("pars", get(paste0(country, "_pars")))
                results = mala_cp(seed_id=seed_id+(10*(i-1)), country, N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5,
                                        pars$s_mu_1, cov_mat_deaths_1[[country]], pars$eps_1,
                                        pars$s_mu_2, cov_mat_deaths_2[[country]], pars$eps_2,
                                        paste0("mala_deaths_", prior_mu, "2_2_interpol", i), pars$cps,
                                        prior_mu, interpol=TRUE)
                assign(paste0("mala_", country, "_deaths_gamma2_2_interpol", i), results, envir=.GlobalEnv)
        }
}
run_gamma2_2_interpol(210, "France")
run_gamma2_2_interpol(211, "Germany")
run_gamma2_2_interpol(212, "Italy")
run_gamma2_2_interpol(213, "Spain")
run_gamma2_2_interpol(214, "US")
run_gamma2_2_interpol(215, "Sweden")
run_gamma2_2_interpol(216, "UK")
run_gamma2_2_interpol(217, "Brazil")
run_gamma2_2_interpol(218, "India")
run_gamma2_2_interpol(219, "China")



# Hyper parameters for mu (using gamma(5,1) prior)
prior_mu = "gamma"
a_mu = 5
b_mu = 1

run_gamma5_1_interpol = function(seed_id, country){
        for (i in 1:5){
                assign("pars", get(paste0(country, "_pars")))
                results = mala_cp(seed_id=seed_id+(10*(i-1)), country, N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, 
                                        pars$s_mu_1, cov_mat_deaths_1[[country]], pars$eps_1, 
                                        pars$s_mu_2, cov_mat_deaths_2[[country]], pars$eps_2, 
                                        paste0("mala_deaths_", prior_mu, "5_1_interpol", i), pars$cps, 
                                        prior_mu, interpol=TRUE)
                assign(paste0("mala_", country, "_deaths_gamma5_1_interpol", i), results, envir=.GlobalEnv)
        }
}
run_gamma5_1_interpol(210, "France")
run_gamma5_1_interpol(211, "Germany")
run_gamma5_1_interpol(212, "Italy")
run_gamma5_1_interpol(213, "Spain")
run_gamma5_1_interpol(214, "US")
run_gamma5_1_interpol(215, "Sweden")
run_gamma5_1_interpol(216, "UK")
run_gamma5_1_interpol(217, "Brazil")
run_gamma5_1_interpol(218, "India")
run_gamma5_1_interpol(219, "China")


# Hyper parameters for mu (using uniform prior)
prior_mu = "uniform"

run_uniform_interpol = function(seed_id, country){
        for (i in 1:5){
                assign("pars", get(paste0(country, "_pars")))
                results = mala_cp(seed_id=seed_id+(10*(i-1)), country, N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, 
                                        pars$s_mu_1, cov_mat_deaths_1[[country]], pars$eps_1, 
                                        pars$s_mu_2, cov_mat_deaths_2[[country]], pars$eps_2, 
                                        paste0("mala_deaths_", prior_mu, "_interpol", i), pars$cps, 
                                        prior_mu, interpol=TRUE)
                assign(paste0("mala_", country, "_deaths_uniform_interpol", i), results, envir=.GlobalEnv)
        }
}

run_uniform_interpol(210, "France")
run_uniform_interpol(211, "Germany")
run_uniform_interpol(212, "Italy")
run_uniform_interpol(213, "Spain")
run_uniform_interpol(214, "US")
run_uniform_interpol(215, "Sweden")
run_uniform_interpol(216, "UK")
run_uniform_interpol(217, "Brazil")
run_uniform_interpol(218, "India")
run_uniform_interpol(219, "China")

