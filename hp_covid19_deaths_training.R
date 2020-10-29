
######## Out-of-sample validation ######## 

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

test_size_train1 = 5
test_size_train2 = 10


max_Ts = list(France=138, Germany=131, Italy=149, Spain=101, US=109,
              Sweden=126, UK=134, Brazil=125, India=81, China=81)


# Function to run through different training/testing combos
oos_preds = function(country, seed_id, prior_desc){
  
  if (prior_desc == "gamma2_2"){
    assign("prior_mu", "gamma", envir = .GlobalEnv)
    assign("a_mu", 2, envir = .GlobalEnv)
    assign("b_mu", 2, envir = .GlobalEnv)
  } else if (prior_desc == "gamma5_1"){
    assign("prior_mu", "gamma", envir = .GlobalEnv)
    assign("a_mu", 5, envir = .GlobalEnv)
    assign("b_mu", 1, envir = .GlobalEnv)
  } else if (prior_desc == "lognormal1_1"){
    assign("prior_mu", "lognormal", envir = .GlobalEnv)
    assign("s_lmu", 1, envir = .GlobalEnv)
    assign("mean_lmu", log(1) - (s_lmu^2)/2, envir = .GlobalEnv)  
  } else if (prior_desc == "lognormal5_1.5"){
    assign("prior_mu", "lognormal", envir = .GlobalEnv)
    assign("s_lmu", 1.5, envir = .GlobalEnv)
    assign("mean_lmu", log(5) - (s_lmu^2)/2, envir = .GlobalEnv)  
  } else if (prior_desc == "uniform"){
    assign("prior_mu", "uniform", envir = .GlobalEnv)
  }
  
  pars = get(paste0(country, "_pars"))
  train_1 = seq(15, pars$cps, test_size_train1)
  cov_mat_1 = get("cov_mat_deaths_1")[[country]]
  cov_mat_2 = get("cov_mat_deaths_2")[[country]]
  get_pred_ints_1 = lapply(train_1, function(x){mala_cp(seed_id=seed_id, country, N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, pars$s_mu_1, cov_mat_1, pars$eps_1,
                                                             pars$s_mu_2, cov_mat_2, pars$eps_2, paste0("mala_deaths_oos_", prior_desc), pars$cps, prior_mu, train=x, test=test_size_train1)})
  
  if (!(country %in% c("India", "Brazil"))){
    train_2 = seq((pars$cps+15), (max_Ts[[country]]-test_size_train2), test_size_train2)
    get_pred_ints_2 = lapply(train_2, function(x){mala_cp(seed_id=(seed_id+10), country, N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, pars$s_mu_1, cov_mat_1, pars$eps_1,
                                                               pars$s_mu_2, cov_mat_2, pars$eps_2, paste0("mala_deaths_oos_", prior_desc), pars$cps, prior_mu, train=x, test=test_size_train2)})
    training_seq = c(train_1, train_2)
    
  } else {
    training_seq = train_1
  }
  
  if (!(country %in% c("India", "Brazil"))){
    return(list(get_pred_ints_1, get_pred_ints_2, train_1, train_2))
  } else {
    return(list(get_pred_ints_1, train_1))
  }
}

# Lognormal(mu=1,s=1) prior
oos_France_lognormal1_1 = oos_preds("France", 110, "lognormal1_1")
oos_Germany_lognormal1_1 = oos_preds("Germany", 111, "lognormal1_1")
oos_Italy_lognormal1_1 = oos_preds("Italy", 112, "lognormal1_1")
oos_Spain_lognormal1_1 = oos_preds("Spain", 113, "lognormal1_1")
oos_US_lognormal1_1 = oos_preds("US", 114, "lognormal1_1")
oos_Sweden_lognormal1_1 = oos_preds("Sweden", 115, "lognormal1_1")
oos_UK_lognormal1_1 = oos_preds("UK", 116, "lognormal1_1")
oos_Brazil_lognormal1_1 = oos_preds("Brazil", 117, "lognormal1_1")
oos_India_lognormal1_1 = oos_preds("India", 118, "lognormal1_1")
oos_China_lognormal1_1 = oos_preds("China", 119, "lognormal1_1")

# Lognormal(mu=5,s=1.5) prior
oos_France_lognormal5_1.5 = oos_preds("France", 110, "lognormal5_1.5")
oos_Germany_lognormal5_1.5 = oos_preds("Germany", 111, "lognormal5_1.5")
oos_Italy_lognormal5_1.5 = oos_preds("Italy", 112, "lognormal5_1.5")
oos_Spain_lognormal5_1.5 = oos_preds("Spain", 113, "lognormal5_1.5")
oos_US_lognormal5_1.5 = oos_preds("US", 114, "lognormal5_1.5")
oos_Sweden_lognormal5_1.5 = oos_preds("Sweden", 115, "lognormal5_1.5")
oos_UK_lognormal5_1.5 = oos_preds("UK", 116, "lognormal5_1.5")
oos_Brazil_lognormal5_1.5 = oos_preds("Brazil", 117, "lognormal5_1.5")
oos_India_lognormal5_1.5 = oos_preds("India", 118, "lognormal5_1.5")
oos_China_lognormal5_1.5 = oos_preds("China", 119, "lognormal5_1.5")

# Gamma(2,2) prior
oos_France_gamma2_2 = oos_preds("France", 110, "gamma2_2")
oos_Germany_gamma2_2 = oos_preds("Germany", 111, "gamma2_2")
oos_Italy_gamma2_2 = oos_preds("Italy", 112, "gamma2_2")
oos_Spain_gamma2_2 = oos_preds("Spain", 113, "gamma2_2")
oos_US_gamma2_2 = oos_preds("US", 114, "gamma2_2")
oos_Sweden_gamma2_2 = oos_preds("Sweden", 115, "gamma2_2")
oos_UK_gamma2_2 = oos_preds("UK", 116, "gamma2_2")
oos_Brazil_gamma2_2 = oos_preds("Brazil", 117, "gamma2_2")
oos_India_gamma2_2 = oos_preds("India", 118, "gamma2_2")
oos_China_gamma2_2 = oos_preds("China", 119, "gamma2_2")

# Gamma (5,1) prior
oos_France_gamma5_1 = oos_preds("France", 110, "gamma5_1")
oos_Germany_gamma5_1 = oos_preds("Germany", 111, "gamma5_1")
oos_Italy_gamma5_1 = oos_preds("Italy", 112, "gamma5_1")
oos_Spain_gamma5_1 = oos_preds("Spain", 113, "gamma5_1")
oos_US_gamma5_1 = oos_preds("US", 114, "gamma5_1")
oos_Sweden_gamma5_1 = oos_preds("Sweden", 115, "gamma5_1")
oos_UK_gamma5_1 = oos_preds("UK", 116, "gamma5_1")
oos_Brazil_gamma5_1 = oos_preds("Brazil", 117, "gamma5_1")
oos_India_gamma5_1 = oos_preds("India", 118, "gamma5_1")
oos_China_gamma5_1 = oos_preds("China", 119, "gamma5_1")

# Uniform prior
oos_France_uniform = oos_preds("France", 110, "uniform")
oos_Germany_uniform = oos_preds("Germany", 111, "uniform")
oos_Italy_uniform = oos_preds("Italy", 112, "uniform")
oos_Spain_uniform = oos_preds("Spain", 113, "uniform")
oos_US_uniform = oos_preds("US", 114, "uniform")
oos_Sweden_uniform = oos_preds("Sweden", 115, "uniform")
oos_UK_uniform = oos_preds("UK", 116, "uniform")
oos_Brazil_uniform = oos_preds("Brazil", 117, "uniform")
oos_India_uniform = oos_preds("India", 118, "uniform")
oos_China_uniform = oos_preds("China", 119, "uniform")

