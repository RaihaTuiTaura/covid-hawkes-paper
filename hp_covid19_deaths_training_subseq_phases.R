
# Re-initialise seed
set.seed(NULL)

# Prepare data and load functions
source('functions/prepare_data_deaths.R')
source('functions/functions_data_prep.R')
source('functions/functions_mcmc_algs.R')
source('functions/functions_model.R')

load("cov_mat_deaths_subseq_v2.RData")
load("mcmc_pars_subseq_v2.RData")

# Hyper parameters for mu
prior_mu = "gamma"
a_mu = 5
b_mu = 1

# Function to run through different training/testing combos
oos_preds = function(country, seed_id, m_start=3){

  pars = get(paste0(country, "_pars"))
  dim=length(pars$cps)
  
  param_init = list()
  param_init$mu = rep(0,dim)
  for (m in 1:dim){
    param_init$mu[m] = 1
  }
  
  param_init$alpha = rep(0,dim)
  for (m in 1:dim){
    param_init$alpha[m] = 1
  }
  
  param_init$beta = rep(0,dim)
  for (m in 1:dim){
    param_init$beta[m] = 0.5
  }
  
  
  breaks = list()
  if (country=="Brazil"){
    breaks$new_end = 74
    breaks$cp_ind = 1
  } else if (country=="China"){
    breaks$new_end = 81
    breaks$cp_ind = 2
  } else if (country=="India"){
    breaks$new_end = 81
    breaks$cp_ind = 1
  } else if (country=="Spain"){
    breaks$new_end = 101
    breaks$cp_ind = 2
  } else if (country=="US"){
    breaks$new_end = 109
    breaks$cp_ind = 2
  } else if (country=="France"){
    breaks$new_end = 141
    breaks$cp_ind = 2
  } else if (country=="Germany"){
    breaks$new_end = 135
    breaks$cp_ind = 2
  } else if (country=="Italy"){
    breaks$new_end = 151
    breaks$cp_ind = 2
  } else if (country=="Sweden"){
    breaks$new_end = 130
    breaks$cp_ind = 2
  } else if (country=="UK"){
    breaks$new_end = 136
    breaks$cp_ind = 2
  } 
  
  test_size=10
  train = lapply((m_start-1):(dim-1), function(x) { if ((pars$cps[(x+1)]-pars$cps[x]-15) < test_size){
                                                    c()
                                                  } else {
                                                    if ((x+1) == breaks$cp_ind){
                                                      seq((pars$cps[x]+15), (breaks$new_end-test_size), test_size)
                                                    } else {
                                                      seq((pars$cps[x]+15), (pars$cps[(x+1)]-test_size), test_size)
                                                    }
                                                  }})
  train = unlist(train)
  input = list(param=param_init, test=test_size)

  par_algo = list()
  par_algo$N = 60000
  par_algo$burnin = 20000
  par_algo$thin = 10
  par_algo$cps = pars$cps
  par_algo$s_mu = pars$s_mu
  par_algo$eps = pars$eps
  par_algo$g_alpha_beta = cov_mat_deaths[[country]]
  
  get_pred_ints = lapply(1:length(train), function(x){input$train = train[x]; mala_cp(seed_id=seed_id, country, input, par_algo, 
                                                                              paste0("mala_deaths_subseq_oos_"), prior_mu, m_start)})
  
  return(list(get_pred_ints, train))

}

oos_France_gamma5_1 = oos_preds("France", 110)
oos_Germany_gamma5_1 = oos_preds("Germany", 111)
oos_Italy_gamma5_1 = oos_preds("Italy", 112)
oos_Spain_gamma5_1 = oos_preds("Spain", 113)
oos_US_gamma5_1 = oos_preds("US", 114)
oos_Sweden_gamma5_1 = oos_preds("Sweden", 115)
oos_UK_gamma5_1 = oos_preds("UK", 116)
oos_Brazil_gamma5_1 = oos_preds("Brazil", 117, m_start=2)
oos_India_gamma5_1 = oos_preds("India", 118, m_start=2)
oos_China_gamma5_1 = oos_preds("China", 119)
