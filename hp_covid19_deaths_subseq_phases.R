
# Re-initialise seed
set.seed(NULL)

# Prepare data and load functions
source('functions/prepare_data_deaths.R')
source('functions/functions_data_prep.R')
source('functions/functions_mcmc_algs.R')
source('functions/functions_model.R')

load("cov_mat_deaths_subseq_v2.RData")
load("mcmc_pars_subseq_v2.RData")


# MCMC parameters  

# Hyper parameters for mu
prior_mu = "gamma"
a_mu = 5
b_mu = 1

run_subseq_waves = function(country, seed_id, m_start=3){
        
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
        
        input = list(param=param_init, train=NULL, test=NULL)
        
        par_algo = list()
        par_algo$N = 60000
        par_algo$burnin = 20000
        par_algo$thin = 10
        par_algo$cps = pars$cps
        par_algo$s_mu = pars$s_mu
        par_algo$eps = pars$eps
        par_algo$g_alpha_beta = cov_mat_deaths[[country]]

        results = mala_cp(seed_id=seed_id, country, input, par_algo, paste0("mala_deaths_subseq_"), prior_mu, m_start)
        
        return(results)
}

mala_France_deaths_gamma5_1=run_subseq_waves("France", 10)
mala_Germany_deaths_gamma5_1=run_subseq_waves("Germany", 11)
mala_Italy_deaths_gamma5_1=run_subseq_waves("Italy", 12)
mala_Spain_deaths_gamma5_1=run_subseq_waves("Spain", 13)
mala_US_deaths_gamma5_1=run_subseq_waves("US", 14)
mala_Sweden_deaths_gamma5_1=run_subseq_waves("Sweden", 15)
mala_UK_deaths_gamma5_1=run_subseq_waves("UK", 16)
mala_Brazil_deaths_gamma5_1=run_subseq_waves("Brazil", 17, m_start = 2)
mala_India_deaths_gamma5_1=run_subseq_waves("India", 18, m_start = 2)
mala_China_deaths_gamma5_1=run_subseq_waves("China", 19)
