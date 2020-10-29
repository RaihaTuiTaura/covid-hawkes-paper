
######## Modelling COVID-19 using discrete-time Hawkes process ######## 

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

mala_France_deaths_lognormal1_1 = mala_cp(seed_id=10, "France",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, France_pars$s_mu_1, cov_mat_deaths_1$France, France_pars$eps_1, France_pars$s_mu_2, cov_mat_deaths_2$France, France_pars$eps_2, paste0("mala_deaths_", prior_mu, "1_1"), France_pars$cps, prior_mu)
mala_Germany_deaths_lognormal1_1 = mala_cp(seed_id=11, "Germany",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Germany_pars$s_mu_1, cov_mat_deaths_1$Germany, Germany_pars$eps_1, Germany_pars$s_mu_2, cov_mat_deaths_2$Germany, Germany_pars$eps_2, paste0("mala_deaths_", prior_mu, "1_1"),Germany_pars$cps, prior_mu)
mala_Italy_deaths_lognormal1_1 = mala_cp(seed_id=12, "Italy",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Italy_pars$s_mu_1, cov_mat_deaths_1$Italy, Italy_pars$eps_1, Italy_pars$s_mu_2, cov_mat_deaths_2$Italy, Italy_pars$eps_2, paste0("mala_deaths_", prior_mu, "1_1"), Italy_pars$cps, prior_mu)
mala_Spain_deaths_lognormal1_1 = mala_cp(seed_id=13, "Spain",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Spain_pars$s_mu_1, cov_mat_deaths_1$Spain, Spain_pars$eps_1, Spain_pars$s_mu_2, cov_mat_deaths_2$Spain, Spain_pars$eps_2, paste0("mala_deaths_", prior_mu, "1_1"), Spain_pars$cps, prior_mu)
mala_US_deaths_lognormal1_1 = mala_cp(seed_id=14, "US",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, US_pars$s_mu_1, cov_mat_deaths_1$US, US_pars$eps_1, US_pars$s_mu_2, cov_mat_deaths_2$US, US_pars$eps_2, paste0("mala_deaths_", prior_mu, "1_1"), US_pars$cps, prior_mu)
mala_Sweden_deaths_lognormal1_1 = mala_cp(seed_id=15, "Sweden",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Sweden_pars$s_mu_1, cov_mat_deaths_1$Sweden, Sweden_pars$eps_1, Sweden_pars$s_mu_2, cov_mat_deaths_2$Sweden, Sweden_pars$eps_2, paste0("mala_deaths_", prior_mu, "1_1"), Sweden_pars$cps, prior_mu)
mala_UK_deaths_lognormal1_1 = mala_cp(seed_id=16, "UK",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, UK_pars$s_mu_1, cov_mat_deaths_1$UK, UK_pars$eps_1, UK_pars$s_mu_2, cov_mat_deaths_2$UK, UK_pars$eps_2, paste0("mala_deaths_", prior_mu, "1_1"), UK_pars$cps, prior_mu)
mala_Brazil_deaths_lognormal1_1 = mala_cp(seed_id=17, "Brazil",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Brazil_pars$s_mu_1, cov_mat_deaths_1$Brazil, Brazil_pars$eps_1, Brazil_pars$s_mu_2, cov_mat_deaths_1$Brazil, Brazil_pars$eps_2, paste0("mala_deaths_", prior_mu, "1_1"), Brazil_pars$cps, prior_mu)
mala_India_deaths_lognormal1_1 = mala_cp(seed_id=18, "India",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, India_pars$s_mu_1, cov_mat_deaths_1$India, India_pars$eps_1, India_pars$s_mu_2, cov_mat_deaths_1$India, India_pars$eps_2, paste0("mala_deaths_", prior_mu, "1_1"), India_pars$cps, prior_mu)
mala_China_deaths_lognormal1_1 = mala_cp(seed_id=19, "China",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, China_pars$s_mu_1, cov_mat_deaths_1$China, China_pars$eps_1, China_pars$s_mu_2, cov_mat_deaths_2$China, China_pars$eps_2, paste0("mala_deaths_", prior_mu, "1_1"), China_pars$cps, prior_mu)


# Hyper parameters for mu (using lognormal(5,1.5) prior)
prior_mu = "lognormal"
s_lmu = 1.5
mean_lmu = log(5) - (s_lmu^2)/2

mala_France_deaths_lognormal5_1.5 = mala_cp(seed_id=10, "France",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, France_pars$s_mu_1, cov_mat_deaths_1$France, France_pars$eps_1, France_pars$s_mu_2, cov_mat_deaths_2$France, France_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1.5"), France_pars$cps, prior_mu)
mala_Germany_deaths_lognormal5_1.5 = mala_cp(seed_id=11, "Germany",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Germany_pars$s_mu_1, cov_mat_deaths_1$Germany, Germany_pars$eps_1, Germany_pars$s_mu_2, cov_mat_deaths_2$Germany, Germany_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1.5"),Germany_pars$cps, prior_mu)
mala_Italy_deaths_lognormal5_1.5 = mala_cp(seed_id=12, "Italy",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Italy_pars$s_mu_1, cov_mat_deaths_1$Italy, Italy_pars$eps_1, Italy_pars$s_mu_2, cov_mat_deaths_2$Italy, Italy_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1.5"), Italy_pars$cps, prior_mu)
mala_Spain_deaths_lognormal5_1.5 = mala_cp(seed_id=13, "Spain",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Spain_pars$s_mu_1, cov_mat_deaths_1$Spain, Spain_pars$eps_1, Spain_pars$s_mu_2, cov_mat_deaths_2$Spain, Spain_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1.5"), Spain_pars$cps, prior_mu)
mala_US_deaths_lognormal5_1.5 = mala_cp(seed_id=14, "US",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, US_pars$s_mu_1, cov_mat_deaths_1$US, US_pars$eps_1, US_pars$s_mu_2, cov_mat_deaths_2$US, US_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1.5"), US_pars$cps, prior_mu)
mala_Sweden_deaths_lognormal5_1.5 = mala_cp(seed_id=15, "Sweden",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Sweden_pars$s_mu_1, cov_mat_deaths_1$Sweden, Sweden_pars$eps_1, Sweden_pars$s_mu_2, cov_mat_deaths_2$Sweden, Sweden_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1.5"), Sweden_pars$cps, prior_mu)
mala_UK_deaths_lognormal5_1.5 = mala_cp(seed_id=16, "UK",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, UK_pars$s_mu_1, cov_mat_deaths_1$UK, UK_pars$eps_1, UK_pars$s_mu_2, cov_mat_deaths_2$UK, UK_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1.5"), UK_pars$cps, prior_mu)
mala_Brazil_deaths_lognormal5_1.5 = mala_cp(seed_id=17, "Brazil",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Brazil_pars$s_mu_1, cov_mat_deaths_1$Brazil, Brazil_pars$eps_1, Brazil_pars$s_mu_2, cov_mat_deaths_1$Brazil, Brazil_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1.5"), Brazil_pars$cps, prior_mu)
mala_India_deaths_lognormal5_1.5 = mala_cp(seed_id=18, "India",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, India_pars$s_mu_1, cov_mat_deaths_1$India, India_pars$eps_1, India_pars$s_mu_2, cov_mat_deaths_1$India, India_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1.5"), India_pars$cps, prior_mu)
mala_China_deaths_lognormal5_1.5 = mala_cp(seed_id=19, "China",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, China_pars$s_mu_1, cov_mat_deaths_1$China, China_pars$eps_1, China_pars$s_mu_2, cov_mat_deaths_2$China, China_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1.5"), China_pars$cps, prior_mu)


# Hyper parameters for mu (using gamma(2,2) prior)
prior_mu = "gamma"
a_mu = 2
b_mu = 2

mala_France_deaths_gamma2_2 = mala_cp(seed_id=10, "France",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, France_pars$s_mu_1, cov_mat_deaths_1$France, France_pars$eps_1, France_pars$s_mu_2, cov_mat_deaths_2$France, France_pars$eps_2, paste0("mala_deaths_", prior_mu, "2_2"), France_pars$cps, prior_mu)
mala_Germany_deaths_gamma2_2 = mala_cp(seed_id=11, "Germany",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Germany_pars$s_mu_1, cov_mat_deaths_1$Germany, Germany_pars$eps_1, Germany_pars$s_mu_2, cov_mat_deaths_2$Germany, Germany_pars$eps_2, paste0("mala_deaths_", prior_mu, "2_2"),Germany_pars$cps, prior_mu)
mala_Italy_deaths_gamma2_2 = mala_cp(seed_id=12, "Italy",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Italy_pars$s_mu_1, cov_mat_deaths_1$Italy, Italy_pars$eps_1, Italy_pars$s_mu_2, cov_mat_deaths_2$Italy, Italy_pars$eps_2, paste0("mala_deaths_", prior_mu, "2_2"), Italy_pars$cps, prior_mu)
mala_Spain_deaths_gamma2_2 = mala_cp(seed_id=13, "Spain",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Spain_pars$s_mu_1, cov_mat_deaths_1$Spain, Spain_pars$eps_1, Spain_pars$s_mu_2, cov_mat_deaths_2$Spain, Spain_pars$eps_2, paste0("mala_deaths_", prior_mu, "2_2"), Spain_pars$cps, prior_mu)
mala_US_deaths_gamma2_2 = mala_cp(seed_id=14, "US",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, US_pars$s_mu_1, cov_mat_deaths_1$US, US_pars$eps_1, US_pars$s_mu_2, cov_mat_deaths_2$US, US_pars$eps_2, paste0("mala_deaths_", prior_mu, "2_2"), US_pars$cps, prior_mu)
mala_Sweden_deaths_gamma2_2 = mala_cp(seed_id=15, "Sweden",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Sweden_pars$s_mu_1, cov_mat_deaths_1$Sweden, Sweden_pars$eps_1, Sweden_pars$s_mu_2, cov_mat_deaths_2$Sweden, Sweden_pars$eps_2, paste0("mala_deaths_", prior_mu, "2_2"), Sweden_pars$cps, prior_mu)
mala_UK_deaths_gamma2_2 = mala_cp(seed_id=16, "UK",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, UK_pars$s_mu_1, cov_mat_deaths_1$UK, UK_pars$eps_1, UK_pars$s_mu_2, cov_mat_deaths_2$UK, UK_pars$eps_2, paste0("mala_deaths_", prior_mu, "2_2"), UK_pars$cps, prior_mu)
mala_Brazil_deaths_gamma2_2 = mala_cp(seed_id=17, "Brazil",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Brazil_pars$s_mu_1, cov_mat_deaths_1$Brazil, Brazil_pars$eps_1, Brazil_pars$s_mu_2, cov_mat_deaths_1$Brazil, Brazil_pars$eps_2, paste0("mala_deaths_", prior_mu, "2_2"), Brazil_pars$cps, prior_mu)
mala_India_deaths_gamma2_2 = mala_cp(seed_id=18, "India",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, India_pars$s_mu_1, cov_mat_deaths_1$India, India_pars$eps_1, India_pars$s_mu_2, cov_mat_deaths_1$India, India_pars$eps_2, paste0("mala_deaths_", prior_mu, "2_2"), India_pars$cps, prior_mu)
mala_China_deaths_gamma2_2 = mala_cp(seed_id=19, "China",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, China_pars$s_mu_1, cov_mat_deaths_1$China, China_pars$eps_1, China_pars$s_mu_2, cov_mat_deaths_2$China, China_pars$eps_2, paste0("mala_deaths_", prior_mu, "2_2"), China_pars$cps, prior_mu)


# Hyper parameters for mu (using gamma(5,1) prior)
prior_mu = "gamma"
a_mu = 5
b_mu = 1

mala_France_deaths_gamma5_1 = mala_cp(seed_id=10, "France",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, France_pars$s_mu_1, cov_mat_deaths_1$France, France_pars$eps_1, France_pars$s_mu_2, cov_mat_deaths_2$France, France_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1"), France_pars$cps, prior_mu)
mala_Germany_deaths_gamma5_1 = mala_cp(seed_id=11, "Germany",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Germany_pars$s_mu_1, cov_mat_deaths_1$Germany, Germany_pars$eps_1, Germany_pars$s_mu_2, cov_mat_deaths_2$Germany, Germany_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1"),Germany_pars$cps, prior_mu)
mala_Italy_deaths_gamma5_1 = mala_cp(seed_id=12, "Italy",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Italy_pars$s_mu_1, cov_mat_deaths_1$Italy, Italy_pars$eps_1, Italy_pars$s_mu_2, cov_mat_deaths_2$Italy, Italy_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1"), Italy_pars$cps, prior_mu)
mala_Spain_deaths_gamma5_1 = mala_cp(seed_id=13, "Spain",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Spain_pars$s_mu_1, cov_mat_deaths_1$Spain, Spain_pars$eps_1, Spain_pars$s_mu_2, cov_mat_deaths_2$Spain, Spain_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1"), Spain_pars$cps, prior_mu)
mala_US_deaths_gamma5_1 = mala_cp(seed_id=14, "US",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, US_pars$s_mu_1, cov_mat_deaths_1$US, US_pars$eps_1, US_pars$s_mu_2, cov_mat_deaths_2$US, US_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1"), US_pars$cps, prior_mu)
mala_Sweden_deaths_gamma5_1 = mala_cp(seed_id=15, "Sweden",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Sweden_pars$s_mu_1, cov_mat_deaths_1$Sweden, Sweden_pars$eps_1, Sweden_pars$s_mu_2, cov_mat_deaths_2$Sweden, Sweden_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1"), Sweden_pars$cps, prior_mu)
mala_UK_deaths_gamma5_1 = mala_cp(seed_id=16, "UK",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, UK_pars$s_mu_1, cov_mat_deaths_1$UK, UK_pars$eps_1, UK_pars$s_mu_2, cov_mat_deaths_2$UK, UK_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1"), UK_pars$cps, prior_mu)
mala_Brazil_deaths_gamma5_1 = mala_cp(seed_id=17, "Brazil",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Brazil_pars$s_mu_1, cov_mat_deaths_1$Brazil, Brazil_pars$eps_1, Brazil_pars$s_mu_2, cov_mat_deaths_1$Brazil, Brazil_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1"), Brazil_pars$cps, prior_mu)
mala_India_deaths_gamma5_1 = mala_cp(seed_id=18, "India",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, India_pars$s_mu_1, cov_mat_deaths_1$India, India_pars$eps_1, India_pars$s_mu_2, cov_mat_deaths_1$India, India_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1"), India_pars$cps, prior_mu)
mala_China_deaths_gamma5_1 = mala_cp(seed_id=19, "China",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, China_pars$s_mu_1, cov_mat_deaths_1$China, China_pars$eps_1, China_pars$s_mu_2, cov_mat_deaths_2$China, China_pars$eps_2, paste0("mala_deaths_", prior_mu, "5_1"), China_pars$cps, prior_mu)


# Hyper parameters for mu (using uniform prior)
prior_mu = "uniform"

mala_France_deaths_uniform = mala_cp(seed_id=10, "France",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, France_pars$s_mu_1, cov_mat_deaths_1$France, France_pars$eps_1, France_pars$s_mu_2, cov_mat_deaths_2$France, France_pars$eps_2, paste0("mala_deaths_", prior_mu), France_pars$cps, prior_mu)
mala_Germany_deaths_uniform = mala_cp(seed_id=11, "Germany",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Germany_pars$s_mu_1, cov_mat_deaths_1$Germany, Germany_pars$eps_1, Germany_pars$s_mu_2, cov_mat_deaths_2$Germany, Germany_pars$eps_2, paste0("mala_deaths_", prior_mu),Germany_pars$cps, prior_mu)
mala_Italy_deaths_uniform = mala_cp(seed_id=12, "Italy",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Italy_pars$s_mu_1, cov_mat_deaths_1$Italy, Italy_pars$eps_1, Italy_pars$s_mu_2, cov_mat_deaths_2$Italy, Italy_pars$eps_2, paste0("mala_deaths_", prior_mu), Italy_pars$cps, prior_mu)
mala_Spain_deaths_uniform = mala_cp(seed_id=13, "Spain",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Spain_pars$s_mu_1, cov_mat_deaths_1$Spain, Spain_pars$eps_1, Spain_pars$s_mu_2, cov_mat_deaths_2$Spain, Spain_pars$eps_2, paste0("mala_deaths_", prior_mu), Spain_pars$cps, prior_mu)
mala_US_deaths_uniform = mala_cp(seed_id=14, "US",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, US_pars$s_mu_1, cov_mat_deaths_1$US, US_pars$eps_1, US_pars$s_mu_2, cov_mat_deaths_2$US, US_pars$eps_2, paste0("mala_deaths_", prior_mu), US_pars$cps, prior_mu)
mala_Sweden_deaths_uniform = mala_cp(seed_id=15, "Sweden",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Sweden_pars$s_mu_1, cov_mat_deaths_1$Sweden, Sweden_pars$eps_1, Sweden_pars$s_mu_2, cov_mat_deaths_2$Sweden, Sweden_pars$eps_2, paste0("mala_deaths_", prior_mu), Sweden_pars$cps, prior_mu)
mala_UK_deaths_uniform = mala_cp(seed_id=16, "UK",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, UK_pars$s_mu_1, cov_mat_deaths_1$UK, UK_pars$eps_1, UK_pars$s_mu_2, cov_mat_deaths_2$UK, UK_pars$eps_2, paste0("mala_deaths_", prior_mu), UK_pars$cps, prior_mu)
mala_Brazil_deaths_uniform = mala_cp(seed_id=17, "Brazil",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, Brazil_pars$s_mu_1, cov_mat_deaths_1$Brazil, Brazil_pars$eps_1, Brazil_pars$s_mu_2, cov_mat_deaths_1$Brazil, Brazil_pars$eps_2, paste0("mala_deaths_", prior_mu), Brazil_pars$cps, prior_mu)
mala_India_deaths_uniform = mala_cp(seed_id=18, "India",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, India_pars$s_mu_1, cov_mat_deaths_1$India, India_pars$eps_1, India_pars$s_mu_2, cov_mat_deaths_1$India, India_pars$eps_2, paste0("mala_deaths_", prior_mu), India_pars$cps, prior_mu)
mala_China_deaths_uniform = mala_cp(seed_id=19, "China",N, burnin, thin, 1, 1, 0.5, 1, 1, 0.5, China_pars$s_mu_1, cov_mat_deaths_1$China, China_pars$eps_1, China_pars$s_mu_2, cov_mat_deaths_2$China, China_pars$eps_2, paste0("mala_deaths_", prior_mu), China_pars$cps, prior_mu)

