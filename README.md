# Description

This repository provides the code accompanying the paper,

Browning R, Sulem D, Mengersen K, Rivoirard V, Rousseau J (2021) Simple discrete-time self-exciting models can describe complex dynamic processes: A case study of COVID-19. PLOS ONE 16(4): e0250015. https://doi.org/10.1371/journal.pone.0250015

# Data

The file 'time_series_covid19_deaths_global.csv' contains the COVID-19 mortality data used for these analysis. The data were obtained from the COVID-19 Data Repository by the Center for Systems Science and Engineering (CSSE) at Johns Hopkins University.

# Code

The R scripts of the form 'hp_covid19_deaths_XX.R' contain the code needed to run the model (for all countries and under multiple prior settings) for multiple types of analyses:

- hp_covid19_deaths.R: the main results from the paper
- hp_covid19_deaths_interpol.R: run the model, removing some observations and interpolate the missing data
- hp_covid19_deaths_lfocv.R: perform leave-future-out cross validation
- hp_covid19_deaths_training.R: perform out of sample validation at multiple time points throughout the time series

Each of these scripts sources the functions required to load and prepare the data, load the MCMC algorithm, and any other bespoke functions required for the particular analysis being performed.

Algorithmic parameters are also read into each of these scripts. For each country, tuning parameters and parameters for data preparation are saved in 'mcmc_pars.RData'. The file 'cov_mat_deaths.RData' contains covariance matrices required for the MALA algorithm, that were obtained from previous MCMC runs using single Metropolis-Hastings updates.

These details relate to the master branch of this repository. However, a similar structure is followed in the branch 'subsequent_phases', that extends the analysis by adding in more recent data (up to 2021).
