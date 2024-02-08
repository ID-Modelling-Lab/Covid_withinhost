# ===================================================================#
# model 4
## V: beta, c
## VG: -
## AG: gamma, omega
# ===================================================================#
library(FME)
library(tidyverse)
library(rstan)
library(readxl)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ===================================================================#
# delta model fit
# ===================================================================#

df.delta <- as.data.frame(read_excel('VL_standata.xlsx'))

## arrange patients by vacc group and age group
df.delta$AnonID <- as.factor(df.delta$AnonID)
df.delta <- df.delta %>% 
  dplyr::mutate_at(c('vaccgroup1', 'agegroup', 'ageiden', 'vaccgroup_range'), as.numeric)
df.delta <- arrange(df.delta, vaccgroup1, agegroup, ageiden, agegroup_range)

# no. of samples in each age group & time since vaccination combination (ag_tsv)
# position of the 1st value of each ag_tsv group
ct_samplesinagtsv <- df.delta %>% 
  dplyr::group_by(vaccgroup1, agegroup, vaccgroup_range, ageiden, agegroup_range) %>% 
  dplyr::count()
pos_samplesinagtsv <- rep(1, nrow(ct_samplesinagtsv))
temp = cumsum(ct_samplesinagtsv$n)
pos_samplesinagtsv[2:(nrow(ct_samplesinagtsv))] <- pos_samplesinagtsv[2:(nrow(ct_samplesinagtsv))] + temp[1:(nrow(ct_samplesinagtsv) - 1)]

# vector of timepoints of interest for stan
t_solve <- sort(as.numeric(unique(as.character(df.delta$DayIP))))

init_func <- function(...) {
  list(beta_par = 3*10^-11,
       gamma_par = array(rep(5, times = length(unique(df.delta$vaccgroup_range)))),
       rho = 1,
       eta = 1,
       c_par = 1,
       omega_par = array(rep(1*10^-4, times = length(unique(df.delta$vaccgroup_range)))),
       theta_mult_g = array(rep(1.0, times = 5)),
       theta_mult_o = array(rep(1.0, times = 5)),
       phi =  2,
       lambda = 5 )}

dfst <- list(L = dim(df.delta)[1],
             VGtsv = length(unique(df.delta$vaccgroup_range)),
             AGtsv = length(unique(df.delta$agegroup_range)),
             l_tsolve = length(t_solve),
             logvirus = df.delta$logvalue,
             t0 = 0,
             ts = as.numeric(df.delta$DayIP),
             ts_ind = as.integer(df.delta$dayip_index),
             t_solve = t_solve,
             pos_VGtsvofAGtsv = ct_samplesinagtsv$vaccgroup_range,
             pos_AGofAGtsv = ct_samplesinagtsv$ageiden,
             pos_agtsv = pos_samplesinagtsv,
             count_agtsv = ct_samplesinagtsv$n,
             rel_tol = 10^-8,
             abs_tol = 10^-6,
             max_num_steps = 10^7,
             which_fixed = 1,
             fixed_value = 1.0)

param_to_save = c("beta_par", "c_par", "theta_g", "theta_o",
                  "gamma_mod", "omega_mod","log_lik", "sumloglike")

rstan_options(disable_march_warning = TRUE)
test = stan("model4.stan",
            data = dfst, 
            chains = 4,
            iter = 3000, 
            warmup = 1500,
            pars = param_to_save,
            include = TRUE,
            save_warmup = FALSE,
            cores = parallel::detectCores(), 
            control = list(adapt_delta = 0.82), 
            init = init_func)

# save test
saveRDS(test, "model4_delta_op.rds")

# ===================================================================#
# omicron model fit
# ===================================================================#

df.omicron <- as.data.frame(read_excel('VL_standata.xlsx'))

## arrange patients by vacc group and age group
df.omicron$AnonID <- as.factor(df.omicron$AnonID)
df.omicron <- df.omicron %>% 
  dplyr::mutate_at(c('vaccgroup1', 'agegroup', 'ageiden', 'vaccgroup_range'), as.numeric)
df.omicron <- arrange(df.omicron, vaccgroup1, agegroup, ageiden, agegroup_range)

# no. of samples in each age group & time since vaccination combination (ag_tsv)
# position of the 1st value of each ag_tsv group
ct_samplesinagtsv <- df.omicron %>% 
  dplyr::group_by(vaccgroup1, agegroup, vaccgroup_range, ageiden, agegroup_range) %>% 
  dplyr::count()
pos_samplesinagtsv <- rep(1, nrow(ct_samplesinagtsv))
temp = cumsum(ct_samplesinagtsv$n)
pos_samplesinagtsv[2:(nrow(ct_samplesinagtsv))] <- pos_samplesinagtsv[2:(nrow(ct_samplesinagtsv))] + temp[1:(nrow(ct_samplesinagtsv) - 1)]

# vector of timepoints of interest for stan
t_solve <- sort(as.numeric(unique(as.character(df.omicron$DayIP))))

init_func <- function(...) {
  list(beta_par = 5*10^-11,
       gamma_par = array(rep(5, times = length(unique(df.omicron$vaccgroup_range)))),
       rho = 1,
       eta = 1,
       c_par = 1,
       omega_par = array(rep(1*10^-4, times = length(unique(df.omicron$vaccgroup_range)))),
       theta_mult_g = array(rep(1.0, times = 5)),
       theta_mult_o = array(rep(1.0, times = 5)),
       phi =  2,
       lambda = 5 )}

dfst <- list(L = dim(df.omicron)[1],
             VGtsv = length(unique(df.omicron$vaccgroup_range)),
             AGtsv = length(unique(df.omicron$agegroup_range)),
             l_tsolve = length(t_solve),
             logvirus = df.omicron$logvalue,
             t0 = 0,
             ts = as.numeric(df.omicron$DayIP),
             ts_ind = as.integer(df.omicron$dayip_index),
             t_solve = t_solve,
             pos_VGtsvofAGtsv = ct_samplesinagtsv$vaccgroup_range,
             pos_AGofAGtsv = ct_samplesinagtsv$ageiden,
             pos_agtsv = pos_samplesinagtsv,
             count_agtsv = ct_samplesinagtsv$n,
             rel_tol = 10^-8,
             abs_tol = 10^-6,
             max_num_steps = 10^7,
             which_fixed = 1,
             fixed_value = 1.0)

param_to_save = c("beta_par", "c_par", "theta_g", "theta_o",
                  "gamma_mod", "omega_mod","log_lik", "sumloglike")

rstan_options(disable_march_warning = TRUE)
test = stan("model4.stan",
            data = dfst, 
            chains = 4,
            iter = 3000, 
            warmup = 1500,
            pars = param_to_save,
            include = TRUE,
            save_warmup = FALSE,
            cores = parallel::detectCores(), 
            control = list(adapt_delta = 0.82), 
            init = init_func)

# save test
saveRDS(test, "model4_omicron_op.rds")
