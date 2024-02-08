# ===================================================================#
# Supplementary figures 10
# Trace plots & R-hat plots for virus neutralisation model
# ===================================================================#
library(tidyverse)
library(gridExtra)
library(rstan)
library(readxl)
library(loo)
library(patchwork)
library(bayesplot)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

model2_delta_op <- readRDS("op_SA_virus_d.rds")
model2_omicron_op <- readRDS("op_SA_virus_o.rds")

post.delta.2 <- rstan::extract(model2_delta_op)
post.omicron.2 <-rstan::extract(model2_omicron_op)
params = c("beta_par", "delta_par", "gamma_par", "omega_mod")

print(model2_delta_op, pars = params, digits = 12)
rstan::check_hmc_diagnostics(model2_delta_op)
print(model2_omicron_op, pars = params, digits = 12)
rstan::check_hmc_diagnostics(model2_omicron_op)

trace.m2.d <- traceplot(model2_delta_op, pars = params)
rhat.m2.d <- mcmc_rhat(rhat(model2_delta_op))
trace.m2.o <- traceplot(model2_omicron_op, pars = params)
rhat.m2.o <- mcmc_rhat(rhat(model2_omicron_op))

ggsave("SuppFig_10_trace_delta.png", 
       plot = (trace.m2.d), 
       width = 18000, height = 7000, units = "px", dpi = 700)
ggsave("SuppFig_10_trace_omicron.png", 
       plot = (trace.m2.o), 
       width = 18000, height = 7000, units = "px", dpi = 700)
ggsave("SuppFig_10_Rhat.png", 
       plot = (rhat.m2.d | rhat.m2.o) + plot_layout(guides = 'collect'), 
       width = 15000, height = 6000, units = "px", dpi = 700)

