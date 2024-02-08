# ===================================================================#
# Supplementary figures 3 - 5
# Trace plots, R-hat plots, posterior distributions
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

# =================================================== #
# metrics for model selection
# =================================================== #

#LL: loglikelihood
#k:the number of estimated parameters in the model
#n: the number of observations used in the model
cal.AIC.c <- function(LL,k,n){
  AIC.c = -2*LL +2*k + 2*k*(k+1)/(n - k -1)
  print(AIC.c)
}
cal.AIC <- function(LL, k){
  AIC = -2*LL + 2*k
  print(AIC)
}
cal.BIC <- function(LL,k,n){
  BIC = -2*LL+k*log(n)
  print(BIC)
}

# =================================================== #
# Supplementary Figure 3-1(a) and 3-1(b)
# Trace plot and R-hat plot for model type 1
# =================================================== #
###### model 1, gamma vg, omega vg
post.delta.1 <- rstan::extract(model1_delta_op)
post.omicron.1 <-rstan::extract(model1_omicron_op)
params = c("beta_par", "c_par", "gamma_par", "omega_par")

print(model1_delta_op, pars = params, digits = 12)
rstan::check_hmc_diagnostics(model1_delta_op)
print(model1_omicron_op, pars = params, digits = 12)
rstan::check_hmc_diagnostics(model1_omicron_op)

trace.m1.d <- traceplot(model1_delta_op, pars = params)
rhat.m1.d <- mcmc_rhat(rhat(model1_delta_op))
trace.m1.o <- traceplot(model1_omicron_op, pars = params)
rhat.m1.o <- mcmc_rhat(rhat(model1_omicron_op))

ggsave("SuppFig_3_1_trace.png", 
       plot = (trace.m1.d | trace.m1.o) + plot_layout(guides = 'collect'), 
       width = 20000, height = 6500, units = "px", dpi = 800)
ggsave("SuppFig_3_1_Rhat.png", 
       plot = (rhat.m1.d | rhat.m1.o) + plot_layout(guides = 'collect'), 
       width = 15000, height = 6000, units = "px", dpi = 800)

# =================================================== #
# Supplementary Figure 3-2(a) and 3-2(b)
# Trace plot and R-hat plot for model type 2
# =================================================== #
###### model 2, gamma vg, omega ag
post.delta.2 <- rstan::extract(model2_delta_op)
post.omicron.2 <-rstan::extract(model2_omicron_op)
params = c("beta_par", "c_par", "gamma_par", "omega_mod")

print(model2_delta_op, pars = params, digits = 12)
rstan::check_hmc_diagnostics(model2_delta_op)
print(model2_omicron_op, pars = params, digits = 12)
rstan::check_hmc_diagnostics(model2_omicron_op)

trace.m2.d <- traceplot(model2_delta_op, pars = params)
rhat.m2.d <- mcmc_rhat(rhat(model2_delta_op))
trace.m2.o <- traceplot(model2_omicron_op, pars = params)
rhat.m2.o <- mcmc_rhat(rhat(model2_omicron_op))

ggsave("SuppFig_3_2_trace_delta.png", 
       plot = (trace.m2.d), 
       width = 18000, height = 7000, units = "px", dpi = 700)
ggsave("SuppFig_3_2_trace_omicron.png", 
       plot = (trace.m2.o), 
       width = 18000, height = 7000, units = "px", dpi = 700)
ggsave("SuppFig_3_2_Rhat.png", 
       plot = (rhat.m2.d | rhat.m2.o) + plot_layout(guides = 'collect'), 
       width = 15000, height = 6000, units = "px", dpi = 700)


# =================================================== #
# Supplementary Figure 3-3(a) and 3-3(b)
# Trace plot and R-hat plot for model type 3
# =================================================== #
###### model 3, gamma ag, omega vg
post.delta.3 <- rstan::extract(model3_delta_op)
post.omicron.3 <-rstan::extract(model3_omicron_op)
params = c("beta_par", "c_par", "gamma_mod", "omega_par")

print(model3_delta_op, pars = params, digits = 12)
rstan::check_hmc_diagnostics(model3_delta_op)
print(model3_omicron_op, pars = params, digits = 12)
rstan::check_hmc_diagnostics(model3_omicron_op)

trace.m3.d <- traceplot(model3_delta_op, pars = params)
rhat.m3.d <- mcmc_rhat(rhat(model3_delta_op))
trace.m3.o <- traceplot(model3_omicron_op, pars = params)
rhat.m3.o <- mcmc_rhat(rhat(model3_omicron_op))

ggsave("SuppFig_3_3_trace_delta.png", 
       plot = trace.m3.d, 
       width = 18000, height = 7000, units = "px", dpi = 700)
ggsave("SuppFig_3_3_trace_omicron.png", 
       plot = trace.m3.o, 
       width = 18000, height = 7000, units = "px", dpi = 700)
ggsave("SuppFig_3_3_Rhat.png", 
       plot = (rhat.m3.d | rhat.m3.o) + plot_layout(guides = 'collect'), 
       width = 15000, height = 6000, units = "px", dpi = 800)


# =================================================== #
# Supplementary Figure 3-4(a) and 3-4(b)
# Trace plot and R-hat plot for model type 4
# =================================================== #
###### model 4, gamma ag, omega vg
post.delta.4 <- rstan::extract(model4_delta_op)
post.omicron.4 <-rstan::extract(model4_omicron_op)
params = c("beta_par", "c_par", "gamma_mod", "omega_mod")

print(model4_delta_op, pars = params, digits = 12)
rstan::check_hmc_diagnostics(model4_delta_op)
print(model4_omicron_op, pars = params, digits = 12)
rstan::check_hmc_diagnostics(model4_omicron_op)

trace.m4.d <- traceplot(model4_delta_op, pars = params)
rhat.m4.d <- mcmc_rhat(rhat(model4_delta_op))
trace.m4.o <- traceplot(model4_omicron_op, pars = params)
rhat.m4.o <- mcmc_rhat(rhat(model4_omicron_op))

ggsave("SuppFig_3_4_trace_delta.png", 
       plot = trace.m4.d, 
       width = 18000, height = 7000, units = "px", dpi = 700)
ggsave("SuppFig_3_4_trace_omicron.png", 
       plot = trace.m4.o, 
       width = 18000, height = 7000, units = "px", dpi = 700)
ggsave("SuppFig_3_4_Rhat.png", 
       plot = (rhat.m4.d | rhat.m4.o) + plot_layout(guides = 'collect'), 
       width = 15000, height = 6000, units = "px", dpi = 800)


# =================================================== #
# Supplementary Figure 4
# Posterior distributions of age modifier values for 
# gamma parameter of model type 4
# =================================================== #
theta.d <- mcmc_areas(as.matrix(model4_delta_op),
                      pars = c("theta_g[2]", "theta_g[3]", "theta_g[4]", "theta_g[5]", "theta_g[6]"),
                      prob = 0.95)
theta.o <- mcmc_areas(as.matrix(model4_omicron_op),
                      pars = c("theta_g[2]", "theta_g[3]", "theta_g[4]", "theta_g[5]", "theta_g[6]"),
                      prob = 0.95)

ggsave("SuppFig_4.png",
       (theta.d + xlim(0,4)) / (theta.o + xlim(0,4)),
       width = 4000, height = 6500, units = "px", dpi = 800)

# =================================================== #
# Supplementary Figure 5
# Posterior distributions of age modifier values for 
# omega parameter of model type 2
# =================================================== #
######## theta age modifier
theta.d <- mcmc_areas(as.matrix(model2_delta_op),
                      pars = c("theta[2]", "theta[3]", "theta[4]", "theta[5]", "theta[6]"),
                      prob = 0.95)

theta.o <- mcmc_areas(as.matrix(model2_omicron_op),
                      pars = c("theta[2]", "theta[3]", "theta[4]", "theta[5]", "theta[6]"),
                      prob = 0.95)

ggsave("SuppFig_5.png",
       (theta.d + xlim(0,4)) / (theta.o + xlim(0,4)),
       width = 4000, height = 6500, units = "px", dpi = 800)

# =================================================== #
# Table 2
# Summary of goodness fit
# =================================================== #
########### model comparison ###########
## model 1, gamma vg, omega vg
# delta: k = 100, omicron: k = 100

AIC <- cal.AIC( LL= median(post.delta.1$sumloglike), k = 100)
AICc <- cal.AIC.c( LL= median(post.delta.1$sumloglike), k = 100, n = dim(df.delta.ss)[1])
BIC <- cal.BIC( LL= median(post.delta.1$sumloglike), k = 100, n = dim(df.delta.ss)[1])

medianLL <- median(post.delta.1$sumloglike)
medianLL

AIC <- cal.AIC( LL= median(post.omicron.1$sumloglike), k = 100)
AICc <- cal.AIC.c( LL= median(post.omicron.1$sumloglike), k = 100, n = dim(df.omicron.ss)[1])
BIC <- cal.BIC( LL= median(post.omicron.1$sumloglike), k = 100, n = dim(df.omicron.ss)[1])

medianLL <- median(post.omicron.1$sumloglike)
medianLL

## model 2, gamma vg, omega ag (age modification adds 5 parameters)
# delta: k = 105, omicron: k = 105
AIC <- cal.AIC( LL= median(post.delta.2$sumloglike), k = 105)
AICc <- cal.AIC.c( LL= median(post.delta.2$sumloglike), k = 105, n = dim(df.delta.ss)[1])
BIC <- cal.BIC( LL= median(post.delta.2$sumloglike), k = 105, n = dim(df.delta.ss)[1])

medianLL <- median(post.delta.2$sumloglike) 
medianLL

AIC <- cal.AIC( LL= median(post.omicron.2$sumloglike), k = 105)
AICc <- cal.AIC.c( LL= median(post.omicron.2$sumloglike), k = 105, n = dim(df.omicron.ss)[1])
BIC <- cal.BIC( LL= median(post.omicron.2$sumloglike), k = 105, n = dim(df.omicron.ss)[1])

medianLL <- median(post.omicron.2$sumloglike) 
medianLL

## model 3, gamma ag, omega vg (age modification adds 5 pars)
# delta: k = 107, omicron: k = 105

AIC <- cal.AIC( LL= median(post.delta.3$sumloglike), k = 105)
AICc <- cal.AIC.c( LL= median(post.delta.3$sumloglike), k = 105, n = dim(df.delta.ss)[1])
BIC <- cal.BIC( LL= median(post.delta.3$sumloglike), k = 105, n = dim(df.delta.ss)[1])

medianLL <- median(post.delta.3$sumloglike) 
medianLL

AIC <- cal.AIC( LL= median(post.omicron.3$sumloglike), k = 105)
AICc <- cal.AIC.c( LL= median(post.omicron.3$sumloglike), k = 105, n = dim(df.omicron.ss)[1])
BIC <- cal.BIC( LL= median(post.omicron.3$sumloglike), k = 105, n = dim(df.omicron.ss)[1])

medianLL <- median(post.omicron.3$sumloglike) 
medianLL

## model 4, gamma ag, omega ag (age modification adds 5 pars)
# delta: k = 112, omicron: k = 110

AIC <- cal.AIC( LL= median(post.delta.4$sumloglike), k = 110)
AICc <- cal.AIC.c( LL= median(post.delta.4$sumloglike), k = 110, n = dim(df.delta.ss)[1])
BIC <- cal.BIC( LL= median(post.delta.4$sumloglike), k = 110, n = dim(df.delta.ss)[1])

medianLL <- median(post.delta.4$sumloglike)
medianLL

AIC <- cal.AIC( LL= median(post.omicron.4$sumloglike), k = 110)
AICc <- cal.AIC.c( LL= median(post.omicron.4$sumloglike), k = 110, n = dim(df.omicron.ss)[1])
BIC <- cal.BIC( LL= median(post.omicron.4$sumloglike), k = 110, n = dim(df.omicron.ss)[1])

medianLL <- median(post.omicron.4$sumloglike)
medianLL

