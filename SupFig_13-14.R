# ===================================================================#
# Supplementary figures 13 & 14
# Posterior distributions for estimated parameters
# Viral load errors sensitivity analysis
# ===================================================================#
library(tidyverse)

VL_d <- readRDS("model2_delta_op.rds")
posterior.d <- rstan::extract(VL_d)

beta_all_d <- data.frame(beta = posterior.d$beta_par, run = rep(1,length(posterior.d$beta_par)))
c_all_d <- data.frame(c = posterior.d$c_par, run = rep(1,length(posterior.d$beta_par)))
gamma_all_d <- data.frame(gamma = posterior.d$gamma_par, run = rep(1,length(posterior.d$beta_par)))
theta_all_d <- data.frame(theta = posterior.d$theta, run = rep(1,length(posterior.d$beta_par)))
omega_all_d <- data.frame(omega = posterior.d$omega_par, run = rep(1,length(posterior.d$beta_par)))

VL_d <- readRDS("output_VL_2_d.rds")
posterior.d <- rstan::extract(VL_d)

temp_beta <- data.frame(beta = posterior.d$beta_par, run = rep(2,length(posterior.d$beta_par)))
temp_c <- data.frame(c = posterior.d$c_par, run = rep(2,length(posterior.d$beta_par)))
temp_gamma <- data.frame(gamma = posterior.d$gamma_par, run = rep(2,length(posterior.d$beta_par)))
temp_theta <- data.frame(theta = posterior.d$theta, run = rep(2,length(posterior.d$beta_par)))
temp_omega <- data.frame(omega = posterior.d$omega_par, run = rep(2,length(posterior.d$beta_par)))

beta_all_d <- rbind(beta_all_d, temp_beta)
c_all_d <- rbind(c_all_d, temp_c)
gamma_all_d <- rbind(gamma_all_d, temp_gamma)
theta_all_d <- rbind(theta_all_d, temp_theta)
omega_all_d <- rbind(omega_all_d, temp_omega)

VL_d <- readRDS("output_VL_05_d.rds")
posterior.d <- rstan::extract(VL_d)

temp_beta <- data.frame(beta = posterior.d$beta_par, run = rep(0.5,length(posterior.d$beta_par)))
temp_c <- data.frame(c = posterior.d$c_par, run = rep(0.5,length(posterior.d$beta_par)))
temp_gamma <- data.frame(gamma = posterior.d$gamma_par, run = rep(0.5,length(posterior.d$beta_par)))
temp_theta <- data.frame(theta = posterior.d$theta, run = rep(0.5,length(posterior.d$beta_par)))
temp_omega <- data.frame(omega = posterior.d$omega_par, run = rep(0.5,length(posterior.d$beta_par)))

beta_all_d <- rbind(beta_all_d, temp_beta)
c_all_d <- rbind(c_all_d, temp_c)
gamma_all_d <- rbind(gamma_all_d, temp_gamma)
theta_all_d <- rbind(theta_all_d, temp_theta)
omega_all_d <- rbind(omega_all_d, temp_omega)

theta_all_d_without1 <- theta_all_d[-c(1)] # remove theta[1] because this is fixed at 1.0

betas <- beta_all_d %>%
  gather(key = 'par', value = 'par_value', -run)

cs <- c_all_d %>%
  gather(key = 'par', value = 'par_value', -run)

gammas <- gamma_all_d %>%
  gather(key = 'par', value = 'par_value', -run)

thetas <- theta_all_d_without1 %>%
  gather(key = 'par', value = 'par_value', -run)

omegas <- omega_all_d %>%
  gather(key = 'par', value = 'par_value', -run)

betas$run <- as.factor(betas$run)
cs$run <- as.factor(cs$run)
gammas$run <- as.factor(gammas$run)
gammas$par <- factor(gammas$par, levels = c('gamma.1', 'gamma.2', 'gamma.3', 'gamma.4', 'gamma.5',
                                            'gamma.6', 'gamma.7', 'gamma.8', 'gamma.9', 'gamma.10',
                                            'gamma.11', 'gamma.12', 'gamma.13', 'gamma.14', 'gamma.15',
                                            'gamma.16', 'gamma.17', 'gamma.18', 'gamma.19', 'gamma.20',
                                            'gamma.21', 'gamma.22', 'gamma.23', 'gamma.24', 'gamma.25',
                                            'gamma.26', 'gamma.27', 'gamma.28', 'gamma.29', 'gamma.30',
                                            'gamma.31', 'gamma.32', 'gamma.33', 'gamma.34', 'gamma.35',
                                            'gamma.36', 'gamma.37', 'gamma.38', 'gamma.39', 'gamma.40',
                                            'gamma.41', 'gamma.42', 'gamma.43', 'gamma.44', 'gamma.45',
                                            'gamma.46', 'gamma.47', 'gamma.48', 'gamma.49'))
thetas$run <- as.factor(thetas$run)
thetas$par <- factor(thetas$par, levels = c('theta.2', 'theta.3', 'theta.4', 'theta.5','theta.6'))

omegas$run <- as.factor(omegas$run)
omegas$par <- factor(omegas$par, levels = c('omega.1', 'omega.2', 'omega.3', 'omega.4', 'omega.5',
                                            'omega.6', 'omega.7', 'omega.8', 'omega.9', 'omega.10',
                                            'omega.11', 'omega.12', 'omega.13', 'omega.14', 'omega.15',
                                            'omega.16', 'omega.17', 'omega.18', 'omega.19', 'omega.20',
                                            'omega.21', 'omega.22', 'omega.23', 'omega.24', 'omega.25',
                                            'omega.26', 'omega.27', 'omega.28', 'omega.29', 'omega.30',
                                            'omega.31', 'omega.32', 'omega.33', 'omega.34', 'omega.35',
                                            'omega.36', 'omega.37', 'omega.38', 'omega.39', 'omega.40',
                                            'omega.41', 'omega.42', 'omega.43', 'omega.44', 'omega.45',
                                            'omega.46', 'omega.47', 'omega.48', 'omega.49'))


betas_plot <- ggplot(data=betas) +
  geom_density(aes(x=par_value, fill = run, colour = run), alpha = 0.5) + 
  theme(legend.position = 'none',
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9)) + 
  labs(x = "Value of Beta",
       y = "Density")

cs_plot <- ggplot(data=cs) +
  geom_density(aes(x=par_value, fill = run, colour = run), alpha = 0.5) + 
  theme(legend.position = 'none',
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9)) + 
  labs(x = "Value of c",
       y = "Density")

gammas_plot <- ggplot(data=gammas) +
  geom_density(aes(x=par_value, fill = run, colour = run), alpha = 0.5) +
  facet_wrap(~par) + 
  theme(legend.position = 'none',
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14)) + 
  labs(x = "Value of gamma",
       y = "Density")

thetas_plot <- ggplot(data=thetas) +
  geom_density(aes(x=par_value, fill = run, colour = run), alpha = 0.5) +
  facet_wrap(~par) + xlim(0.1,2.5) + 
  theme(legend.position = 'bottom',
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 10),
        axis.title.x = element_text(size = 11)) + 
  labs(x = "Value of theta",
       y = "Density")

omegas_plot <- ggplot(data=omegas) +
  geom_density(aes(x=log10(par_value), fill = run, colour = run), alpha = 0.5) +
  facet_wrap(~par) + 
  theme(legend.position = 'none',
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14)) + 
  labs(x = "Log10(value of omega)",
       y = "Density")

ggsave("betas_cs_d.png", 
       plot = wrap_elements(betas_plot | cs_plot ), 
       width = 5500, height = 2000, units = "px", dpi = 700)
ggsave("thetas_d.png",
       thetas_plot,
       width = 5000, height = 2500, units = "px", dpi = 700)
ggsave("gammas_d.png",
       gammas_plot,
       width = 15000, height = 9000, units = "px", dpi = 700)
ggsave("omegas_d.png",
       omegas_plot,
       width = 15000, height = 9000, units = "px", dpi = 700)

############### Omicron
VL_o <- readRDS("model2_omicron_op.rds")
posterior.o <- rstan::extract(VL_o)

beta_all_o <- data.frame(beta = posterior.o$beta_par, run = rep(1,length(posterior.o$beta_par)))
c_all_o <- data.frame(c = posterior.o$c_par, run = rep(1,length(posterior.o$beta_par)))
gamma_all_o <- data.frame(gamma = posterior.o$gamma_par, run = rep(1,length(posterior.o$beta_par)))
theta_all_o <- data.frame(theta = posterior.o$theta, run = rep(1,length(posterior.o$beta_par)))
omega_all_o <- data.frame(omega = posterior.o$omega_par, run = rep(1,length(posterior.o$beta_par)))

VL_o <- readRDS("output_VL_2_o.rds")
posterior.o <- rstan::extract(VL_o)

temp_beta <- data.frame(beta = posterior.o$beta_par, run = rep(2,length(posterior.o$beta_par)))
temp_c <- data.frame(c = posterior.o$c_par, run = rep(2,length(posterior.o$beta_par)))
temp_gamma <- data.frame(gamma = posterior.o$gamma_par, run = rep(2,length(posterior.o$beta_par)))
temp_theta <- data.frame(theta = posterior.o$theta, run = rep(2,length(posterior.o$beta_par)))
temp_omega <- data.frame(omega = posterior.o$omega_par, run = rep(2,length(posterior.o$beta_par)))

beta_all_o <- rbind(beta_all_o, temp_beta)
c_all_o <- rbind(c_all_o, temp_c)
gamma_all_o <- rbind(gamma_all_o, temp_gamma)
theta_all_o <- rbind(theta_all_o, temp_theta)
omega_all_o <- rbind(omega_all_o, temp_omega)

VL_o <- readRDS("output_VL_05_o.rds")
posterior.o <- rstan::extract(VL_o)

temp_beta <- data.frame(beta = posterior.o$beta_par, run = rep(0.5,length(posterior.o$beta_par)))
temp_c <- data.frame(c = posterior.o$c_par, run = rep(0.5,length(posterior.o$beta_par)))
temp_gamma <- data.frame(gamma = posterior.o$gamma_par, run = rep(0.5,length(posterior.o$beta_par)))
temp_theta <- data.frame(theta = posterior.o$theta, run = rep(0.5,length(posterior.o$beta_par)))
temp_omega <- data.frame(omega = posterior.o$omega_par, run = rep(0.5,length(posterior.o$beta_par)))

beta_all_o <- rbind(beta_all_o, temp_beta)
c_all_o <- rbind(c_all_o, temp_c)
gamma_all_o <- rbind(gamma_all_o, temp_gamma)
theta_all_o <- rbind(theta_all_o, temp_theta)
omega_all_o <- rbind(omega_all_o, temp_omega)

theta_all_o_without1 <- theta_all_o[-c(1)]

betas <- beta_all_o %>%
  gather(key = 'par', value = 'par_value', -run)

cs <- c_all_o %>%
  gather(key = 'par', value = 'par_value', -run)

gammas <- gamma_all_o %>%
  gather(key = 'par', value = 'par_value', -run)

thetas <- theta_all_o_without1 %>%
  gather(key = 'par', value = 'par_value', -run)

omegas <- omega_all_o %>%
  gather(key = 'par', value = 'par_value', -run)

betas$run <- as.factor(betas$run)
cs$run <- as.factor(cs$run)
gammas$run <- as.factor(gammas$run)
gammas$par <- factor(gammas$par, levels = c('gamma.1', 'gamma.2', 'gamma.3', 'gamma.4', 'gamma.5',
                                            'gamma.6', 'gamma.7', 'gamma.8', 'gamma.9', 'gamma.10',
                                            'gamma.11', 'gamma.12', 'gamma.13', 'gamma.14', 'gamma.15',
                                            'gamma.16', 'gamma.17', 'gamma.18', 'gamma.19', 'gamma.20',
                                            'gamma.21', 'gamma.22', 'gamma.23', 'gamma.24', 'gamma.25',
                                            'gamma.26', 'gamma.27', 'gamma.28', 'gamma.29', 'gamma.30',
                                            'gamma.31', 'gamma.32', 'gamma.33', 'gamma.34', 'gamma.35',
                                            'gamma.36', 'gamma.37', 'gamma.38', 'gamma.39', 'gamma.40',
                                            'gamma.41', 'gamma.42', 'gamma.43', 'gamma.44', 'gamma.45',
                                            'gamma.46', 'gamma.47', 'gamma.48', 'gamma.49'))
thetas$run <- as.factor(thetas$run)
thetas$par <- factor(thetas$par, levels = c('theta.2', 'theta.3', 'theta.4', 'theta.5','theta.6'))

omegas$run <- as.factor(omegas$run)
omegas$par <- factor(omegas$par, levels = c('omega.1', 'omega.2', 'omega.3', 'omega.4', 'omega.5',
                                            'omega.6', 'omega.7', 'omega.8', 'omega.9', 'omega.10',
                                            'omega.11', 'omega.12', 'omega.13', 'omega.14', 'omega.15',
                                            'omega.16', 'omega.17', 'omega.18', 'omega.19', 'omega.20',
                                            'omega.21', 'omega.22', 'omega.23', 'omega.24', 'omega.25',
                                            'omega.26', 'omega.27', 'omega.28', 'omega.29', 'omega.30',
                                            'omega.31', 'omega.32', 'omega.33', 'omega.34', 'omega.35',
                                            'omega.36', 'omega.37', 'omega.38', 'omega.39', 'omega.40',
                                            'omega.41', 'omega.42', 'omega.43', 'omega.44', 'omega.45',
                                            'omega.46', 'omega.47', 'omega.48', 'omega.49'))


betas_plot <- ggplot(data=betas) +
  geom_density(aes(x=par_value, fill = run, colour = run), alpha = 0.5) + 
  theme(legend.position = 'none',
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9)) + 
  labs(x = "Value of Beta",
       y = "Density")

cs_plot <- ggplot(data=cs) +
  geom_density(aes(x=par_value, fill = run, colour = run), alpha = 0.5) + 
  theme(legend.position = 'none',
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9)) + 
  labs(x = "Value of c",
       y = "Density")

gammas_plot <- ggplot(data=gammas) +
  geom_density(aes(x=par_value, fill = run, colour = run), alpha = 0.5) +
  facet_wrap(~par) + 
  theme(legend.position = 'none',
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14)) + 
  labs(x = "Value of gamma",
       y = "Density")

thetas_plot <- ggplot(data=thetas) +
  geom_density(aes(x=par_value, fill = run, colour = run), alpha = 0.5) +
  facet_wrap(~par) + xlim(0.3,3.2) + 
  theme(legend.position = 'none',
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 10),
        axis.title.x = element_text(size = 11)) + 
  labs(x = "Value of theta",
       y = "Density")

omegas_plot <- ggplot(data=omegas) +
  geom_density(aes(x=log10(par_value), fill = run, colour = run), alpha = 0.5) +
  facet_wrap(~par) + 
  theme(legend.position = 'none',
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14)) + 
  labs(x = "Log10(value of omega)",
       y = "Density")

ggsave("betas_cs_o.png", 
       plot = wrap_elements(betas_plot | cs_plot ), 
       width = 5500, height = 2000, units = "px", dpi = 700)
ggsave("thetas_o.png",
       thetas_plot,
       width = 5000, height = 2500, units = "px", dpi = 700)
ggsave("gammas_o.png",
       gammas_plot,
       width = 15000, height = 9000, units = "px", dpi = 700)
ggsave("omegas_o.png",
       omegas_plot,
       width = 15000, height = 9000, units = "px", dpi = 700)
