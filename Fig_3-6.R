# ===================================================================#
# Figures 3 - 6 of main manuscript
# ===================================================================#
library(FME)
library(tidyverse)
library(gridExtra)
library(rstan)
library(readxl)
library(patchwork)
library(ggh4x)
library(bayesplot)

# load rds files (model fits)
model2_delta_op <- readRDS("model2_delta_op.rds")
model2_omicron_op <- readRDS("model2_omicron_op.rds")

# load rdata file for vaccination regime and age group labelling
load("labels.RData")

post.delta.2 <- rstan::extract(model2_delta_op)
post.omicron.2 <-rstan::extract(model2_omicron_op)
params = c("beta_par", "c_par", "gamma_par", "omega_mod")

# ===================================================================#
# Figure 3A - 3D
# Plots for Gamma parameter (infected cell clearance by immunity)
# ===================================================================#
# delta
model2_d_gamma <- data.frame(gamma = post.delta.2$gamma_par)
colnames(model2_d_gamma) <- vgtsv_label_d

model2_d_gamma_g <- model2_d_gamma %>%
  gather(key = 'vg', value = 'par_val')

model2_d_gamma_g[c('vg_label', 'tsv_label')] <- stringr::str_split_fixed(as.character(model2_d_gamma_g$vg), ', ', 2)
model2_d_gamma_g$vg_label <- factor(model2_d_gamma_g$vg_label, 
                                    levels = c("Unvaccinated", "Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3"))
model2_d_gamma_g$tsv_label <- factor(model2_d_gamma_g$tsv_label, 
                                     levels = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"))

# omicron
model2_o_gamma <- data.frame(gamma = post.omicron.2$gamma_par)
colnames(model2_o_gamma) <- vgtsv_label_o

model2_o_gamma_g <- model2_o_gamma %>%
  gather(key = 'vg', value = 'par_val')

model2_o_gamma_g[c('vg_label', 'tsv_label')] <- str_split_fixed(as.character(model2_o_gamma_g$vg), ', ', 2)
model2_o_gamma_g$vg_label <- factor(model2_o_gamma_g$vg_label, 
                                    levels = c("Unvaccinated", "Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3"))
model2_o_gamma_g$tsv_label <- factor(model2_o_gamma_g$tsv_label, 
                                     levels = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"))

# ===========================================#
# Figure 3A
# Distribution for gamma parameter
# Delta dataset
# ===========================================#
model2.d.g.vg.plot <- ggplot(model2_d_gamma_g, aes(x=tsv_label, y = par_val, fill = vg_label)) +
  geom_violin(width = 0.75)+
  theme(plot.title = element_text(hjust=0.5,face="bold",size=14,colour = "blue"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        legend.position = "none",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5,size=12))+
  labs(fill = "Vaccine group", 
       x = "Days since last vaccine dose",
       y = "Difference between infected cell clearance by immune response and unvaccinated group's median value, per immune cell per day",
       title = "Delta VOC infections") +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange",
               color="black",
               alpha = 0.5,
               size = 0.15,
               position = position_dodge(0.75)) + 
  ylim(1,10)

# ===========================================#
# Figure 3B
# Median values for each vaccination regime
# Delta dataset
# ===========================================#
model2_d_gamma_med <- model2_d_gamma_g %>% group_by(vg, vg_label, tsv_label) %>% summarise(median = median(par_val))
model2_d_gamma_med$vg_label <- factor(model2_d_gamma_med$vg_label, levels = c("Unvaccinated", "Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3"))
model2_d_gamma_med$tsv_label <- factor(model2_d_gamma_med$tsv_label, levels = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"))

model2.d.g.vg.med <- ggplot(model2_d_gamma_med, aes(x=tsv_label, y=median, group = vg_label)) +
  geom_point(aes(shape = vg_label, colour = vg_label), size = 3, stroke = 2) +
  geom_line(aes(colour = vg_label), alpha = 0.4, linewidth = 2) +
  theme(plot.title = element_text(hjust=0.5,face="bold",size=14,colour = "blue"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        legend.position = "none",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 9.5))+
  labs(shape = "Vaccine group",
       colour = "Vaccine group") +
  labs(x = "Days since last vaccine dose",
       y = "Infected cell clearance by immune response and unvaccinated group's median value, per immune cell per day",
       title = "Delta VOC infections") +
  ylim(1,10)

# ===========================================#
# Figure 3C
# Distribution for gamma parameter
# Omicron dataset
# ===========================================#
model2.o.g.vg.plot <- ggplot(model2_o_gamma_g, aes(x=tsv_label, y = par_val, fill = vg_label)) +
  geom_violin(width = 0.75)+
  theme(plot.title = element_text(hjust=0.5,face="bold",size=14,colour = "blue"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 12))+
  labs(fill = "Vaccine group", 
       x = "Days since last vaccine dose",
       y = "Difference between infected cell clearance by immune response and unvaccinated group's median value, per immune cell per day",
       title = "Omicron VOC infections") +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange",
               color="black",
               alpha = 0.5,
               size = 0.15,
               position = position_dodge(0.75)) + 
  ylim(1,10)

# ===========================================#
# Figure 3D
# Median values for each vaccination regime
# Omicron dataset
# ===========================================#
model2_o_gamma_med <- model2_o_gamma_g %>% group_by(vg, vg_label, tsv_label) %>% summarise(median = median(par_val))
model2_o_gamma_med$vg_label <- factor(model2_o_gamma_med$vg_label, levels = c("Unvaccinated", "Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3"))
model2_o_gamma_med$tsv_label <- factor(model2_o_gamma_med$tsv_label, levels = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"))

model2.o.g.vg.med <- ggplot(model2_o_gamma_med, aes(x=tsv_label, y=median, group = vg_label)) +
  geom_point(aes(shape = vg_label, colour = vg_label), size = 3, stroke = 2) +
  geom_line(aes(colour = vg_label), alpha = 0.4, linewidth = 2) +
  theme(plot.title = element_text(hjust=0.5,face="bold",size=14,colour = "blue"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 9.5))+
  labs(shape = "Vaccine group",
       colour = "Vaccine group") +
  labs(x = "Days since last vaccine dose",
       y = "Infected cell clearance by immune response and unvaccinated group's median value, per immune cell per day",
       title = "Omicron VOC infections") + 
  ylim(1,10)

# ===========================================#
# Compiled Figure 3
# ===========================================#
byvg <- ((model2.d.g.vg.plot + model2.d.g.vg.med + model2.o.g.vg.plot + model2.o.g.vg.med) + plot_layout(widths = c(2,1))) + plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 10)) & 
  scale_fill_discrete(limits = c("Unvaccinated", "Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3")) & 
  scale_shape_manual(values=1:nlevels(model2_o_gamma_med$vg_label)) & 
  scale_color_discrete(limits = c("Unvaccinated", "Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3"))

ggsave("gamma_fig3_paper_scale_v2.png", 
       plot = wrap_elements(byvg) + 
         labs(tag = "Infected cell clearance by immune response, per immune cell per day") +
         theme(
           plot.tag = element_text(size = rel(1), angle = 90, face = "bold"),
           plot.tag.position = "left"
         ), 
       width = 15000, height = 7000, units = "px", dpi = 800)

# ===================================================================#
# Figure 4A
# Plots for Omega parameter (growth rate of immunity)
# Delta dataset
# ===================================================================#
# delta
model2_d_omega <- data.frame(gamma = post.delta.2$omega_mod)
colnames(model2_d_omega) <- ag_label_d

model2_d_omega_g <- model2_d_omega %>%
  gather(key = 'agtsv', value = 'par_val')

model2_d_omega_g[c('vg_label', 'ag_label', 'tsv_label')] <- str_split_fixed(as.character(model2_d_omega_g$agtsv), ', ', 3)
model2_d_omega_g$vg_label <- factor(model2_d_omega_g$vg_label, 
                                    levels = c("Unvaccinated", "Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3"))
model2_d_omega_g$ag_label <- factor(model2_d_omega_g$ag_label, 
                                    levels = c("0 - 4 yrs", "5 - 11 yrs", "12 - 17 yrs", "18 - 39 yrs", "40 - 60 yrs", "60+ yrs"))
model2_d_omega_g$tsv_label <- factor(model2_d_omega_g$tsv_label, 
                                     levels = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"))

design <- "
ADGI
BEH#
CF##
"

# creating vaccinated groups' plots
model2.d.o.agbyvg.vac.plot <- ggplot(filter(model2_d_omega_g, vg_label != "Unvaccinated"), aes(x=tsv_label, y = log10(par_val), fill = ag_label))+
  stat_summary(fun=median,
               geom="line",
               aes(group = ag_label, colour = ag_label),
               linewidth = 1.2,
               alpha = 0.7) +
  guides(colour = "none")+
  geom_violin(position = "identity", alpha = 0.7)+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        plot.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(face = "bold", angle = 90, hjust = 0.5, size = 12))+
  labs(fill = "Age group", 
       x = "Days since last vaccine dose",
       y = "log10(value) of immune proliferation (omega), per immune cell per day",
       title = "Delta infections") +  
  ggh4x::facet_manual(~ vg_label, design = design)+
  stat_summary(fun.data=mean_sdl,
               geom="pointrange",
               color="black",
               alpha = 0.4,
               size = 0.15,
               position = "identity") +
  ylim(-6.2,-3.35)

# creating unvaccinated group's plots
model2.d.o.agbyvg.unvac.plot <- ggplot(filter(model2_d_omega_g, vg_label == "Unvaccinated"), aes(x=tsv_label, y = log10(par_val), fill = ag_label))+
  geom_violin(position = "identity", alpha = 0.7, width = 0.5)+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(face = "bold", angle = 90, hjust = 0.5, size = 12))+
  labs(fill = "Age group", 
       x = "Age group",
       y = "log10(value) of immune proliferation (omega), per immune cell per day",
       title = "Delta infections") +  
  facet_wrap(~ vg_label, nrow = 1) +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange",
               color="black",
               alpha = 0.4,
               size = 0.15,
               position = "identity") +
  ylim(-6.2,-3.35)

byag.acrossvg <- (model2.d.o.agbyvg.vac.plot + 
                    ((model2.d.o.agbyvg.unvac.plot / plot_spacer() / guide_area()) +
                       plot_layout(guides = 'collect',
                                   heights = c(0.969, 0.5, 1.5)))) + plot_layout(widths = c(7, 1)) & scale_fill_discrete(limits = c("0 - 4 yrs", "5 - 11 yrs", "12 - 17 yrs", "18 - 39 yrs", "40 - 60 yrs", "60+ yrs")) & scale_colour_discrete(limits = c("0 - 4 yrs", "5 - 11 yrs", "12 - 17 yrs", "18 - 39 yrs", "40 - 60 yrs", "60+ yrs"))

ggsave("omega_fig4a_delta.png", 
       plot = wrap_elements(byag.acrossvg) + 
         labs(tag = "Log10(Immune proliferation), per immune cell per day") +
         theme(
           plot.tag = element_text(size = rel(1), angle = 90, face = "bold"),
           plot.tag.position = "left"
         ), 
       width = 11000, height = 8000, units = "px", dpi = 800)

# ===================================================================#
# Figure 4B
# Plots for Omega parameter (growth rate of immunity)
# Omicron dataset
# ===================================================================#
model2_o_omega <- data.frame(gamma = post.omicron.2$omega_mod)
colnames(model2_o_omega) <- ag_label_o

model2_o_omega_g <- model2_o_omega %>%
  gather(key = 'agtsv', value = 'par_val')

model2_o_omega_g[c('vg_label', 'ag_label', 'tsv_label')] <- str_split_fixed(as.character(model2_o_omega_g$agtsv), ', ', 3)
model2_o_omega_g$vg_label <- factor(model2_o_omega_g$vg_label, 
                                    levels = c("Unvaccinated", "Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3"))
model2_o_omega_g$ag_label <- factor(model2_o_omega_g$ag_label, 
                                    levels = c("0 - 4 yrs", "5 - 11 yrs", "12 - 17 yrs", "18 - 39 yrs", "40 - 60 yrs", "60+ yrs"))
model2_o_omega_g$tsv_label <- factor(model2_o_omega_g$tsv_label, 
                                     levels = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"))

## omicron
design <- "
ADGI
BEHJ
CF##
"
model2.o.o.agbyvg.vac.plot <- ggplot(filter(model2_o_omega_g, vg_label != "Unvaccinated"), aes(x=tsv_label, y = log10(par_val), fill = ag_label))+
  stat_summary(fun=median,
               geom="line",
               aes(group = ag_label, colour = ag_label),
               linewidth = 1.2,
               alpha = 0.7) +
  guides(colour = "none")+
  geom_violin(position = "identity", alpha = 0.7)+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        plot.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(face = "bold", angle = 90, hjust = 0.5, size = 12))+
  labs(fill = "Age group", 
       x = "Days since last vaccine dose",
       y = "log10(value) of immune proliferation (omega), per immune cell per day",
       title = "Omicron infections") +  
  ggh4x::facet_manual(~ vg_label, design = design)+
  stat_summary(fun.data=mean_sdl,
               geom="pointrange",
               color="black",
               alpha = 0.4,
               size = 0.15,
               position = "identity") +
  ylim(-6.2,-2.8)
model2.o.o.agbyvg.vac.plot

model2.o.o.agbyvg.unvac.plot <- ggplot(filter(model2_o_omega_g, vg_label == "Unvaccinated"), aes(x=tsv_label, y = log10(par_val), fill = ag_label))+
  geom_violin(position = "identity", alpha = 0.7, width = 0.5)+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(face = "bold", angle = 90, hjust = 0.5, size = 12))+
  labs(fill = "Age group", 
       x = "Age group",
       y = "log10(value) of immune proliferation (omega), per immune cell per day",
       title = "Omicron infections") +  
  facet_wrap(~ vg_label, nrow = 1) +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange",
               color="black",
               alpha = 0.4,
               size = 0.15,
               position = "identity") +
  ylim(-6.2,-2.8)
model2.o.o.agbyvg.unvac.plot

byag.acrossvg <- (model2.o.o.agbyvg.vac.plot + 
                    ((model2.o.o.agbyvg.unvac.plot / plot_spacer() / guide_area()) +
                       plot_layout(guides = 'collect',
                                   heights = c(0.969, 0.5, 1.5)))) + plot_layout(widths = c(7, 1)) & scale_fill_discrete(limits = c("0 - 4 yrs", "5 - 11 yrs", "12 - 17 yrs", "18 - 39 yrs", "40 - 60 yrs", "60+ yrs")) & scale_colour_discrete(limits = c("0 - 4 yrs", "5 - 11 yrs", "12 - 17 yrs", "18 - 39 yrs", "40 - 60 yrs", "60+ yrs"))

ggsave("omega_fig4b_omicron.png",
       plot = wrap_elements(byag.acrossvg) + 
         labs(tag = "Log10(Immune proliferation), per immune cell per day") +
         theme(
           plot.tag = element_text(size = rel(1), angle = 90, face = "bold"),
           plot.tag.position = "left"
         )
       ,
       width = 11000, height = 8000, units = "px", dpi = 800)

# ===================================================================#
# Figure 5
# Delta vs Omicron
# Gamma parameter (infected cell clearance by immunity)
# ===================================================================#
model2_d_gamma_g$variant <- c("Delta")
model2_o_gamma_g$variant <- c("Omicron")
model2_gamma_g <- rbind(model2_d_gamma_g, model2_o_gamma_g)
model2_gamma_g$tsv_label <- factor(model2_gamma_g$tsv_label, levels = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"))
model2_gamma_g$vg_label <- factor(model2_gamma_g$vg_label, levels = c("Unvaccinated", "Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3"))

design <- "
ABC
DEF
GH#
IJ#
"

model2.g.plot <- ggplot(filter(model2_gamma_g, vg_label != "Unvaccinated"), aes(x=tsv_label, y = par_val, fill = variant))+
  geom_violin(width = 0.7)+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        # plot.title = element_text(hjust=0.5,face="bold",size=12,colour = "blue"),
        plot.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(face = "bold", angle = 90, hjust = 0.5, size = 12))+
  labs(fill = "VOC", 
       x = "Days since last vaccine dose",
       y = "Infected cell clearance by immune response per immune cell per day") +  
  ggh4x::facet_manual(~ vg_label, design = design)+
  stat_summary(fun.data=mean_sdl,
               geom="pointrange",
               color="black",
               alpha = 0.5,
               size = 0.15,
               position = position_dodge(0.7)) +
  ylim(1,10)
model2.g.plot

model2.g.unvac.plot <- ggplot(filter(model2_gamma_g, vg_label == "Unvaccinated"), aes(x=variant, y = par_val, fill = variant))+
  geom_violin(width = 0.7)+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        # plot.title = element_text(hjust=0.5,face="bold",size=12,colour = "blue"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 12))+
  labs(fill = "VOC", 
       x = "VOC",
       y = "Infected cell clearance by immune response per immune cell per day") +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange",
               color="black",
               alpha = 0.5,
               size = 0.15,
               position = position_dodge(0.7)) +
  facet_grid(~vg_label) +
  ylim(1,10)
model2.g.unvac.plot

byvar.acrossvg <- (model2.g.plot + 
                     ((model2.g.unvac.plot / plot_spacer() / guide_area()) + 
                        plot_layout(guides = 'collect',
                                    heights = c(0.912, 1.5, 1.5))) + plot_layout(widths = c(5,1)) & scale_fill_discrete(limits = c("Delta", "Omicron")) )

ggsave("gamma_figure_notnormvariants.png", 
       plot =  wrap_elements(byvar.acrossvg) + 
         labs(tag = "Infected cell clearance by immune response per immune cell per day") +
         theme(
           plot.tag = element_text(size = rel(1), angle = 90),
           plot.tag.position = "left"
         ),
       width = 13000, height = 7000, units = "px", dpi = 800)

# ===================================================================#
# Figure 6
# Delta vs Omicron
# Omega parameter (growth rate of immunity)
# ===================================================================#
model2_d_omega_g$variant <- c("Delta")
model2_o_omega_g$variant <- c("Omicron")
model2_omega_g <- rbind(model2_d_omega_g, model2_o_omega_g)
model2_omega_g$vg_label <- factor(model2_omega_g$vg_label, levels = c("Unvaccinated", "Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3"))

model2.o.vac.plot <- ggplot(filter(model2_omega_g, vg_label != "Unvaccinated"), aes(x=tsv_label, y = log10(par_val), fill = variant))+
  stat_summary(fun=median,
               geom="line",
               aes(group = variant, colour = variant),
               linewidth = 1.2,
               alpha = 0.7) +
  guides(colour = "none")+
  geom_violin(position = "identity", alpha = 0.8)+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        # plot.title = element_text(hjust=0.5,face="bold",size=9,colour = "blue"),
        plot.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(face = "bold", angle = 90, hjust = 0.5, size = 12))+
  labs(fill = "Variants",
       x = "Days since last vaccine dose",
       y = "log10(immune proliferation) per immune cell per day",
       title = "log10(value) of immune proliferation, over time since vaccination for each infecting VOC, across vaccine group and age group") +  
  facet_grid(vars(vg_label), vars(ag_label))
model2.o.vac.plot

model2.o.unvac.plot <- ggplot(filter(model2_omega_g, vg_label == "Unvaccinated"), aes(x=tsv_label, y = log10(par_val), fill = variant))+
  geom_violin(position = "identity", alpha = 0.7, width = 0.3)+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_blank(),
        axis.text = element_text(face = "bold", angle = 90, hjust = 0.5, size = 12))+
  labs(fill = "Variants", 
       x = "Age group",
       y = "log10(value) of immune proliferation (omega), per immune cell per day",
       title = "Omicron infections") +  
  facet_grid(vars(vg_label), vars(ag_label))
model2.o.unvac.plot

byvar.acrossvgag <- (model2.o.unvac.plot / model2.o.vac.plot + 
                       plot_layout(guides = "collect",
                                   heights = c(1,8))) &
  scale_fill_discrete(limits = c("Delta", "Omicron")) &
  scale_color_discrete(limits = c("Delta", "Omicron"))

ggsave("omega_fig6_variants.png", 
       plot = wrap_elements(byvar.acrossvgag) + 
         labs(tag = "Log10(Immune proliferation), per immune cell per day") +
         theme(
           plot.tag = element_text(size = rel(1), angle = 90, face = "bold"),
           plot.tag.position = "left"
         )
       ,
       width = 8500, height = 12000, units = "px", dpi = 800)
