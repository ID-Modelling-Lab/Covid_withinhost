# ===================================================================#
# Figure 2 of main manuscript
# Viral load output based on median parameter values from fit 
# ===================================================================#
library(FME)
library(tidyverse)
library(gridExtra)
library(rstan)
library(readxl)
library(patchwork)
library(ggh4x)
library(ggnewscale)

# load rds files (model fits)
model2_delta_op <- readRDS("model2_delta_op.rds")
model2_omicron_op <- readRDS("model2_omicron_op.rds")

# load rdata file for vaccination regime and age group labelling
load("labels.RData")

post.delta.2 <- rstan::extract(model2_delta_op)
post.omicron.2 <-rstan::extract(model2_omicron_op)
params = c("beta_par", "c_par", "gamma_par", "omega_mod")

# setting up ode model
sc2sys <- function(t, state, para) {
  with(as.list(c(state, para)), {
    dT = 0.14*8*10^7 - 0.14*Tr - beta * V * Tr
    dI = beta * V * Tr - 0.14 * I - gamma * Z * I
    dV = 1.12*10^4 * I - c * V
    dZ = omega * I * Z
    list(c(dT, dI, dV, dZ))
  } ) }
state = c(Tr = 8*10^7,
          I = 0,
          V = 1,
          Z = 10)
times = seq(0, 20, by = 0.1)

# store median posterior run
med.delta.ode <- list()
med.omicron.ode <- list()

# ===========================================#
# Figure 2A
# Delta dataset
# ===========================================#
for (i in 1:127) {
  para.d = c(beta = median(post.delta.2$beta_par),
             c = median(post.delta.2$c_par),
             gamma = median(post.delta.2$gamma_par[,ct_samplesinagtsv.d$vaccgroup_range[i]]),
             omega = median(post.delta.2$omega_mod[,i]))
  solved = as.data.frame(ode(y = state, times = times, func = sc2sys, parms = para.d))
  med.delta.ode[[i]] <- solved
}

viral.output.d <- bind_cols(bind_rows(lapply(med.delta.ode, '[', 1)), bind_rows(lapply(med.delta.ode, '[', 4))) 
viral.output.d$vgagtsv <- rep(ag_label_d, each = 201)
viral.output.d[c('vg_label', 'ag_label', 'tsv_label')] <- str_split_fixed(as.character(viral.output.d$vgagtsv), ', ', 3)
viral.output.d$vg_label <- factor(viral.output.d$vg_label, levels = c("Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3", "Unvaccinated"))
viral.output.d$ag_label <- factor(viral.output.d$ag_label, levels = c("0 - 4 yrs", "5 - 11 yrs", "12 - 17 yrs", "18 - 39 yrs", "40 - 60 yrs", "60+ yrs"))
viral.output.d$tsv_label <- factor(viral.output.d$tsv_label, levels = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"))

palette <- rep(c())

viral.output.d.vac <- filter(viral.output.d, vg_label != "Unvaccinated") %>%
  mutate(mrna = case_when(vg_label == 'Pfizer 2' ~ '2-dose mRNA',
                          vg_label == 'Pfizer 3' ~ '3-dose mRNA',
                          vg_label == 'Moderna 2' ~ '2-dose mRNA',
                          vg_label == 'Moderna 3' ~ '3-dose mRNA',
                          vg_label == 'PPM' ~ '3-dose mRNA',
                          vg_label == 'MMP' ~ '3-dose mRNA',
                          vg_label == 'Sinovac 2' ~ '2-dose non-mRNA',
                          vg_label == 'Sinovac 3' ~ '3-dose non-mRNA',
                          vg_label == 'Sinopharm 2' ~ '2-dose non-mRNA',
                          vg_label == 'Sinopharm 3' ~ '3-dose non-mRNA'
  ))

# creating the unvaccinated reference groups for plotting
## reference unvaccinated line for 2-dose mrna groups & separating out the unvaccinate 0-4 and 5-11 yrs age grps
temp <- filter(viral.output.d, vg_label == "Unvaccinated")
temp$mrna <- c('2-dose mRNA')
temp$mrna[temp$ag_label == "0 - 4 yrs" | temp$ag_label == "5 - 11 yrs"] <- c('Unvaccinated')
viral.output.d.mrnasplit <- rbind(viral.output.d.vac, temp)

## reference unvaccinated line for 3-dose mrna groups
temp <- filter(viral.output.d, vg_label == "Unvaccinated")
temp$mrna <- c('3-dose mRNA')
temp <- filter(temp, (ag_label == '18 - 39 yrs' | ag_label == '40 - 60 yrs' | ag_label == '60+ yrs') )
viral.output.d.mrnasplit <- rbind(viral.output.d.mrnasplit, temp)

## reference unvaccinated line for 2-dose non-mrna groups
temp <- filter(viral.output.d, vg_label == "Unvaccinated")
temp$mrna <- c('2-dose non-mRNA')
temp <- filter(temp, (ag_label == '18 - 39 yrs' | ag_label == '40 - 60 yrs' | ag_label == '60+ yrs') )
viral.output.d.mrnasplit <- rbind(viral.output.d.mrnasplit, temp)

## reference unvaccinated line for 3-dose non-mrna groups
temp <- filter(viral.output.d, vg_label == "Unvaccinated")
temp$mrna <- c('3-dose non-mRNA')
temp <- filter(temp, (ag_label == '18 - 39 yrs' | ag_label == '40 - 60 yrs') )
viral.output.d.mrnasplit <- rbind(viral.output.d.mrnasplit, temp)

viral.output.d.mrnasplit$mrna <- factor(viral.output.d.mrnasplit$mrna, levels = c("Unvaccinated", "2-dose non-mRNA", "3-dose non-mRNA", "2-dose mRNA", "3-dose mRNA"))

design <- "
ABC#
DEFG
HIJK
L#MN
"

viral.output.d.plot <- ggplot(data = NULL, aes(x=time, y = log10(V))) + 
  facet_manual(~ag_label + mrna,
               scales = "free_x",
               design = design) +
  geom_line(data = filter(viral.output.d.mrnasplit, vg_label == "Sinopharm 2"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 1) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#C7EEFF", "#ABD4EF", "#8EB9DF", "#729FCF","#5584C0","#396AB0","#1C4FA0","#003590", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.d.mrnasplit, vg_label == "Sinovac 2"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 0.75) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#FFC7C7", "#EDABAB", "#DB8F8F", "#C97373","#B85656","#A63A3A","#941E1E","#820202", "#000000"))  +
  new_scale_colour() +
  geom_line(data = filter(viral.output.d.mrnasplit, vg_label == "Sinovac 3"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 1) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#FFC7C7", "#EDABAB", "#DB8F8F", "#C97373","#B85656","#A63A3A","#941E1E","#820202", "#000000"))  +
  new_scale_colour() +
  geom_line(data = filter(viral.output.d.mrnasplit, vg_label == "Moderna 2"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 1) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#AEF78C", "#96E178", "#7ECB64", "#66B550","#4E9F3C","#368928","#1E7314","#065D00", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.d.mrnasplit, vg_label == "Pfizer 2"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 0.8) +
  scale_color_manual(breaks =  c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#FFD666", "#FAC357", "#F5B049", "#F09D3A","#EA892C","#E5761D","#E0630F","#DB5000", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.d.mrnasplit, vg_label == "Moderna 3"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 1) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#AEF78C", "#96E178", "#7ECB64", "#66B550","#4E9F3C","#368928","#1E7314","#065D00", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.d.mrnasplit, vg_label == "Pfizer 3"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 0.87) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#FFD666", "#FAC357", "#F5B049", "#F09D3A","#EA892C","#E5761D","#E0630F","#DB5000", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.d.mrnasplit, vg_label == "PPM"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 0.87) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#DEAAF1", "#C892DC","#B279C6","#9C61B1","#85499C","#6F3187","#591871","#43005C","#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.d.mrnasplit, vg_label == "MMP"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 0.8) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#ACCFF2", "#93B2DD", "#7B96C7", "#6279B2","#4A5C9C","#313F87","#192371","#00065C", "#000000")) +
  geom_line(data = filter(viral.output.d.mrnasplit, vg_label == "Unvaccinated"),
            linewidth = 1.1, alpha = 0.8) + 
  theme_bw()+
  theme(axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        plot.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 12),
        strip.text.x = element_text(face = "bold", size = 11)) +
  labs(x = "Day of infection",
       y = "log10(virus)")

ggsave("viraloutput_delta_mrnasplit_dosesplit.png", 
       plot = wrap_elements(viral.output.d.plot),
       width = 8000, height = 7500, units = "px", dpi = 700)

# ===========================================#
# Figure 2A
# Omicron dataset
# ===========================================#
for (i in 1:128) {
  para.o = c(beta = median(post.omicron.2$beta_par),
             c = median(post.omicron.2$c_par),
             gamma = median(post.omicron.2$gamma_par[,ct_samplesinagtsv.o$vaccgroup_range[i]]),
             omega = median(post.omicron.2$omega_mod[,i]))
  solved = as.data.frame(ode(y = state, times = times, func = sc2sys, parms = para.o))
  med.omicron.ode[[i]] <- solved
}

viral.output.o <- bind_cols(bind_rows(lapply(med.omicron.ode, '[', 1)), bind_rows(lapply(med.omicron.ode, '[', 4))) 
viral.output.o$vgagtsv <- rep(ag_label_o, each = 201)
viral.output.o[c('vg_label', 'ag_label', 'tsv_label')] <- str_split_fixed(as.character(viral.output.o$vgagtsv), ', ', 3)
viral.output.o$vg_label <- factor(viral.output.o$vg_label, levels = c("Pfizer 2", "Pfizer 3", "PPM", "Moderna 2", "Moderna 3", "MMP", "Sinovac 2", "Sinovac 3", "Sinopharm 2", "Sinopharm 3", "Unvaccinated"))
viral.output.o$ag_label <- factor(viral.output.o$ag_label, levels = c("0 - 4 yrs", "5 - 11 yrs", "12 - 17 yrs", "18 - 39 yrs", "40 - 60 yrs", "60+ yrs"))
viral.output.o$tsv_label <- factor(viral.output.o$tsv_label, levels = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"))

palette <- rep(c())

viral.output.o.vac <- filter(viral.output.o, vg_label != "Unvaccinated") %>%
  mutate(mrna = case_when(vg_label == 'Pfizer 2' ~ '2-dose mRNA',
                          vg_label == 'Pfizer 3' ~ '3-dose mRNA',
                          vg_label == 'Moderna 2' ~ '2-dose mRNA',
                          vg_label == 'Moderna 3' ~ '3-dose mRNA',
                          vg_label == 'PPM' ~ '3-dose mRNA',
                          vg_label == 'MMP' ~ '3-dose mRNA',
                          vg_label == 'Sinovac 2' ~ '2-dose non-mRNA',
                          vg_label == 'Sinovac 3' ~ '3-dose non-mRNA',
                          vg_label == 'Sinopharm 2' ~ '2-dose non-mRNA',
                          vg_label == 'Sinopharm 3' ~ '3-dose non-mRNA'
  ))

# creating the unvaccinated reference groups for plotting
## reference unvaccinated line for 2-dose mrna groups & separating out the unvaccinated 0-4 yrs age grp
temp <- filter(viral.output.o, vg_label == "Unvaccinated")
temp$mrna <- c('2-dose mRNA')
temp$mrna[temp$ag_label == "0 - 4 yrs"] <- c('Unvaccinated')
viral.output.o.mrnasplit <- rbind(viral.output.o.vac, temp)

## reference unvaccinated line for 3-dose mrna groups
temp <- filter(viral.output.o, vg_label == "Unvaccinated")
temp$mrna <- c('3-dose mRNA')
temp <- filter(temp, (ag_label == '12 - 17 yrs' | ag_label == '18 - 39 yrs' | ag_label == '40 - 60 yrs' | ag_label == '60+ yrs') )
viral.output.o.mrnasplit <- rbind(viral.output.o.mrnasplit, temp)

## reference unvaccinated line for 2-dose non-mrna groups
temp <- filter(viral.output.o, vg_label == "Unvaccinated")
temp$mrna <- c('2-dose non-mRNA')
temp <- filter(temp, (ag_label == '18 - 39 yrs' | ag_label == '40 - 60 yrs' | ag_label == '60+ yrs') )
viral.output.o.mrnasplit <- rbind(viral.output.o.mrnasplit, temp)

## reference unvaccinated line for 3-dose non-mrna groups
temp <- filter(viral.output.o, vg_label == "Unvaccinated")
temp$mrna <- c('3-dose non-mRNA')
temp <- filter(temp, (ag_label == '18 - 39 yrs' | ag_label == '40 - 60 yrs' | ag_label == '60+ yrs') )
viral.output.o.mrnasplit <- rbind(viral.output.o.mrnasplit, temp)

viral.output.o.mrnasplit$mrna <- factor(viral.output.o.mrnasplit$mrna, levels = c("Unvaccinated", "2-dose non-mRNA", "3-dose non-mRNA", "2-dose mRNA", "3-dose mRNA"))

design <- "
ABCD
EFGH
IJKL
MNOP
"

viral.output.o.plot <- ggplot(data = NULL, aes(x=time, y = log10(V))) + 
  facet_manual(~ag_label + mrna,
               scales = "free_x",
               design = design) +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "Sinopharm 2"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 1) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#C7EEFF", "#ABD4EF", "#8EB9DF", "#729FCF","#5584C0","#396AB0","#1C4FA0","#003590", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "Sinovac 2"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 0.75) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#FFC7C7", "#EDABAB", "#DB8F8F", "#C97373","#B85656","#A63A3A","#941E1E","#820202", "#000000"))  +
  new_scale_colour() +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "Sinopharm 3"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 1) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#C7EEFF", "#ABD4EF", "#8EB9DF", "#729FCF","#5584C0","#396AB0","#1C4FA0","#003590", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "Sinovac 3"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 1) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#FFC7C7", "#EDABAB", "#DB8F8F", "#C97373","#B85656","#A63A3A","#941E1E","#820202", "#000000"))  +
  new_scale_colour() +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "Moderna 2"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 1) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#AEF78C", "#96E178", "#7ECB64", "#66B550","#4E9F3C","#368928","#1E7314","#065D00", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "Pfizer 2"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 0.8) +
  scale_color_manual(breaks =  c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#FFD666", "#FAC357", "#F5B049", "#F09D3A","#EA892C","#E5761D","#E0630F","#DB5000", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "Moderna 3"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 1) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#AEF78C", "#96E178", "#7ECB64", "#66B550","#4E9F3C","#368928","#1E7314","#065D00", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "Pfizer 3"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 0.87) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#FFD666", "#FAC357", "#F5B049", "#F09D3A","#EA892C","#E5761D","#E0630F","#DB5000", "#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "PPM"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 0.87) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#DEAAF1", "#C892DC","#B279C6","#9C61B1","#85499C","#6F3187","#591871","#43005C","#000000")) +
  new_scale_colour() +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "MMP"), 
            aes(colour = tsv_label),
            linewidth = 1.35, alpha = 0.8) +
  scale_color_manual(breaks = c("(0,7]", "(7,14]", "(14,28]", "(28,60]", "(60,90]", "(90,180]", "(180,270]", "(270,360]", "unvax"),
                     values = c("#ACCFF2", "#93B2DD", "#7B96C7", "#6279B2","#4A5C9C","#313F87","#192371","#00065C", "#000000")) +
  geom_line(data = filter(viral.output.o.mrnasplit, vg_label == "Unvaccinated"),
            linewidth = 1.1, alpha = 0.8) + 
  theme_bw()+
  theme(axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        plot.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(face = "bold", angle = 0, hjust = 0.5, size = 12),
        strip.text.x = element_text(face = "bold", size = 11)) +
  labs(x = "Day of infection",
       y = "log10(virus)")

ggsave("viraloutput_omicron_mrnasplit_dosesplit.png", 
       plot = wrap_elements(viral.output.o.plot),
       width = 8000, height = 7500, units = "px", dpi = 700)
