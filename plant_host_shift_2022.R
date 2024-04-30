##################################################################
#                                                                #
#   Plant-host shift, spatial persistence, and the viability     # 
#            of an invasive insect population                    #
#                                                                #        
#          Isabelle Bueno Silva, Blake McGrane-Corrigan,         #
#             Oliver Mason, Rafael de Andrade Moral              #  
#                 and Wesley Augusto Conde Godoy                 #
#                                                                #
#  https://www.biorxiv.org/content/10.1101/2021.09.20.461112v5   #
#                                                                #
##################################################################

# Clear the workspace
rm(list=ls())

require(tidyverse)
require(JM)
require(ggplot2)
require(gridExtra)
require(patchwork)
require(hnp)
require(coefplot2)
require(latticeExtra)
require(survival)
require(lme4)

####################################
#                                  #
# Generate n complementary colours #
#                                  #
####################################

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colours <- gg_color_hue(7)

#############################################
#                                           # 
# Calculate F_i and S_i for dispersal model #
#                                           #
#############################################

total <- read.csv("Total.csv")
total %>%
  group_by(treatment) %>%
  filter(egg_laying_days > 0) %>%
  summarise("Fecundity" = mean(na.omit(total_eggs/egg_laying_days)),
            "Survival" = mean(viable/ifelse(total_eggs == 0, 1, total_eggs)))

########################################
#                                      #
# Simulations and Bifurcation daysgrams #
#                                      #
########################################

# empirical values for R S and F are:
# S_1 = 0.335, R_1 = 0.5, F_1 = 9.43
# S_2 = 0.216, R_2 = 0.5, F_2 = 21.7

########################################################
# Plot of trajectories with varying initial conditions #
#             with dispersal rates fixed               #
########################################################

time = 20; d1 = 0.6; d2 = 0.4; alpha = 0.0064;
x1 <- x2 <- numeric(time)

# R_iF_S_i values got from experiments 
RFS1 <- 1.5795
RFS2 <- 2.3436

# Create a sequence of initial conditions between 0 and 30
initialk1_grid <- seq(0,150, length.out = 20)
initialk2_grid <- seq(0,150, length.out = 20)

# Every possible combination of initial conditions
initialk <-  expand.grid(initialk1 = initialk1_grid, 
                         initialk2 = initialk2_grid)

# Trajectories combined into matrix with columns the trajectories for each combination
traj_initial1k <- traj_initial2k <- matrix(NA, ncol = nrow(initialk), nrow = time)

# For each combination simulate model and add trajectory to matrix traj
for(i in 1:nrow(initialk)) {
  
  # Initialise
  x1[1] <- initialk[i,]$initialk1
  x2[1] <- initialk[i,]$initialk2
  
  # Iterate difference equations
  for(t in 2:time) {
    
    g1 <- RFS1 * x1[t-1] * exp(-alpha * x1[t-1])
    g2 <- RFS2 * x2[t-1] * exp(-alpha * x2[t-1])
    x1[t] <- (1 - d1) * g1 + d2 * g2 
    x2[t] <- (1 - d2) * g2 + d1 * g1
  }
  
  # Add trajectory to matrix
  traj_initial1k[,i] <- c(x1)
  traj_initial2k[,i] <- c(x2)
  
}

sim1.initial1k <- as.data.frame(traj_initial1k)
sim1.initial2k <- as.data.frame(traj_initial2k)

# Take some of these trajectories to plot
sim1.initial1k <- sim1.initial1k[, seq(5, 400, 23)]
sim1.initial2k <- sim1.initial2k[, seq(5, 400, 23)]

sim1.initial1k$index <- factor("Patch 1")
sim1.initial1k$Generation <- 1:time

sim1.initial1k <- sim1.initial1k %>% 
  pivot_longer(cols = starts_with("V"),
               values_to = "Trajectory")

sim1.initial2k$index <- factor("Patch 2")
sim1.initial2k$Generation <- 1:time

sim1.initial2k <- sim1.initial2k %>% 
  pivot_longer(cols = starts_with("V"),
               values_to = "Trajectory")

plotx1 <- sim1.initial1k %>%
  ggplot(aes(x = Generation, 
             y = Trajectory, 
             group = name)) +
  theme_bw() +
  geom_line(size = 0.2, 
            colour = colours[2]) +
  ylab("Population Size")+
  xlab("Generation")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "#BDBDBD", 
                                    fill = NA, 
                                    linetype = "dotted"),
        aspect.ratio = 0.8, 
        axis.text = element_text(colour = 1, 
                                 size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "#BDBDBD", 
                                             fill = NA, 
                                             linetype = "dotted"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        legend.text = element_text(face = "plain",
                                   hjust = 0, 
                                   size = 10))

plotx2 <- sim1.initial2k %>%
  ggplot(aes(x = Generation, 
             y = Trajectory, 
             group = name)) +
  theme_bw() +
  geom_line(size = 0.2, 
            colour = colours[3]) +
  ylab("Population Size")+
  xlab("Generation")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "#BDBDBD", 
                                    fill = NA, 
                                    linetype = "dotted"),
        aspect.ratio = 0.8, 
        axis.text = element_text(colour = 1, 
                                 size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "#BDBDBD", 
                                             fill = NA, 
                                             linetype = "dotted"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        legend.text = element_text(face = "plain",
                                   hjust = 0, 
                                   size = 10))

cowplot::plot_grid(plotx1, 
                   plotx2,
                   nrow = 1,
                   labels = "",
                   label_colour = "#70976c",
                   align = "v")

########################################################
# Plot model trajectories with varying dispersal rates #
########################################################

time = 75; initial_x1 = 10; initial_x2 = 10; alpha = 0.0064;
x1 <- x2 <- numeric(time)

# Create a sequence of mu_i values between 0.01 and 0.05
mu1_grid <- round(seq(0.01,0.05, length.out = 25),5)
mu2_grid <- round(seq(0.01,0.05, length.out = 25),5)

# Get every possible combination of these mu_i values
mu <- expand.grid(muk1 = mu1_grid, 
                  muk2 = mu2_grid)

# Create matrix with columns the trajectories for each combination
traj_mu1 <- traj_mu2 <- matrix(NA, ncol = nrow(mu), nrow = time)

# R_iF_iS_i values got from experiments 
RFS1 <- 1.5795
RFS2 <- 2.3436

# Minimum of 10% dispersal permitted between patches
r1 <- 0.9
r2 <- 0.9

# For each combination of mu_i values simulate dispersal model
for(i in 1:nrow(mu)) {
  
  # Initialise
  x1[1] <- initial_x1
  x2[1] <- initial_x2
  
  # Iterate difference equations
  for(t in 2:time) {
    
    g1 <- RFS1 * x1[t-1] * exp(-alpha * x1[t-1])
    g2 <- RFS2 * x2[t-1] * exp(-alpha * x2[t-1])
    D1 <- 1-r1 * exp(-mu[i,]$muk1 * x1[t-1])
    D2 <- 1-r2 * exp(-mu[i,]$muk2 * x2[t-1])
    x1[t] <- (1 - D1) * g1 + D2 * g2 
    x2[t] <- (1 - D2) * g2 + D1 * g1
  }
  
  # Add trajectory to trajectory matrix
  traj_mu1[,i] <- c(x1)
  traj_mu2[,i] <- c(x2)
  
}

sim1.mu1 <- as.data.frame(traj_mu1)
sim1.mu2 <- as.data.frame(traj_mu2)

# Take some of these trajectories to plot
sim1.mu1 <- sim1.mu1[, round(seq(5, 625, length.out = 5))]
sim1.mu2 <- sim1.mu2[, round(seq(5, 625, length.out = 5))]

sim1.mu1$index <- factor("Patch 1")
sim1.mu1$Generation <- 1:time
sim1.mu1 <- sim1.mu1 %>% 
  pivot_longer(cols = starts_with("V"), 
               values_to = "Trajectory")
sim1.mu2$index <- factor("Patch 2")
sim1.mu2$Generation <- 1:time

sim1.mu2 <- sim1.mu2 %>% 
  pivot_longer(cols = starts_with("V"), 
               values_to = "Trajectory")

plotmu1 <- sim1.mu1 %>%
  ggplot(aes(x = Generation, 
             y = Trajectory, 
             group = name)) +
  theme_bw() +
  geom_line(size = 0.2, 
            colour = colours[2]) +
  ylab("Population Size") +
  xlab("Generation") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "#BDBDBD", 
                                    fill = NA, 
                                    linetype = "dotted"),
        aspect.ratio = 0.8, 
        axis.text = element_text(colour = 1,
                                 size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "#BDBDBD", 
                                             fill = NA, 
                                             linetype = "dotted"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        legend.text = element_text(face = "plain",
                                   hjust = 0, size = 10))

plotmu2 <- sim1.mu2 %>%
  ggplot(aes(x = Generation, 
             y = Trajectory, 
             group = name)) +
  theme_bw() +
  geom_line(size = 0.2, 
            colour = colours[3]) +
  ylab("Population Size") +
  xlab("Generation") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "#BDBDBD", 
                                    fill = NA, 
                                    linetype = "dotted"),
        aspect.ratio = 0.8, 
        axis.text = element_text(colour = 1, 
                                 size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "#BDBDBD", 
                                             fill = NA, 
                                             linetype = "dotted"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        legend.text = element_text(face = "plain",
                                   hjust = 0, size = 10))

b <- cowplot::plot_grid(plotmu1, plotmu2, 
                   nrow = 1, 
                   labels = "", 
                   label_colour = "#70976c", 
                   align = "v")

###########################################################
# Bifurcation daysgrams for mu_i and  R_iF_iS_i parameters #
###########################################################

patch_dynamics <- function(time, x1_init, x2_init, mu1, mu2, a, b, alpha, r1, r2) {
  
  # Initialise
  x1 <- x2 <- numeric(time)
  x1[1] <- x1_init
  x2[1] <- x2_init
  
  # Iterate difference equations
  for(t in 2:(time)) {
    
    g1 <- a*x1[t-1]*exp(-alpha*x1[t-1])
    g2 <- b*x2[t-1]*exp(-alpha*x2[t-1])
    D1 <- 1 - r1*exp(-mu1 * x1[t-1])
    D2 <- 1 - r2*exp(-mu2 * x2[t-1])
    x1[t] <- (1 - D1) * g1 + D2 * g2 
    x2[t] <- (1 - D2) * g2 + D1 * g1
  }
  
  ret <- tibble(x1, x2, time = 1:time)
  
  return(ret)
}

# Bifurcation plot for mu_i parameters
bifurcation_mu <- function(total_time = 1000, bif_time = 100, 
                           x1_init = 10, x2_init = 10, 
                           a = 20, b = 24, alpha = 0.0064,
                           r1 = 0.9, r2 = 0.9, 
                           mu_range = c(0.001, 0.1), 
                           # increase mu_length for better resolution
                           mu_length = 150) {
  
  sim <- list()
  # Initialise index
  i <- 1
  j <- 1
  for(mu_1 in seq(from = mu_range[1], 
                  to = mu_range[2], 
                  length = mu_length)) {
    for(mu_2 in seq(from = mu_range[1], 
                    to = mu_range[2], 
                    length = mu_length)) {
      # Simulate model for each combination of mu_i in mu_range
      sim[[i]] <- patch_dynamics(time = total_time,  
                                 x1_init = x1_init, x2_init = x2_init, 
                                 a = a, b = b, alpha = alpha, mu1 = mu_1, 
                                 mu2 = mu_2, r1 = r1, r2 =r2)
      # Take out last (total_time - bif_time) population values
      sim[[i]] <- sim[[i]] %>%
        filter(time > (total_time - bif_time)) %>%
        dplyr::select(-time) %>%
        mutate(mu_1 = mu_1,
               mu_2 = mu_2) %>%
        pivot_longer(cols = 1:2,
                     names_to = "Patch",
                     values_to = "N")
      i <- i + 1
    }
    # Print iterations
    cat("Iteration ", 
        j, 
        " out of ", 
        mu_length, "\n")
    j <- j + 1
  }
  
  ret <- do.call(rbind, sim)
  return(ret)
}

# bifurcation_mu() takes sufficient time to run
# Check ".RData" file to create bifurcation plot

# Run bifurcation function:
bif_mu <- bifurcation_mu()

# Save RData file
save(bif_mu, file = "bif_mu.RData", compress = T)

# Note that the RData file below will produce an image of lesser resolution
# To increase resolution make mu_length larger
load("bif_mu.RData")

# Return unique population values to ensure we have reached equilibrium or limit cycle
bif_mu_tidy <- bif_mu %>%
  group_by(mu_1, mu_2, Patch) %>%
  filter(Patch == "x2") %>%
  summarise(n = length(unique(round(N, 4)))) %>%
  mutate(n_bin = cut(n, breaks = c(-Inf, 1, 2, 4, 8, 16, 32, 100)))

# Periods of limit cycles
levels(bif_mu_tidy$n_bin) <- c("1","2","4","8","16","32","> 32")

plot_bif_mu <- bif_mu_tidy %>%
  ggplot(aes(x = mu_1, 
             y = mu_2, 
             fill = n_bin)) +
  facet_wrap(~ Patch) +
  theme_bw() +
  geom_tile() +
  scale_fill_brewer(name = "Periodicity",
                    palette = "Greens") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "#BDBDBD", 
                                    fill = NA, 
                                    linetype = "dotted"),
        aspect.ratio = 1, 
        axis.text = element_text(colour = 1, 
                                 size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "#BDBDBD", 
                                             fill = NA, 
                                             linetype = "dotted"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "right",
        legend.text = element_text(face = "plain",
                                   hjust = 0, 
                                   size = 10)) +
  xlab(expression(mu[1])) +
  ylab(expression(mu[2]))

RFS_func <- function(time, x1_init, x2_init, alpha, mu1, mu2, a, b, r1, r2){
  x1 <- x2 <- NULL
  # Initialise 
  x1[1] <- x1_init
  x2[1] <- x2_init
  
  ## Drawing values for density
  for(t in 2:time) {
    
    g1 <- a*x1[t-1]*exp(-alpha*x1[t-1])
    g2 <- b*x2[t-1]*exp(-alpha*x2[t-1])
    D1 <- 1 - r1*exp(-mu1 * x1[t-1])
    D2 <- 1 - r2*exp(-mu2 * x2[t-1])
    x1[t] <- (1 - D1) * g1 + D2 * g2 
    x2[t] <- (1 - D2) * g2 + D1 * g1
  }
  
  out <- tibble(x1, x2, time = 1:time)
  
  return(out)
  
}

# Bifurcation plot for RFS
bifurcation_RFS <- function(total_time = 1000, bif_time = 100, x1_init = 10, 
                            x2_init = 10, alpha = 0.0064, mu1 = 0.2, mu2 = 0.3,
                            r1 = 0.9, r2 = 0.9, RFS_range = c(0.01, 25), 
                            # change RFS_length to increase resolution
                            RFS_length = 150) {
  
  sim <- list()
  # Initialise index
  i <- 1
  j <- 1
  for(RFS_1 in seq(from = RFS_range[1], to = RFS_range[2], length = RFS_length)) {
    for(RFS_2 in seq(from = RFS_range[1], to = RFS_range[2], length = RFS_length)) {
      # Simulate model for each combination of R_iF_iS_i in RFS_range
      sim[[i]] <- RFS_func(time = total_time,  x1_init = x1_init,
                           x2_init = x2_init, alpha = alpha, 
                           mu = mu1, mu2 = mu2, a = RFS_1, b = RFS_2,
                           r1 = r1, r2 = r2)
      # Take out last (total_time - bif_time) population values
      sim[[i]] <- sim[[i]] %>%
        filter(time > (total_time - bif_time)) %>%
        dplyr::select(-time) %>%
        mutate(RFS_1 = RFS_1, 
               RFS_2 = RFS_2) %>%
        pivot_longer(cols = 1:2,
                     names_to = "Patch",
                     values_to = "N")
      i <- i + 1
    }
    # Print iterations
    cat("Iteration ", j, " out of ", RFS_length, "\n")
    j <- j + 1
  }
  
  ret <- do.call(rbind, sim)
  return(ret)
}

# bifurcation_mu() takes sufficient time to run
# Check ".RData" file to create bifurcation plot

# Run bifurcation function:
bif_RFS <- bifurcation_RFS()
save(bif_RFS, file = "bif_RFS.RData")

load("bif_RFS.RData")

# Return unique population values to ensure we have reached equilibrium or limit cycle
bif_RFS_tidy <- bif_RFS %>%
  group_by(RFS_1, RFS_2, Patch) %>%
  filter(Patch == "x2") %>%
  summarise(n = length(unique(round(N, 4)))) %>%
  mutate(n_bin = cut(n, breaks = c(-Inf, 1, 2, 4, 8, 16, 32, 100)))

# Periods of limit cycles
levels(bif_RFS_tidy$n_bin) <- c("1","2","4","8","16","32","> 32")

plot_bif_RFS <- bif_RFS_tidy %>%
  ggplot(aes(x = RFS_1, y = RFS_2, fill = n_bin)) +
  facet_wrap(~ Patch) +
  theme_bw() +
  geom_tile() +
  scale_fill_brewer(name = "Periodicity",
                    palette = "Greens") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "#BDBDBD", fill = NA, linetype = "dotted"),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "#BDBDBD", fill = NA, linetype = "dotted"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        legend.text = element_text(face = "plain",
                                   hjust = 0, size = 10)) +
  xlab(expression("R"[1]*"F"[1]*"S"[1])) +
  ylab(expression("R"[2]*"F"[2]*"S"[2]))

plot_bif_RFS + plot_bif_mu


############################################################################################
# Take values from R_iF_iS_i bifurcation plot and plot trajectories of the dispersal model # 
#                 corresponding to different limit cycle periodicities                     #
############################################################################################

# Let mu1 = 0.2 and mu2 = 0.3 so we have moderate dispersal rates

# Period 1 - nonnegative equilibrium (lightest green)
bifsim1 <- RFS_func(time = 50, x1_init = 10, 
                    x2_init = 10, alpha = 0.0064, mu1 = 0.2, mu2 = 0.3,
                    r1 = 0.9, r2 = 0.9, a = 5, b = 6)

# Period 2 - oscillation around equilibrium of period 2 (second darkest green)
bifsim2 <- RFS_func(time = 50, x1_init = 10, 
                    x2_init = 10, alpha = 0.0064, mu1 = 0.2, mu2 = 0.3,
                    r1 = 0.9, r2 = 0.9, a = 12.5, b = 5)

# Period 4 - oscillation around equilibrium of period 4 (third darkest green)
bifsim3 <- RFS_func(time = 50, x1_init = 10, 
                    x2_init = 10, alpha = 0.0064, mu1 = 0.2, mu2 = 0.3,
                    r1 = 0.9, r2 = 0.9, a = 17, b = 5.5)

# Period > 32 - oscillation around equilibrium of period greater than 32 (deterministic chaos)
# (darkest green)
bifsim4 <- RFS_func(time = 50, x1_init = 10, 
                    x2_init = 10, alpha = 0.0064, mu1 = 0.2, mu2 = 0.3,
                    r1 = 0.9, r2 = 0.9, a = 20, b = 20)

sim_bif <- rbind(bifsim1, bifsim2, bifsim3, bifsim4)

sim_bif <- sim_bif %>%
  rename(Strawberry = x1,
         Raspberry = x2,
         Generation = time)

sim_bif$panel <- factor(rep(c(
  "Period 1",
  "Period 2",
  "Period 4",
  "Period > 32"), each = 50),
  levels = c(
    "Period 1",
    "Period 2",
    "Period 4",
    "Period > 32")
)
sim_final <- sim_bif %>%
  pivot_longer(1:2,
               names_to = "Patch",
               values_to = "N")

sim_final %>%
  ggplot(aes(x = Generation, y = N, col = Patch)) +
  theme_bw() +
  geom_line() +
  scale_colour_manual(values=c(colours[2],colours[3])) +
  facet_wrap(~ panel, scales = "free_y") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "#BDBDBD", fill = NA, linetype = "dotted"),
        legend.spacing.y = unit(0, "mm"),
        aspect.ratio = 0.5, axis.text = element_text(colour = 1, size = 12),
        strip.text = element_text(face = "plain", color = "black",
                                  hjust = 0, size = 8),
        strip.background = element_rect(fill = "white", color = "#BDBDBD", linetype = "dotted"),
        legend.box.background = element_rect(colour = "#BDBDBD", fill = NA, linetype = "dotted"),
        legend.title = element_blank(),
        legend.text = element_text(face = "plain", hjust = 0, size = 10),
        legend.position = "none",
        legend.background = element_rect(fill = "white", color = "#BDBDBD", linetype = "dotted")) +
  ylab("Population Size")

########################
#                      #
# Statistical analyses #
#                      #
########################

###########################
# Quasi-Poisson modelling #
###########################

# Standard error
se <- function(x) sd(x)/sqrt(length(x))

total$block <- as.factor(total$block)
total$generation <- as.factor(total$generation)
total$treatment <- as.factor(total$treatment)
total$nf <- ifelse(total$generation == "F1", 10, 5)

# Fit quasi-Poisson GLM with total_eggs as response
# Block (intercept) and covariate treatment, including an offset for if F1 or F2 
# We take into account each generation by multiplying it by treatment
fit1 <- glm(total_eggs ~ block + generation * treatment + offset(log(nf)),
            family = quasipoisson, data = total)

# Half normal plot daysgnostic
hnp(fit1)

# F test for model fit 
anova(fit1, test = "F")

# Plot confidence intervals for model fit1
coefplot2(update(fit1, . ~ 0 + generation:treatment),
          intercept = T)

# Fit quasi-Poisson GLM with total_eggs as response for F1 population only
# Block (intercept) and covariate treatment (simple linear model)
fit1.2 <- glm(total_eggs ~ block + treatment, family = quasipoisson,
              data = subset(total, generation == "F1"))

# F test for model fit
anova(fit1.2, test = "F")

# Fit quasi-Poisson GLM with total_eggs as response for F2 population only
# Block (intercept) and covariate treatment (simple linear model)
fit1.3 <- glm(total_eggs ~ block + treatment, family = quasipoisson,
              data = subset(total, generation == "F2"))

# F test for model fit
anova(fit1.3, test = "F")

# Fit quasi-Poisson GLM with total_eggs as response for strawberry population
# Block (intercept) and covariate generation (simple linear model + offset)
fit1.4 <- glm(total_eggs ~ block + generation + offset(log(nf)),
              family = quasipoisson,
              data = subset(total, treatment == "strawberry"))

# F test for model fit
anova(fit1.4, test = "F")

# Fit quasi-Poisson GLM with total_eggs as response for raspberry population
# Block (intercept) and covariate generation (simple linear model + offset)
fit1.5 <- glm(total_eggs ~ block + generation + offset(log(nf)),
              family = quasipoisson,
              data = subset(total, treatment == "raspberry"))

# F test for model fit
anova(fit1.5, test = "F")

# Get mean and SE for each model with generation and treatment
round(cbind(with(total, tapply(total_eggs/nf, generation:treatment, mean)),
            with(total, tapply(total_eggs/nf, generation:treatment, se))), 4)

# Plot parameter estimates with confidence intervals
bwplot(total_eggs/nf ~ generation:treatment, data = total)

# Fit quasi-Poisson GLM with egg_laying_day as response using whole dataset
fit1 <- glm(egg_laying_days ~ block + generation * treatment,
            family = quasipoisson, data = total)

hnp(fit1)

anova(fit1, test = "F")

round(cbind(with(total, tapply(egg_laying_days, generation, mean)),
            with(total, tapply(egg_laying_days, generation, se))), 4)

round(cbind(with(total, tapply(egg_laying_days, treatment, mean)),
            with(total, tapply(egg_laying_days, treatment, se))), 4)

##################################################
# Survival and Cox-proportional hazards analyses #
##################################################

# Survival analysis with survival_days as response for total population
# and generation and treatment as covariates
km1 <- survfit(Surv(survival_days) ~ generation + treatment, data = total)
plot(km1)

# Cox-proportional hazards model with survival_days as response for total population
# Block (intercept) and covariate treatment taking into account generation
cox1 <- coxph(Surv(survival_days) ~ block + generation * treatment, data = total)
anova(cox1)

# Chi-squared test for model fit
cox.zph(cox1)

print(km1, rmean = "individual")
print(update(km1, . ~ generation), rmean = "individual")
print(update(km1, . ~ treatment), rmean = "individual")

# Fit quasi-Poisson GLM with viable eggs as response for individuals with positive egg count
# Block (intercept) and covariate generation (simple linear model) taking into account generation
total2 <- subset(total, total_eggs > 0)
fit1 <- glm(cbind(viable, total_eggs - viable) ~ block + generation * treatment,
            family = quasibinomial, data = total2)

hnp(fit1)

anova(fit1, test = "F")

round(cbind(with(total2, tapply(viable/total_eggs*100, generation, mean)),
            with(total2, tapply(viable/total_eggs*100, generation, se))), 4)

round(cbind(with(total2, tapply(viable/total_eggs*100, treatment, mean)),
            with(total2, tapply(viable/total_eggs*100, treatment, se))), 4)

####################################
# Longitudinal and joint modelling #
####################################

eggs <-  read.csv("eggs.csv")

eggs$nf <- ifelse(eggs$generation == "F1", 10, 5)
eggs$block <- as.factor(eggs$block)
eggs$replicate <- as.factor(eggs$replicate)
eggs$generation <- as.factor(eggs$generation)
eggs$treatment <- as.factor(eggs$treatment)

eggs$id <- as.factor(as.character(with(eggs, generation:treatment:block:replicate)))
levels(eggs$id) <- 1:36
survival <- with(eggs, tapply(day, id, max))
survival[is.na(survival)] <- 0
eggs$survival <- unlist(sapply(survival, function(x) replicate(x, x)))
eggs.id <- eggs[cumsum(survival),]

# Longitudinal model
fit1 <- lmer(log(eggs_per_day + 1) ~ block + generation * treatment * poly(day, 2) +
               (1 | id), data = eggs)
hnp(fit1, verb = TRUE)

# Joint model
eggs.id$survival <- eggs.id$survival + .1
fitLME <- lme(log(eggs_per_day + 1) ~ block + generation * treatment * day +
                offset(log(nf)),
              random = ~ 1 | id, data = eggs)
fitSURV <- coxph(Surv(survival) ~ block + generation * treatment,
                 data = eggs.id, x = TRUE)

fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "day",
                       method = "Cox-PH-GH")
fitLME2 <- lme(log(eggs_per_day + 1) ~ block + generation * treatment + day +
                 offset(log(nf)),
               random = ~ 1 | id, data = eggs)
fitJOINT2 <- jointModel(fitLME2, fitSURV, timeVar = "day",
                        method = "Cox-PH-GH")
anova(fitJOINT2, fitJOINT)

fitLME3 <- lme(log(eggs_per_day + 1) ~ block + generation + treatment + day +
                 offset(log(nf)),
               random = ~ 1 | id, data = eggs)
fitJOINT3 <- jointModel(fitLME3, fitSURV, timeVar = "day",
                        method = "Cox-PH-GH")
anova(fitJOINT3, fitJOINT2)

fitSURV2 <- coxph(Surv(survival) ~ block + generation + treatment + offset(log(nf)),
                  data = eggs.id, x = TRUE)
fitJOINT4 <- jointModel(fitLME2, fitSURV2, timeVar = "day",
                        method = "Cox-PH-GH")
anova(fitJOINT4, fitJOINT2)

fitSURV3 <- coxph(Surv(survival) ~ block + treatment + offset(log(nf)),
                  data = eggs.id, x = TRUE)

fitJOINT5 <- jointModel(fitLME2, fitSURV3, timeVar = "day",
                        method = "Cox-PH-GH")
anova(fitJOINT5, fitJOINT4)
summary(fitJOINT5)

#############################################
# Spaghetti plot of eggs per female vs days #
#############################################

eggs <- eggs[,c(9,1,3,4,5,6,7,8)]
eggs.id <- eggs.id[,c(9,1,3,4,5,6,7,8)]
eggs$treatment2 <- eggs$treatment
levels(eggs$treatment2) <- c("Raspberry","Strawberry")

ggplot(data = eggs, mapping = aes(x = day, y = eggs_per_day/nf, group = id)) + theme_bw() +
  geom_line(lwd = .3, col = colours[3]) +
  scale_x_continuous(name = "Time (days)") +
  scale_y_continuous(name = "Number of eggs laid per female") +
  facet_wrap(~ generation + treatment2) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "#BDBDBD", fill = NA, linetype = "dotted"),
        legend.spacing.y = unit(0, "mm"),
        aspect.ratio = 0.5, axis.text = element_text(colour = 1, size = 12),
        strip.text = element_text(face = "plain", color = "black",
                                  hjust = 0, size = 8),
        strip.background = element_rect(fill = "white", color = "#BDBDBD", linetype = "dotted"),
        legend.box.background = element_rect(colour = "#BDBDBD", fill = NA, linetype = "dotted"),
        legend.title = element_blank(),
        legend.text = element_text(face = "plain", hjust = 0, size = 10),
        legend.position = "none",
        legend.background = element_rect(fill = "white", color = "#BDBDBD", linetype = "dotted"))
