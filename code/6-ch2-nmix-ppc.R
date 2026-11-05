## =================================================
##
## Title: 6-ch2-nmix-ppc.R
##
## Author: Jasmine Williamson
## Date Created: 11/5/2026
##
## Description: Posterior predictive check for age-structured N-mixture model.
##   Simulates replicated datasets from the posterior and compares summary
##   statistics to observed data. Two checks:
##     1. Total captures by removal pass (tests removal trend)
##     2. Total captures by age class (tests age structure fit)
##
## =================================================


## settings ---------------------------------------------------

  rm(list = ls())
  setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/abundance-ch2")

  library(ggplot2)
  library(dplyr)

  dataset <- "oss"   # or "enes"
  n_sim   <- 500     # posterior draws to simulate from

  out_dir <- "data"
  fig_dir <- "figures"
  K       <- 3


## load MCMC samples ------------------------------------------

  a       <- readRDS(file.path(out_dir, paste0("mcmc_list_notrt_", dataset, ".rds")))
  samples <- do.call(rbind, lapply(a, as.matrix))
  draw_idx <- sample(1:nrow(samples), n_sim)


## load and prep data (mirrors model script sections 1-5) -----

  counts_file <- switch(dataset,
    oss  = "data/pass-level-counts-o.csv",
    enes = "data/pass-level-counts-e.csv"
  )
  counts.o    <- read.csv(counts_file,
                           colClasses = c(stand = "character", obs = "character"))
  counts.o$subplot <- as.integer(counts.o$subplot)
  counts.o$pass    <- as.integer(counts.o$pass)

  sals        <- read.csv("data/sals.complete.csv",
                           colClasses = c(stand = "character", subplot = "character",
                                          pass = "character", age_class = "factor",
                                          recap = "character", spp = "factor"))
  subplot.dat <- read.csv("data/subplot.complete.new.csv",
                           colClasses = c(stand = "character", subplot = "integer"))
  env2        <- read.csv("data/env_subset_corr2.csv")

  sites <- counts.o %>%
    distinct(site_id, stand, trt, year, subplot) %>%
    arrange(site_id, subplot) %>%
    mutate(site_idx = row_number())
  I <- nrow(sites)

  year_levels  <- sort(unique(sites$year))
  stand_levels <- sort(unique(sites$stand))
  year         <- match(sites$year,  year_levels)
  stand        <- match(sites$stand, stand_levels)
  nstand       <- length(stand_levels)


## build observed y array [I, K, 3] ---------------------------

  spp_code   <- switch(dataset, oss = "OSS", enes = "ENES")
  age_levels <- c("J", "SA", "A")

  sals_spp <- sals %>%
    filter(spp == spp_code, recap != "1", age_class %in% age_levels) %>%
    mutate(subplot = as.integer(subplot), pass = as.integer(pass))

  sal_counts <- sals_spp %>%
    group_by(site_id, subplot, pass, age_class) %>%
    summarise(count = n(), .groups = "drop")

  y <- array(0L, dim = c(I, K, 3), dimnames = list(NULL, NULL, age_levels))

  for (r in seq_len(nrow(sal_counts))) {
    si <- which(sites$site_id == sal_counts$site_id[r] &
                sites$subplot  == sal_counts$subplot[r])
    k  <- sal_counts$pass[r]
    a  <- match(as.character(sal_counts$age_class[r]), age_levels)
    if (length(si) == 1 && !is.na(k) && k >= 1 && k <= K && !is.na(a))
      y[si, k, a] <- sal_counts$count[r]
  }

  cat("Total captures (known age):", sum(y), "\n")
  cat("Captures by pass:", apply(y, 2, sum), "\n")
  cat("Captures by age:", apply(y, 3, sum), "\n")


## covariates -------------------------------------------------

  site_covs <- sites %>%
    left_join(subplot.dat %>% select(site_id, subplot, canopy_cov, dwd_count, dwd_cov,
                                     decay_cl, char_cl, temp, soil_moist_avg),
              by = c("site_id", "subplot")) %>%
    left_join(env2 %>% select(site_id, days_since_rain, avg_volume), by = "site_id")

  zs <- function(x) as.numeric(scale(x))

  canopy_cov      <- zs(site_covs$canopy_cov)
  dwd_count       <- zs(site_covs$dwd_count)
  dwd_cov         <- zs(site_covs$dwd_cov)
  decay_cl        <- zs(site_covs$decay_cl);    decay_cl[is.na(decay_cl)] <- 0
  char_cl         <- zs(site_covs$char_cl);     char_cl[is.na(char_cl)]   <- 0
  avg_volume      <- zs(site_covs$avg_volume)
  temp_raw        <- site_covs$temp
  temp            <- zs(temp_raw)
  temp_sq         <- zs(temp_raw^2)
  soil_moist      <- zs(site_covs$soil_moist_avg)
  days_since_rain <- zs(site_covs$days_since_rain)


## back-calculate unmonitored parameters ----------------------

  # mu.p and eps.p from p_age_baseline (eps.p[1] = 0 by constraint)
  mu_p_samp  <- qlogis(samples[draw_idx, "p_age_baseline[1]"])
  eps_p2_samp <- qlogis(samples[draw_idx, "p_age_baseline[2]"]) - mu_p_samp
  eps_p3_samp <- qlogis(samples[draw_idx, "p_age_baseline[3]"]) - mu_p_samp

  # phi0[2] and phi0[3] from pi_age_baseline (phi0[1] = 0 by constraint)
  phi0_2_samp <- log(samples[draw_idx, "pi_age_baseline[2]"] /
                     samples[draw_idx, "pi_age_baseline[1]"])
  phi0_3_samp <- log(samples[draw_idx, "pi_age_baseline[3]"] /
                     samples[draw_idx, "pi_age_baseline[1]"])


## simulate replicated datasets --------------------------------

  sim_by_pass <- matrix(NA_real_, n_sim, K)
  sim_by_age  <- matrix(NA_real_, n_sim, 3)

  for (s in 1:n_sim) {
    d <- draw_idx[s]

    # abundance
    b0       <- samples[d, "beta0"]
    b.canopy <- samples[d, "beta.canopy"]
    b.dwd    <- samples[d, "beta.dwd.count"]
    b.decay  <- samples[d, "beta.decay"]
    b.char   <- samples[d, "beta.char"]
    b.vol    <- samples[d, "beta.vol"]
    b.year2  <- samples[d, "beta.year[2]"]

    # age composition
    g.dwdcov <- c(0, samples[d, "gam.dwdcov[2]"], samples[d, "gam.dwdcov[3]"])
    g.decay  <- c(0, samples[d, "gam.decay[2]"],  samples[d, "gam.decay[3]"])
    g.soil   <- c(0, samples[d, "gam.soil[2]"],   samples[d, "gam.soil[3]"])

    # detection (observer RE set to 0 — averages to zero across observers)
    mu_p   <- mu_p_samp[s]
    eps_p  <- c(0, eps_p2_samp[s], eps_p3_samp[s])
    b.temp  <- samples[d, "beta.temp"]
    b.temp2 <- samples[d, "beta.temp2"]
    b.soil  <- samples[d, "beta.soil"]
    b.days  <- samples[d, "beta.days"]

    y_rep <- array(0L, dim = c(I, K, 3))

    # draw stand random effects from their estimated distribution
    sigma.stand   <- samples[d, "sigma.stand"]
    alpha.stand_s <- rnorm(nstand, 0, sigma.stand)

    for (i in 1:I) {

      # expected abundance
      yr_eff <- if (year[i] == 1) 0 else b.year2
      lam <- exp(b0 + alpha.stand_s[stand[i]] + yr_eff +
                 b.canopy*canopy_cov[i] + b.dwd*dwd_count[i] +
                 b.decay*decay_cl[i]    + b.char*char_cl[i] +
                 b.vol*avg_volume[i])

      # age composition via softmax
      phi  <- c(0,
                phi0_2_samp[s] + g.dwdcov[2]*dwd_cov[i] + g.decay[2]*decay_cl[i] + g.soil[2]*soil_moist[i],
                phi0_3_samp[s] + g.dwdcov[3]*dwd_cov[i] + g.decay[3]*decay_cl[i] + g.soil[3]*soil_moist[i])
      ep   <- exp(phi)
      pi   <- ep / sum(ep)

      N_age <- rpois(3, lam * pi)

      # removal sampling across passes
      p_det <- plogis(mu_p + eps_p +
                      b.temp*temp[i] + b.temp2*temp_sq[i] +
                      b.soil*soil_moist[i] + b.days*days_since_rain[i])

      for (a in 1:3) {
        avail <- N_age[a]
        for (k in 1:K) {
          caught        <- rbinom(1, avail, p_det[a])
          y_rep[i,k,a] <- caught
          avail         <- avail - caught
        }
      }
    }

    sim_by_pass[s, ] <- apply(y_rep, 2, sum)
    sim_by_age[s, ]  <- apply(y_rep, 3, sum)
  }


## observed statistics ----------------------------------------

  obs_by_pass <- apply(y, 2, sum)
  obs_by_age  <- apply(y, 3, sum)


## plots ------------------------------------------------------

  ## 1. captures by removal pass
  pass_long <- data.frame(
    pass = factor(rep(1:K, each = n_sim)),
    sim  = c(sim_by_pass)
  )
  obs_pass <- data.frame(pass = factor(1:K), obs = obs_by_pass)

  p_pass <- ggplot(pass_long, aes(x = pass, y = sim)) +
    geom_boxplot(fill = "gray80", outlier.size = 0.5, width = 0.5) +
    geom_point(data = obs_pass, aes(x = pass, y = obs),
               color = "firebrick", size = 3, shape = 18) +
    labs(
      title   = paste("PPC: Captures by Removal Pass —", toupper(dataset)),
      subtitle = "Red = observed; boxes = posterior predictive",
      x       = "Pass",
      y       = "Total captures"
    ) +
    theme_classic() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13),
          plot.title = element_text(face = "bold", size = 13))

  p_pass

  ggsave(file.path(fig_dir, paste0("ppc-by-pass-", dataset, ".png")),
         p_pass, width = 4, height = 4, dpi = 300)


  ## 2. captures by age class
  age_long <- data.frame(
    age = factor(rep(c("J","SA","A"), each = n_sim), levels = c("J","SA","A")),
    sim = c(sim_by_age)
  )
  obs_age <- data.frame(
    age = factor(c("J","SA","A"), levels = c("J","SA","A")),
    obs = obs_by_age
  )

  p_age <- ggplot(age_long, aes(x = age, y = sim)) +
    geom_boxplot(fill = "gray80", outlier.size = 0.5, width = 0.5) +
    geom_point(data = obs_age, aes(x = age, y = obs),
               color = "firebrick", size = 3, shape = 18) +
    scale_x_discrete(labels = c("J" = "Juvenile", "SA" = "Sub-adult", "A" = "Adult")) +
    labs(
      title    = paste("PPC: Captures by Age Class —", toupper(dataset)),
      subtitle = "Red = observed; boxes = posterior predictive",
      x        = NULL,
      y        = "Total captures"
    ) +
    theme_classic() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13),
          plot.title = element_text(face = "bold", size = 13))

  p_age

  ggsave(file.path(fig_dir, paste0("ppc-by-age-", dataset, ".png")),
         p_age, width = 5, height = 4, dpi = 300)


  cat("PPC complete for", toupper(dataset), "\n")
  cat("Observed by pass:", obs_by_pass, "\n")
  cat("Sim median by pass:", round(apply(sim_by_pass, 2, median)), "\n")
  cat("Observed by age:", obs_by_age, "\n")
  cat("Sim median by age:", round(apply(sim_by_age, 2, median)), "\n")
