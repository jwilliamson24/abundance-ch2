## =================================================
##
## Title: 8-ch2-nmix-trt-lambda.R
##
## Author: Jasmine Williamson
## Date Created: 11/12/2026
##
## Description: Post-hoc treatment abundance comparison using Josh's approach:
##   For each treatment, compute the mean value of each habitat covariate
##   across all sites within that treatment (z-scored, same as model).
##   Then, for each posterior sample, predict lambda at those treatment-mean
##   covariate values and draw N ~ Poisson(lambda).
##   This gives expected abundance at a "typical site" in each treatment,
##   controlling for habitat composition within treatments.
##
## =================================================


## settings ---------------------------------------------------

  rm(list = ls())
  setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/abundance-ch2")

  library(ggplot2)
  library(dplyr)

  out_dir <- "data"
  fig_dir <- "figures"

  treatment_levels <- c("UU", "BS", "BU", "HB", "HU")
  n_sim <- 1000   # posterior draws to use


## load site covariate data -----------------------------------
## (mirrors model script section 5 — same z-scoring)

  counts.o    <- read.csv("data/pass-level-counts-o.csv",
                           colClasses = c(stand = "character", obs = "character"))
  counts.o$subplot <- as.integer(counts.o$subplot)

  subplot.dat <- read.csv("data/subplot.complete.new.csv",
                           colClasses = c(stand = "character", subplot = "integer"))
  env2        <- read.csv("data/env_subset_corr2.csv")

  sites <- counts.o %>%
    distinct(site_id, stand, trt, year, subplot) %>%
    arrange(site_id, subplot) %>%
    mutate(site_idx = row_number())

  site_covs <- sites %>%
    left_join(subplot.dat %>% select(site_id, subplot, canopy_cov, dwd_count,
                                     decay_cl, char_cl),
              by = c("site_id", "subplot")) %>%
    left_join(env2 %>% select(site_id, avg_volume), by = "site_id")

  zs <- function(x) as.numeric(scale(x))

  site_covs$canopy_z   <- zs(site_covs$canopy_cov)
  site_covs$dwd_z      <- zs(site_covs$dwd_count)
  site_covs$decay_z    <- zs(site_covs$decay_cl);  site_covs$decay_z[is.na(site_covs$decay_z)]  <- 0
  site_covs$char_z     <- zs(site_covs$char_cl);   site_covs$char_z[is.na(site_covs$char_z)]    <- 0
  site_covs$vol_z      <- zs(site_covs$avg_volume)


## treatment-mean z-scored covariates -------------------------

  trt_means <- site_covs %>%
    group_by(trt) %>%
    summarise(
      canopy = mean(canopy_z, na.rm = TRUE),
      dwd    = mean(dwd_z,    na.rm = TRUE),
      decay  = mean(decay_z,  na.rm = TRUE),
      char   = mean(char_z,   na.rm = TRUE),
      vol    = mean(vol_z,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(factor(trt, levels = treatment_levels))

  cat("Treatment-mean z-scored covariates:\n")
  print(trt_means)


## predict lambda for each treatment and species --------------

  predict_trt_lambda <- function(samples, species_label) {

    draw_idx <- sample(1:nrow(samples), n_sim)

    results <- lapply(1:5, function(t) {
      trt_name <- treatment_levels[t]
      tm       <- trt_means[trt_means$trt == trt_name, ]

      lam_draws <- sapply(draw_idx, function(d) {
        exp(
          samples[d, "beta0"] +
          samples[d, "beta.canopy"]    * tm$canopy +
          samples[d, "beta.dwd.count"] * tm$dwd    +
          samples[d, "beta.decay"]     * tm$decay  +
          samples[d, "beta.char"]      * tm$char   +
          samples[d, "beta.vol"]       * tm$vol
        )
      })

      data.frame(
        treatment = trt_name,
        Mean      = mean(lam_draws),
        Median    = median(lam_draws),
        LCI       = quantile(lam_draws, 0.025),
        UCI       = quantile(lam_draws, 0.975),
        species   = species_label,
        row.names = NULL
      )
    })

    do.call(rbind, results)
  }


## load samples and run ---------------------------------------

  samples_oss  <- do.call(rbind, lapply(
    readRDS(file.path(out_dir, "mcmc_list_notrt_oss.rds")), as.matrix))
  samples_enes <- do.call(rbind, lapply(
    readRDS(file.path(out_dir, "mcmc_list_notrt_enes.rds")), as.matrix))

  preds_oss  <- predict_trt_lambda(samples_oss,  "OSS")
  preds_enes <- predict_trt_lambda(samples_enes, "ENES")

  both_preds <- rbind(preds_oss, preds_enes)
  both_preds$treatment <- factor(both_preds$treatment,
                                  levels = c("UU", "BU", "HU", "HB", "BS"))

  cat("\nPredicted N per site by treatment:\n")
  print(both_preds[, c("treatment","species","Median","LCI","UCI")])


## plot -------------------------------------------------------

  p_trt <- ggplot(both_preds, aes(x = treatment, y = Median, color = species)) +
    geom_point(position = position_dodge(0.5), size = 4) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI),
                  width = 0.1, linewidth = 0.8,
                  position = position_dodge(0.5)) +
    scale_color_manual(
      values = c("ENES" = "#6091C2", "OSS" = "gray20"),
      labels = c("ENES" = "Ensatina", "OSS" = "Oregon Slender Salamander")
    ) +
    scale_x_discrete(labels = c(
      "UU" = "Control",
      "BU" = "Burn-only",
      "HU" = "Harvest",
      "HB" = "Harvest-Burn",
      "BS" = "Burn-Salvage"
    )) +
    labs(
      title = "Predicted Abundance by Treatment (at Treatment-Mean Habitat)",
      y     = "Expected abundance per subplot (λ)",
      x     = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x            = element_text(size = 15),
      axis.text.y            = element_text(size = 16),
      axis.title.y           = element_text(size = 16),
      legend.title           = element_blank(),
      legend.text            = element_text(size = 14),
      legend.position        = "inside",
      legend.position.inside = c(0.75, 0.85),
      panel.grid.major.y     = element_blank(),
      panel.grid.major.x     = element_blank(),
      panel.border           = element_rect(color = "gray40", fill = NA, linewidth = 0.5),
      axis.ticks             = element_line(color = "gray40", linewidth = 0.5),
      axis.ticks.length      = unit(0.1, "cm"),
      plot.title             = element_text(hjust = 0.5, face = "bold", size = 14)
    )

  p_trt

  ggsave(
    filename = file.path(fig_dir, "both-spp-trt-lambda.png"),
    plot     = p_trt,
    width    = 8,
    height   = 5,
    bg       = "white",
    dpi      = 300
  )

  cat("Saved: both-spp-trt-lambda.png\n")
