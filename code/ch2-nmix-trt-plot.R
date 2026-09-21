## =================================================
##
## Title: ch2-nmix-trt-plot.R
##
## Author: Jasmine Williamson
## Date Created: 9/21/2026
##
## Description: Treatment effects on abundance — posterior Ntrt[t] by species.
##              Analogous to 4-both-occu-preds-treatment-plot.R from Ch. 1.
##              Ntrt[t] is the total N summed across all sites within treatment t.
##
## Results from nmix-model-oss-ch2.R.
##
## =================================================

## settings ---------------------------------------------------

  rm(list=ls())
  setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git")

  library(ggplot2)
  library(dplyr)

  out_dir <- "abundance-ch2/data"
  fig_dir <- "abundance-ch2/figures"

  # treatment order matches NIMBLE model: t = 1..5
  treatment_levels <- c("UU", "BS", "BU", "HB", "HU")


##### Site counts per treatment (same for both species) -----------------------

  site_csv  <- read.csv(file.path(out_dir, "pass-level-counts-o.csv"))
  n_per_trt <- site_csv %>%
    distinct(site_id, trt) %>%
    count(trt) %>%
    arrange(factor(trt, levels = treatment_levels)) %>%
    pull(n)   # order: UU BS BU HB HU = t 1..5


##### Helper: extract mean N per site from Ntrt posteriors -------------------

  extract_ntrt <- function(samples, species) {
    ntrt_cols <- paste0("Ntrt[", 1:5, "]")
    # divide each posterior draw by n_sites to get mean N per site
    ntrt_mat  <- samples[, ntrt_cols] / rep(n_per_trt, each = nrow(samples))
    means <- colMeans(ntrt_mat)
    cis   <- apply(ntrt_mat, 2, quantile, c(0.025, 0.975))
    data.frame(
      treatment = treatment_levels,
      Mean      = means,
      LCI       = cis[1, ],
      UCI       = cis[2, ],
      species   = species,
      row.names = NULL
    )
  }


##### Load posteriors for both species ----------------------------------------

  a_oss        <- readRDS(file.path(out_dir, "mcmc_list_oss.rds"))
  samples_oss  <- do.call(rbind, lapply(a_oss, as.matrix))

  a_enes       <- readRDS(file.path(out_dir, "mcmc_list_enes.rds"))
  samples_enes <- do.call(rbind, lapply(a_enes, as.matrix))

  oss_preds  <- extract_ntrt(samples_oss,  "OSS")
  enes_preds <- extract_ntrt(samples_enes, "ENES")

  both_preds <- rbind(oss_preds, enes_preds)
  both_preds$treatment <- factor(both_preds$treatment,
                                  levels = c("UU", "BU", "HU", "HB", "BS"))

  both_preds[ , c("Mean","LCI","UCI")] <- round(both_preds[ , c("Mean","LCI","UCI")], 1)
  print(both_preds)


##### Combined plot -----------------------------------------------------------

  p_trt <- ggplot(both_preds, aes(x = treatment, y = Mean, color = species)) +
    geom_point(position = position_dodge(0.5), size = 4) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI),
                  width = 0.1, linewidth = 0.8,
                  position = position_dodge(0.5)) +
    scale_color_manual(
      values = c("ENES" = "#6091C2", "OSS" = "gray20"),
      labels = c("ENES" = "Ensatina",
                 "OSS"  = "Oregon Slender Salamander")
    ) +
    scale_x_discrete(labels = c(
      "UU" = "Control",
      "BU" = "Burn-only",
      "HU" = "Harvest",
      "HB" = "Harvest-Burn",
      "BS" = "Burn-Salvage"
    )) +
    labs(
      title = "Estimated Total Abundance by Treatment and Species",
      y     = "Estimated N per subplot",
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
      strip.text             = element_text(size = 15, face = "bold"),
      plot.title             = element_text(hjust = 0.5, face = "bold", size = 16),
      panel.grid.major.y     = element_blank(),
      panel.grid.major.x     = element_blank(),
      panel.border           = element_rect(color = "gray40", fill = NA, linewidth = 0.5),
      axis.ticks             = element_line(color = "gray40", linewidth = 0.5),
      axis.ticks.length      = unit(0.1, "cm")
    )

  p_trt

  ggsave(
    filename = file.path(fig_dir, "both-spp-trt-abundance.png"),
    plot     = p_trt,
    width    = 8,
    height   = 5,
    bg       = "white",
    dpi      = 500
  )

  cat("Saved:", file.path(fig_dir, "both-spp-trt-abundance.png"), "\n")
