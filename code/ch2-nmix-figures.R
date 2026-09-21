## =================================================
##
## Title: ch2-nmix-figures.R
##
## Author: Jasmine Williamson
## Date Created: 9/21/2026
##
## Description: Results and figures from n-mixture model for ch2.
##
## Results from script nmix-model-oss-ch2.R.
##
##
## =================================================


## settings ---------------------------------------------------

  rm(list=ls())
  setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git")
  
  library(ggplot2)
  library(dplyr)

  ## load saved results
  out_dir       <- "/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/abundance-ch2/data"
  dataset       <- "enes"   # or "enes"
  a             <- readRDS(file.path(out_dir, paste0("mcmc_list_",     dataset, ".rds")))
  chain_samples <- readRDS(file.path(out_dir, paste0("chain_samples_", dataset, ".rds")))
  samples       <- do.call(rbind, lapply(a, as.matrix))

  fig_dir <- "/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/abundance-ch2/figures"


##### Posterior Coefficient Estimates ----------------------------------------

  params_to_plot <- c(
    # Abundance (lambda)
    "beta.trt[2]", "beta.trt[3]", "beta.trt[4]", "beta.trt[5]",
    "beta.year[2]",
    "beta.canopy", "beta.dwd.count", "beta.decay", "beta.char",
    "beta.fwd", "beta.veg", "beta.vol",
    # Age composition (pi_age)
    "gam.canopy[2]", "gam.canopy[3]",
    "gam.dwdcov[2]", "gam.dwdcov[3]",
    "gam.soil[2]",   "gam.soil[3]",
    "gam.fwd[2]",    "gam.fwd[3]",
    # Detection (p)
    "beta.temp", "beta.temp2", "beta.soil", "beta.days", "beta.jul"
  )

  coef_df <- data.frame(
    Parameter = params_to_plot,
    Mean      = colMeans(samples[, params_to_plot]),
    LCI       = apply(samples[, params_to_plot], 2, quantile, 0.025),
    UCI       = apply(samples[, params_to_plot], 2, quantile, 0.975),
    row.names = NULL
  )

  coef_df <- coef_df %>%
    mutate(
      Submodel = case_when(
        Parameter %in% c("beta.trt[2]", "beta.trt[3]", "beta.trt[4]", "beta.trt[5]",
                          "beta.year[2]", "beta.canopy", "beta.dwd.count", "beta.decay",
                          "beta.char", "beta.fwd", "beta.veg", "beta.vol")
                                           ~ "Abundance (λ)",
        Parameter %in% c("gam.canopy[2]", "gam.canopy[3]", "gam.dwdcov[2]", "gam.dwdcov[3]",
                          "gam.soil[2]",   "gam.soil[3]",   "gam.fwd[2]",    "gam.fwd[3]")
                                           ~ "Age Composition (π)",
        Parameter %in% c("beta.temp", "beta.temp2", "beta.soil", "beta.days", "beta.jul")
                                           ~ "Detection (p)"
      ),
      Submodel = factor(Submodel, levels = c("Abundance (λ)",
                                              "Age Composition (π)",
                                              "Detection (p)")),
      Parameter = dplyr::recode(Parameter,
        "beta.trt[2]"    = "Burn Salvage (BS)",
        "beta.trt[3]"    = "Burn (BU)",
        "beta.trt[4]"    = "Harvest Burn (HB)",
        "beta.trt[5]"    = "Harvest (HU)",
        "beta.year[2]"   = "Year (2024)",
        "beta.canopy"    = "Canopy Cover",
        "beta.dwd.count" = "DWD Count",
        "beta.decay"     = "Decay Class",
        "beta.char"      = "Char Class",
        "beta.fwd"       = "FWD Cover",
        "beta.veg"       = "Veg Cover",
        "beta.vol"       = "DWD Volume",
        "gam.canopy[2]"  = "Canopy Cover (SA)",
        "gam.canopy[3]"  = "Canopy Cover (A)",
        "gam.dwdcov[2]"  = "DWD Cover (SA)",
        "gam.dwdcov[3]"  = "DWD Cover (A)",
        "gam.soil[2]"    = "Soil Moisture (SA)",
        "gam.soil[3]"    = "Soil Moisture (A)",
        "gam.fwd[2]"     = "FWD Cover (SA)",
        "gam.fwd[3]"     = "FWD Cover (A)",
        "beta.temp"      = "Air Temp (linear)",
        "beta.temp2"     = "Air Temp (quadratic)",
        "beta.soil"      = "Soil Moisture",
        "beta.days"      = "Days Since Rain",
        "beta.jul"       = "Julian Date"
      )
    )

  p_coef <- ggplot(coef_df, aes(x = Mean, y = reorder(Parameter, Mean))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2) +
    geom_point(size = 2) +
    facet_wrap(~ Submodel, scales = "free_y", ncol = 1) +
    labs(
      title = paste("Posterior Coefficient Estimates —", toupper(dataset)),
      x     = "Posterior Mean (95% CI)",
      y     = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12),
      axis.title.x = element_text(size = 13),
      strip.text   = element_text(size = 13, face = "bold"),
      plot.title   = element_text(face = "bold", size = 13)
    )

  p_coef

  ggsave(
    filename = file.path(fig_dir, paste0("coeff-plot-", dataset, ".png")),
    plot     = p_coef,
    width    = 7,
    height   = 10,
    dpi      = 300
  )

  cat("Saved:", file.path(fig_dir, paste0("coeff-plot-", dataset, ".png")), "\n")


##### Age Composition by Treatment (pi_age_baseline) -------------------------

  trt_labels  <- c("UU", "BS", "BU", "HB", "HU")
  age_labels  <- c("J", "SA", "A")

  age_rows <- lapply(1:5, function(t) {
    lapply(1:3, function(a) {
      col  <- paste0("pi_age_baseline[", t, ", ", a, "]")
      samp <- samples[, col]
      data.frame(
        treatment = trt_labels[t],
        age_class = age_labels[a],
        Mean      = mean(samp),
        LCI       = quantile(samp, 0.025),
        UCI       = quantile(samp, 0.975)
      )
    })
  })
  age_df <- do.call(rbind, do.call(c, age_rows))
  rownames(age_df) <- NULL

  age_df$treatment <- factor(age_df$treatment, levels = c("UU","BU","HU","HB","BS"))
  age_df$age_class <- factor(age_df$age_class, levels = c("J","SA","A"))

  p_age <- ggplot(age_df, aes(x = treatment, y = Mean, color = age_class,
                               group = age_class)) +
    geom_point(position = position_dodge(0.4), size = 3) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI),
                  width = 0.15, linewidth = 0.7,
                  position = position_dodge(0.4)) +
    scale_color_manual(
      values = c("J" = "#F4A261", "SA" = "#2A9D8F", "A" = "#264653"),
      labels = c("J" = "Juvenile", "SA" = "Sub-adult", "A" = "Adult")
    ) +
    scale_x_discrete(labels = c(
      "UU" = "Control", "BU" = "Burn-only",
      "HU" = "Harvest", "HB" = "Harvest-Burn", "BS" = "Burn-Salvage"
    )) +
    labs(
      title  = paste("Age Composition by Treatment —", toupper(dataset)),
      y      = "Proportion",
      x      = NULL,
      color  = "Age Class"
    ) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12),
      axis.title.y = element_text(size = 13),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 12),
      plot.title   = element_text(face = "bold", size = 13)
    )

  p_age

  ggsave(
    filename = file.path(fig_dir, paste0("age-comp-trt-", dataset, ".png")),
    plot     = p_age,
    width    = 7,
    height   = 4.5,
    dpi      = 300
  )

  cat("Saved:", file.path(fig_dir, paste0("age-comp-trt-", dataset, ".png")), "\n")


##### Detection Probability by Age Class (p_age_baseline) --------------------

  det_df <- data.frame(
    age_class = c("J", "SA", "A"),
    Mean      = sapply(1:3, function(a) mean(samples[, paste0("p_age_baseline[", a, "]")])),
    LCI       = sapply(1:3, function(a) quantile(samples[, paste0("p_age_baseline[", a, "]")], 0.025)),
    UCI       = sapply(1:3, function(a) quantile(samples[, paste0("p_age_baseline[", a, "]")], 0.975))
  )
  det_df$age_class <- factor(det_df$age_class, levels = c("J", "SA", "A"))

  p_det <- ggplot(det_df, aes(x = age_class, y = Mean)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.1, linewidth = 0.8) +
    scale_x_discrete(labels = c("J" = "Juvenile", "SA" = "Sub-adult", "A" = "Adult")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = paste("Baseline Detection Probability by Age Class —", toupper(dataset)),
      y     = "Detection probability",
      x     = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(size = 13),
      axis.text.y  = element_text(size = 12),
      axis.title.y = element_text(size = 13),
      plot.title   = element_text(face = "bold", size = 13)
    )

  p_det

  ggsave(
    filename = file.path(fig_dir, paste0("det-prob-age-", dataset, ".png")),
    plot     = p_det,
    width    = 4,
    height   = 4,
    dpi      = 300
  )

  cat("Saved:", file.path(fig_dir, paste0("det-prob-age-", dataset, ".png")), "\n")

