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
  dataset       <- "oss"   # or "enes"
  a             <- readRDS(file.path(out_dir, paste0("mcmc_list_notrt_",     dataset, ".rds")))
  chain_samples <- readRDS(file.path(out_dir, paste0("chain_samples_notrt_", dataset, ".rds")))
  samples       <- do.call(rbind, lapply(a, as.matrix))

  fig_dir <- "/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/abundance-ch2/figures"


##### Posterior Coefficient Estimates ----------------------------------------

  params_to_plot <- c(
    # Abundance (lambda)
    "beta.year[2]",
    "beta.canopy", "beta.dwd.count", "beta.decay", "beta.char",
    "beta.vol",
    # Age composition (pi_age)
    "gam.dwdcov[2]", "gam.dwdcov[3]",
    "gam.decay[2]",  "gam.decay[3]",
    "gam.soil[2]",   "gam.soil[3]",
    # Detection (p)
    "beta.temp", "beta.temp2", "beta.soil", "beta.days"
  )

  coef_df <- data.frame(
    Parameter = params_to_plot,
    Mean      = colMeans(samples[, params_to_plot]),
    LCI       = apply(samples[, params_to_plot], 2, quantile, 0.025),
    UCI       = apply(samples[, params_to_plot], 2, quantile, 0.975),
    LCI75     = apply(samples[, params_to_plot], 2, quantile, 0.125),
    UCI75     = apply(samples[, params_to_plot], 2, quantile, 0.875),
    row.names = NULL
  )

  abund_params <- c("beta.year[2]", "beta.canopy", "beta.dwd.count",
                    "beta.decay", "beta.char", "beta.vol")
  age_params   <- c("gam.dwdcov[2]", "gam.dwdcov[3]",
                    "gam.decay[2]",  "gam.decay[3]",
                    "gam.soil[2]",   "gam.soil[3]")
  det_params   <- c("beta.temp", "beta.temp2", "beta.soil", "beta.days")

  coef_df <- coef_df %>%
    mutate(
      Submodel = case_when(
        Parameter %in% abund_params ~ "Abundance (λ)",
        Parameter %in% age_params   ~ "Age Composition (π)",
        Parameter %in% det_params   ~ "Detection (p)"
      ),
      Submodel = factor(Submodel, levels = c("Abundance (λ)",
                                              "Age Composition (π)",
                                              "Detection (p)")),
      Parameter = dplyr::recode(Parameter,
        "beta.year[2]"   = "Year (2024)",
        "beta.canopy"    = "Canopy Cover",
        "beta.dwd.count" = "DWD Count",
        "beta.decay"     = "Decay Class",
        "beta.char"      = "Char Class",
        "beta.vol"       = "DWD Volume",
        "gam.dwdcov[2]"  = "DWD Cover (SA)",
        "gam.dwdcov[3]"  = "DWD Cover (A)",
        "gam.decay[2]"   = "Decay Class (SA)",
        "gam.decay[3]"   = "Decay Class (A)",
        "gam.soil[2]"    = "Soil Moisture (SA)",
        "gam.soil[3]"    = "Soil Moisture (A)",
        "beta.temp"      = "Air Temp (linear)",
        "beta.temp2"     = "Air Temp (quadratic)",
        "beta.soil"      = "Soil Moisture (det)",
        "beta.days"         = "Days Since Rain"
      )
    )

  # abundance + detection: sort by posterior mean (most negative at bottom)
  abund_order <- coef_df %>% filter(Submodel == "Abundance (λ)") %>% arrange(Mean) %>% pull(Parameter)
  det_order   <- coef_df %>% filter(Submodel == "Detection (p)") %>% arrange(Mean) %>% pull(Parameter)

  # age comp: map each SA/A label back to its covariate group
  age_group_map <- c(
    "DWD Cover (SA)"    = "DWD Cover",     "DWD Cover (A)"     = "DWD Cover",
    "Decay Class (SA)"  = "Decay Class",   "Decay Class (A)"   = "Decay Class",
    "Soil Moisture (SA)"= "Soil Moisture", "Soil Moisture (A)" = "Soil Moisture"
  )
  age_sub <- coef_df %>% filter(Submodel == "Age Composition (π)")
  age_sub$group <- age_group_map[age_sub$Parameter]
  # sort groups by their average mean, then list SA below A within each group
  group_order <- age_sub %>%
    group_by(group) %>% summarise(gmean = mean(Mean), .groups = "drop") %>%
    arrange(gmean) %>% pull(group)
  age_order <- unlist(lapply(group_order, function(g) c(paste0(g, " (SA)"), paste0(g, " (A)"))))

  # set factor levels across all panels (facets use global level order)
  coef_df$Parameter <- factor(coef_df$Parameter,
                               levels = c(det_order, age_order, abund_order))

  p_coef <- ggplot(coef_df, aes(x = Mean, y = Parameter)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = LCI,   xmax = UCI),   height = 0.2, linewidth = 0.6) +
    geom_errorbarh(aes(xmin = LCI75, xmax = UCI75), height = 0,   linewidth = 1.6) +
    geom_point(size = 2) +
    facet_wrap(~ Submodel, scales = "free_y", ncol = 1) +
    coord_cartesian(xlim = c(min(coef_df$LCI) - 0.2, max(coef_df$UCI) + 0.2)) +
    labs(
      title = paste("Posterior Coefficient Estimates —", toupper(dataset)),
      x     = "Posterior Mean (75% / 95% CI)",
      y     = "Parameter"
    ) +
    theme(
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      legend.position = "none",
      strip.text   = element_text(size = 13, face = "bold"),
      plot.title   = element_text(face = "bold", size = 13)
    )

  p_coef

  ggsave(
    filename = file.path(fig_dir, paste0("coeff-plot-notrt-", dataset, ".png")),
    plot     = p_coef,
    width    = 7,
    height   = 10,
    dpi      = 300
  )


  

##### Age Composition — Overall Baseline (pi_age_baseline) -------------------

  age_labels <- c("J", "SA", "A")

  age_base_df <- data.frame(
    age_class = age_labels,
    Mean      = sapply(1:3, function(a) mean(samples[, paste0("pi_age_baseline[", a, "]")])),
    LCI       = sapply(1:3, function(a) quantile(samples[, paste0("pi_age_baseline[", a, "]")], 0.025)),
    UCI       = sapply(1:3, function(a) quantile(samples[, paste0("pi_age_baseline[", a, "]")], 0.975))
  )
  age_base_df$age_class <- factor(age_base_df$age_class, levels = c("J", "SA", "A"))

  p_age_base <- ggplot(age_base_df, aes(x = age_class, y = Mean)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.1, linewidth = 0.8) +
    scale_x_discrete(labels = c("J" = "Juvenile", "SA" = "Sub-adult", "A" = "Adult")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = paste("Baseline Age Composition —", toupper(dataset)),
      y     = "Proportion",
      x     = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(size = 13),
      axis.text.y  = element_text(size = 12),
      axis.title.y = element_text(size = 13),
      plot.title   = element_text(face = "bold", size = 13)
    )

  p_age_base

  ggsave(
    filename = file.path(fig_dir, paste0("age-comp-baseline-notrt-", dataset, ".png")),
    plot     = p_age_base,
    width    = 4,
    height   = 4,
    dpi      = 300
  )

  cat("Saved:", file.path(fig_dir, paste0("age-comp-baseline-", dataset, ".png")), "\n")


##### Age Composition by Treatment (pi_age_baseline) — if treatment in model --
if (FALSE) {  # treatment not in current model; pi_age_baseline is 1D — skip this block

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

} # end if (FALSE) — age comp by treatment block



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


##### Marginal Effects: Predicted Abundance vs. Decay Class -------------------

  # Load original covariate data for axis back-transformation
  subplot.dat <- read.csv(file.path(out_dir, "subplot.complete.new.csv"))
  env2        <- read.csv(file.path(out_dir, "env_subset_corr2.csv"))

  decay_mean <- mean(subplot.dat$decay_cl, na.rm = TRUE)
  decay_sd   <- sd(subplot.dat$decay_cl,   na.rm = TRUE)
  dwd_mean   <- mean(subplot.dat$dwd_cov,  na.rm = TRUE)
  dwd_sd     <- sd(subplot.dat$dwd_cov,    na.rm = TRUE)

  # z-scored sequence spanning observed range
  decay_seq_z    <- seq(-2, 2, length.out = 100)
  decay_seq_orig <- decay_seq_z * decay_sd + decay_mean

  # predicted lambda at each decay value (all other covariates at their mean = 0)
  pred_lambda_mat <- sapply(decay_seq_z, function(dc) {
    exp(samples[, "beta0"] + samples[, "beta.decay"] * dc)
  })

  pred_abund_df <- data.frame(
    decay_orig = decay_seq_orig,
    Median     = apply(pred_lambda_mat, 2, median),
    LCI        = apply(pred_lambda_mat, 2, quantile, 0.025),
    UCI        = apply(pred_lambda_mat, 2, quantile, 0.975)
  )

  p_abund_pred <- ggplot(pred_abund_df, aes(x = decay_orig, y = Median)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), fill = "gray70", alpha = 0.5) +
    geom_line(linewidth = 1) +
    labs(
      title = paste("Predicted Abundance vs. Decay Class —", toupper(dataset)),
      x     = "Decay Class",
      y     = "Expected abundance (λ)"
    ) +
    theme_classic() +
    theme(
      axis.text  = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.title = element_text(face = "bold", size = 13)
    )

  p_abund_pred

  ggsave(
    filename = file.path(fig_dir, paste0("pred-abund-decay-", dataset, ".png")),
    plot     = p_abund_pred,
    width    = 5, height = 4, dpi = 300
  )


##### Marginal Effects: Predicted Age Composition vs. DWD Cover ---------------

  # Back-calculate phi0[2] and phi0[3] from monitored pi_age_baseline
  # (phi0[1] = 0 by definition; phi0[a] = log(pi_base[a] / pi_base[1]))
  pi1_base <- samples[, "pi_age_baseline[1]"]
  pi2_base <- samples[, "pi_age_baseline[2]"]
  pi3_base <- samples[, "pi_age_baseline[3]"]
  phi0_2   <- log(pi2_base / pi1_base)
  phi0_3   <- log(pi3_base / pi1_base)

  gam2 <- samples[, "gam.dwdcov[2]"]
  gam3 <- samples[, "gam.dwdcov[3]"]

  dwd_seq_z    <- seq(-2, 2, length.out = 100)
  dwd_seq_orig <- dwd_seq_z * dwd_sd + dwd_mean

  age_pred_list <- lapply(seq_along(dwd_seq_z), function(j) {
    dw   <- dwd_seq_z[j]
    phi1 <- 0
    phi2 <- phi0_2 + gam2 * dw
    phi3 <- phi0_3 + gam3 * dw
    denom <- exp(phi1) + exp(phi2) + exp(phi3)
    pi1 <- exp(phi1) / denom
    pi2 <- exp(phi2) / denom
    pi3 <- exp(phi3) / denom
    data.frame(
      dwd_orig = dwd_seq_orig[j],
      J_med    = median(pi1),
      J_lci    = as.numeric(quantile(pi1, 0.025)),
      J_uci    = as.numeric(quantile(pi1, 0.975)),
      SA_med   = median(pi2),
      SA_lci   = as.numeric(quantile(pi2, 0.025)),
      SA_uci   = as.numeric(quantile(pi2, 0.975)),
      A_med    = median(pi3),
      A_lci    = as.numeric(quantile(pi3, 0.025)),
      A_uci    = as.numeric(quantile(pi3, 0.975))
    )
  })
  age_pred_wide <- do.call(rbind, age_pred_list)

  age_pred_long <- rbind(
    data.frame(dwd_orig = age_pred_wide$dwd_orig, age_class = "J",
               Median = age_pred_wide$J_med, LCI = age_pred_wide$J_lci, UCI = age_pred_wide$J_uci),
    data.frame(dwd_orig = age_pred_wide$dwd_orig, age_class = "SA",
               Median = age_pred_wide$SA_med, LCI = age_pred_wide$SA_lci, UCI = age_pred_wide$SA_uci),
    data.frame(dwd_orig = age_pred_wide$dwd_orig, age_class = "A",
               Median = age_pred_wide$A_med, LCI = age_pred_wide$A_lci, UCI = age_pred_wide$A_uci)
  )
  age_pred_long$age_class <- factor(age_pred_long$age_class, levels = c("J", "SA", "A"))

  p_age_pred <- ggplot(age_pred_long,
                       aes(x = dwd_orig, y = Median, color = age_class, fill = age_class)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(
      values = c("J" = "#F4A261", "SA" = "#2A9D8F", "A" = "#264653"),
      labels = c("J" = "Juvenile", "SA" = "Sub-adult", "A" = "Adult")
    ) +
    scale_fill_manual(
      values = c("J" = "#F4A261", "SA" = "#2A9D8F", "A" = "#264653"),
      labels = c("J" = "Juvenile", "SA" = "Sub-adult", "A" = "Adult")
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = paste("Predicted Age Composition vs. DWD Cover —", toupper(dataset)),
      x     = "DWD Cover (%)",
      y     = "Proportion",
      color = "Age Class",
      fill  = "Age Class"
    ) +
    theme_classic() +
    theme(
      axis.text    = element_text(size = 12),
      axis.title   = element_text(size = 13),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 12),
      plot.title   = element_text(face = "bold", size = 13)
    )

  p_age_pred

  ggsave(
    filename = file.path(fig_dir, paste0("pred-age-dwd-", dataset, ".png")),
    plot     = p_age_pred,
    width    = 6, height = 4, dpi      = 300
  )



