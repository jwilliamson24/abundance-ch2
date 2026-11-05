## =================================================
##
## Title: 7-ch2-nmix-posterior-dist.R
##
## Author: Jasmine Williamson
## Date Created: 11/5/2026
##
## Description: Posterior distribution density plots for abundance submodel
##   covariate effects, both species side by side. Mirrors style of
##   Chapter 1 occupancy model posterior distribution plots.
##
## =================================================


## settings ---------------------------------------------------

  rm(list = ls())
  setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/abundance-ch2")

  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)

  out_dir <- "data"
  fig_dir <- "figures"


## load samples -----------------------------------------------

  samples_oss  <- do.call(rbind, lapply(
    readRDS(file.path(out_dir, "mcmc_list_notrt_oss.rds")), as.matrix))
  samples_enes <- do.call(rbind, lapply(
    readRDS(file.path(out_dir, "mcmc_list_notrt_enes.rds")), as.matrix))


## colors -----------------------------------------------------

  param_colors <- c(
    "Canopy Cover" = "#69995D",
    "DWD Count"    = "#62B6CB",
    "Decay Class"  = "#EA7317",
    "Char Class"   = "#E63946",
    "DWD Volume"   = "steelblue"
  )


#### OSS -----------------------------------------------------------------------

  posterior_data_O <- as.data.frame(samples_oss) %>%
    select(beta.canopy, beta.dwd.count, beta.decay, beta.char, beta.vol) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      parameter_clean = case_when(
        parameter == "beta.canopy"    ~ "Canopy Cover",
        parameter == "beta.dwd.count" ~ "DWD Count",
        parameter == "beta.decay"     ~ "Decay Class",
        parameter == "beta.char"      ~ "Char Class",
        parameter == "beta.vol"       ~ "DWD Volume"
      )
    ) %>%
    group_by(parameter_clean) %>%
    mutate(mean_val = mean(value)) %>%
    ungroup() %>%
    mutate(parameter_clean = reorder(parameter_clean, mean_val))

  plot_densities_O <- posterior_data_O %>%
    group_by(parameter_clean) %>%
    do(data.frame(x = density(.$value, n = 512)$x,
                  y = density(.$value, n = 512)$y)) %>%
    mutate(
      y_pos   = as.numeric(parameter_clean),
      y_final = y_pos + (y / max(y)) * 0.8
    ) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_hline(yintercept = 1:length(levels(posterior_data_O$parameter_clean)),
               color = "gray20", alpha = 0.7, linewidth = 0.5) +
    geom_ribbon(aes(x = x, ymin = y_pos, ymax = y_final,
                    group = parameter_clean, fill = parameter_clean),
                alpha = 0.7) +
    scale_fill_manual(values = param_colors, guide = "none") +
    scale_y_continuous(
      breaks = 1:length(levels(posterior_data_O$parameter_clean)),
      labels = levels(posterior_data_O$parameter_clean)
    ) +
    labs(title = "OSS", x = "Effect Size") +
    theme_minimal() +
    theme(
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(vjust = 0, size = 14),
      axis.text.x        = element_text(size = 12),
      axis.title.x       = element_text(size = 14),
      plot.title         = element_text(face = "bold", size = 15)
    )

  print(plot_densities_O)


#### ENES -----------------------------------------------------------------------

  posterior_data_E <- as.data.frame(samples_enes) %>%
    select(beta.canopy, beta.dwd.count, beta.decay, beta.char, beta.vol) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      parameter_clean = case_when(
        parameter == "beta.canopy"    ~ "Canopy Cover",
        parameter == "beta.dwd.count" ~ "DWD Count",
        parameter == "beta.decay"     ~ "Decay Class",
        parameter == "beta.char"      ~ "Char Class",
        parameter == "beta.vol"       ~ "DWD Volume"
      )
    ) %>%
    group_by(parameter_clean) %>%
    mutate(mean_val = mean(value)) %>%
    ungroup() %>%
    mutate(parameter_clean = reorder(parameter_clean, mean_val))

  plot_densities_E <- posterior_data_E %>%
    group_by(parameter_clean) %>%
    do(data.frame(x = density(.$value, n = 512)$x,
                  y = density(.$value, n = 512)$y)) %>%
    mutate(
      y_pos   = as.numeric(parameter_clean),
      y_final = y_pos + (y / max(y)) * 0.8
    ) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_hline(yintercept = 1:length(levels(posterior_data_E$parameter_clean)),
               color = "gray20", alpha = 0.7, linewidth = 0.5) +
    geom_ribbon(aes(x = x, ymin = y_pos, ymax = y_final,
                    group = parameter_clean, fill = parameter_clean),
                alpha = 0.7) +
    scale_fill_manual(values = param_colors, guide = "none") +
    scale_y_continuous(
      breaks = 1:length(levels(posterior_data_E$parameter_clean)),
      labels = levels(posterior_data_E$parameter_clean)
    ) +
    labs(title = "ENES", x = "Effect Size") +
    theme_minimal() +
    theme(
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(vjust = 0, size = 14),
      axis.text.x        = element_text(size = 12),
      axis.title.x       = element_text(size = 14),
      plot.title         = element_text(face = "bold", size = 15)
    )

  print(plot_densities_E)


## combine ----------------------------------------------------

  p <- plot_densities_O + plot_densities_E +
    plot_annotation(
      tag_levels = "A",
      title      = "Posterior Distributions of Covariate Effects on Abundance"
    )

  p
  
  ggsave(
    filename = file.path(fig_dir, "both-posterior-dist-abundance.png"),
    plot     = p,
    width    = 12,
    height   = 8,
    dpi      = 300,
    bg       = "white"
  )

  cat("Saved: both-posterior-dist-abundance.png\n")
