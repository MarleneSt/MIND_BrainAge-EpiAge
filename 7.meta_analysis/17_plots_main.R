#main plots for BrainAge-EpiAge project
# last updated: 07/07/2026 - MS 


#### library ----
library(tidyverse)
library(stringr)
library(metafor)
library(ggplot2)
library(pheatmap)
library(forcats)
library(tidyr)
library(purrr)
library(Polychrome)
library(scales)
library(patchwork)

#### define palettes ----
# [REPRODUCIBILITY] epi model palette (copied from section 2C) defined here, because the
# first (unused) epi violin plot uses it before it is defined further down
model_palette <- c(
  "PCBrainAge" = "#222222",
  "PCGrimAge" = "#f3c300",
  "CorticalClock" = "#875692",
  "cAge" = "#8db600",
  "AltumAge" = "#a1caf1",
  "Horvath2013" = "#be0032",
  "DamAge" = "#c2b280",
  "AdaptAge" = "#848482",
  "Hannum" = "#008856",
  "PhenoAge" = "#e68fac",
  "skinHorvath" = "#0067a5",
  "PedBE" = "#f99379",
  "Wu" = "#604e97",
  "ZhangBLUP" = "#f6a600",
  "ZhangEN" = "#b3446c",
  "Bohlin" = "#dcd300",
  "Knight" = "#654522",
  "EPIC" = "#2b3d26",
  "DunedinPACE" = "#e25822"
)

#### define important folders ----
home <- "data" # [REPRODUCIBILITY] relative path
data_dir <- file.path(home)
results_dir <- "results" # [REPRODUCIBILITY] relative path
plots_dir <- "results/plots" # [REPRODUCIBILITY] relative path
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE) # [REPRODUCIBILITY] create folder on a fresh run


#### load and prep data ----
perf_brain <- readRDS(file.path(results_dir, "perf_brain.rds"))
perf_epi <- readRDS(file.path(results_dir, "perf_epi.rds"))
mod2_eff2 <- readRDS(file.path(results_dir, "mod2_eff2.rds"))

perf_epi<- perf_epi %>%
  mutate(
    model = recode(
      model,
      "Bohlin_conv" = "Bohlin",
      "Knight_conv" = "Knight",
      "EPIC_conv" = "EPIC"
    )
  )


perf_brain_sub <- perf_brain %>% filter(!model %in% c("DunedinPACNI")) %>% droplevels()

perf_epi_sub   <- perf_epi   %>% filter(!model %in% c("DunedinPACE", "DNAmTL")) %>% droplevels()
perf_epi_sub <- perf_epi_sub |>
  mutate(
    log_wMAE = log(wMAE_test),
    # Delta method: Var(log(X)) ≈ (SE_X / X)^2
    log_wMAE_var = (wMAE_SE_boot / wMAE_test)^2
  )

# [REPRODUCIBILITY] log unweighted MAE, needed for the unweighted MAE plots below
# (copied from 04_meta_unweightedMAE.R, lines 84-89; previously only in memory from that script)
perf_epi_sub <- perf_epi_sub |>
  mutate(
    log_MAE = log(MAE),
    # Delta method: Var(log(X)) ≈ (SE_X / X)^2
    log_MAE_var = (MAE_SE_boot / MAE)^2
  )


# rename gestational clocks to not have _conv as suffix 
perf_epi_sub <- perf_epi_sub %>%
  mutate(
    model = recode(
      model,
      "Bohlin_conv" = "Bohlin",
      "Knight_conv" = "Knight",
      "EPIC_conv" = "EPIC"
    )
  )

perf_epi_sub_nogest   <- perf_epi_sub   %>% filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>% droplevels()

#### violin boxplot age ----

##### violin boxplot epi -----
# includes meta estimates 

perf_epi_sub <- perf_epi_sub %>%
  mutate(
    age_bin = case_when(
      mean_age < 2  ~ "Infancy\n<2y",
      mean_age < 6  ~ "Early childhood\n2–5y",
      mean_age < 12 ~ "Middle childhood\n6–11y",
      mean_age < 18 ~ "Adolescence\n12–17y",
      TRUE ~ "Adult\n18+y"
    ),
    age_bin = factor(
      age_bin,
      levels = c(
        "Infancy\n<2y",
        "Early childhood\n2–5y",
        "Middle childhood\n6–11y",
        "Adolescence\n12–17y",
        "Adult\n18+y"
      )
    )
  )


epi_age_bin_counts <- perf_epi_sub %>%
  group_by(age_bin) %>%
  summarise(
    n_cohorts = n_distinct(cohort),
    n_timepoints = n_distinct(timepoint),
    n_cohort_timepoints = n_distinct(paste(cohort, timepoint)),
    n_rows = n(),
    .groups = "drop"
  )

# # A tibble: 5 × 5
# age_bin                   n_cohorts n_timepoints n_cohort_timepoints n_rows
# <fct>                         <int>        <int>               <int>  <int>
#   1 "Infancy\n<2y"                    5           10                  10    168
# 2 "Early childhood\n2–5y"           2            2                   2     29
# 3 "Middle childhood\n6–11y"        10           18                  18    260
# 4 "Adolescence\n12–17y"             7           14                  14    196
# 5 "Adult\n18+y"                     1            1                   1     15

meta_results <- perf_epi_sub %>%
  filter(
    !is.na(log_wMAE),
    !is.na(log_wMAE_var),
    log_wMAE_var > 0
  ) %>%
  group_by(age_bin, model) %>%
  group_modify(~{
    
    dat <- .x
    
    if (nrow(dat) == 1) {
      return(tibble(
        pooled_log_mae = dat$log_wMAE[1],
        ci_lb_log = NA_real_,
        ci_ub_log = NA_real_,
        pooled_mae = exp(dat$log_wMAE[1]),
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        k = 1,
        meta_model = "single estimate"
      ))
    }
    
    res <- rma.mv(
      yi = log_wMAE,
      V = log_wMAE_var,
      random = ~ 1 | cohort/timepoint,
      data = dat,
      method = "REML"
    )
    
    tibble(
      pooled_log_mae = res$b[1, 1],
      ci_lb_log = res$ci.lb,
      ci_ub_log = res$ci.ub,
      pooled_mae = exp(res$b[1, 1]),
      ci_lb = exp(res$ci.lb),
      ci_ub = exp(res$ci.ub),
      k = res$k,
      meta_model = "random-effects"
    )
  }) %>%
  ungroup()

pd <- position_dodge(width = 0.7)

png(
  filename = 
    file.path(plots_dir, "plot_violin_perf_epi_age_brackets.png"),
  width = 2340,
  height = 1950,
  res = 300
)
ggplot(perf_epi_sub, aes(x = age_bin, y = log_wMAE)) +
  geom_violin(fill = "grey88", colour = NA, trim = FALSE) +
  
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = NA
  ) +
  
  geom_errorbar(
    data = meta_results,
    aes(
      x = age_bin,
      ymin = ci_lb_log,
      ymax = ci_ub_log,
      colour = model,
      group = model
    ),
    position = pd,
    width = 0.08,
    alpha = 0.5,
    na.rm = TRUE,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = meta_results,
    aes(
      x = age_bin,
      y = pooled_log_mae,
      colour = model,
      group = model
    ),
    position = pd,
    alpha = 0.7,
    size = 2.4,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(values = model_palette) +
  
  labs(
    x = "Developmental age group",
    y = "log(weighted MAE)",
    title = "Epigenetic clock performance across age groups",
    colour = "Model"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank()
  )

dev.off()

##### USED violin boxplot epi more brackets (2C) -----
# includes meta estimates 

model_palette <- c(
  "PCBrainAge" = "#222222",
  "PCGrimAge" = "#f3c300",
  "CorticalClock" = "#875692",
  "cAge" = "#8db600",
  "AltumAge" = "#a1caf1",
  "Horvath2013" = "#be0032",
  "DamAge" = "#c2b280",
  "AdaptAge" = "#848482",
  "Hannum" = "#008856",
  "PhenoAge" = "#e68fac",
  "skinHorvath" = "#0067a5",
  "PedBE" = "#f99379",
  "Wu" = "#604e97",
  "ZhangBLUP" = "#f6a600",
  "ZhangEN" = "#b3446c",
  "Bohlin" = "#dcd300",
  "Knight" = "#654522",
  "EPIC" = "#2b3d26",
  "DunedinPACE" = "#e25822"
)

perf_epi_sub <- perf_epi_sub %>%
  mutate(
    age_bin = case_when(
      mean_age < 0.25  ~ "birth\n<3m",
      mean_age < 3  ~ "Infancy\n<3y",
      mean_age < 6  ~ "Early childhood\n3–5y",
      mean_age < 10 ~ "Middle childhood\n6–9y",
      mean_age < 15 ~ "Late childhood\n10–14y",
      mean_age < 18 ~ "Adolescence\n15–17y",
      TRUE ~ "Young adulthood\n18+y"
    ),
    age_bin = factor(
      age_bin,
      levels = c(
        "birth\n<3m",
        "Infancy\n<3y",
        "Early childhood\n3–5y",
        "Middle childhood\n6–9y",
        "Late childhood\n10–14y",
        "Adolescence\n15–17y",
        "Young adulthood\n18+y"
      )
    )
  )

epi_age_bin_counts <- perf_epi_sub %>%
  group_by(age_bin) %>%
  summarise(
    k = n_distinct(paste(cohort, timepoint)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("k=", k))

meta_results <- perf_epi_sub %>%
  filter(
    !is.na(log_wMAE),
    !is.na(log_wMAE_var),
    log_wMAE_var > 0
  ) %>%
  group_by(age_bin, model) %>%
  group_modify(~{
    
    dat <- .x
    
    if (nrow(dat) == 1) {
      return(tibble(
        pooled_log_mae = dat$log_wMAE[1],
        ci_lb_log = NA_real_,
        ci_ub_log = NA_real_,
        pooled_mae = exp(dat$log_wMAE[1]),
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        k = 1,
        meta_model = "single estimate"
      ))
    }
    
    res <- rma.mv(
      yi = log_wMAE,
      V = log_wMAE_var,
      random = ~ 1 | cohort/timepoint,
      data = dat,
      method = "REML"
    )
    
    tibble(
      pooled_log_mae = res$b[1, 1],
      ci_lb_log = res$ci.lb,
      ci_ub_log = res$ci.ub,
      pooled_mae = exp(res$b[1, 1]),
      ci_lb = exp(res$ci.lb),
      ci_ub = exp(res$ci.ub),
      k = res$k,
      meta_model = "random-effects"
    )
  }) %>%
  ungroup()

pd <- position_dodge(width = 0.7)

# plot went into axis labels originally so had to adapt coordinates based on CIs
y_min <- min(
  perf_epi_sub$log_wMAE,
  meta_results$ci_lb_log,
  meta_results$pooled_log_mae,
  na.rm = TRUE
)

y_max <- max(
  perf_epi_sub$log_wMAE,
  meta_results$ci_ub_log,
  meta_results$pooled_log_mae,
  na.rm = TRUE
)

y_pad <- 0.08 * (y_max - y_min)


png(
  filename = 
    file.path(plots_dir, "plot_violin_perf_epi_age_brackets_morebrackets_k.png"),
  width = 2340,
  height = 1950,
  res = 300
)
ggplot(perf_epi_sub, aes(x = age_bin, y = log_wMAE)) +
  
  geom_hline(
    yintercept = 0,
    colour = "grey92",
    linewidth = 0.5
  ) +
  
  geom_violin(fill = "grey85", colour = NA, trim = FALSE) +

  
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = NA
  ) +
  
  geom_errorbar(
    data = meta_results,
    aes(
      x = age_bin,
      ymin = ci_lb_log,
      ymax = ci_ub_log,
      colour = model,
      group = model
    ),
    position = pd,
    width = 0.08,
    alpha = 0.45,
    na.rm = TRUE,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = meta_results,
    aes(
      x = age_bin,
      y = pooled_log_mae,
      colour = model,
      group = model
    ),
    position = pd,
    alpha = 0.6,
    size = 2.4,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(values = model_palette) +
  
  geom_text(
    data = epi_age_bin_counts,
    aes(
      x = age_bin,
      y = y_min - 0.35 * y_pad,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = 1,
    size = 3.5
  ) +
  
  coord_cartesian(
    ylim = c(
      y_min - y_pad,
      y_max + y_pad
    ),
    clip = "off"
  ) +

  
  labs(
    x = "Developmental age group",
    y = "log(weighted MAE)",
    title = "Epigenetic clock performance across age groups",
    colour = "Model"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  )
dev.off()

##### USED violin boxplot epi more brackets no GEST -----
# includes meta estimates 

perf_epi_sub_nogest   <- perf_epi_sub   %>% filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>% droplevels()


perf_epi_sub_nogest <- perf_epi_sub_nogest %>%
  mutate(
    age_bin = case_when(
      mean_age < 0.25  ~ "birth\n<3m",
      mean_age < 3  ~ "Infancy\n<3y",
      mean_age < 6  ~ "Early childhood\n3–5y",
      mean_age < 10 ~ "Middle childhood\n6–9y",
      mean_age < 15 ~ "Late childhood\n10–14y",
      mean_age < 18 ~ "Adolescence\n15–17y",
      TRUE ~ "Young adulthood\n18+y"
    ),
    age_bin = factor(
      age_bin,
      levels = c(
        "birth\n<3m",
        "Infancy\n<3y",
        "Early childhood\n3–5y",
        "Middle childhood\n6–9y",
        "Late childhood\n10–14y",
        "Adolescence\n15–17y",
        "Young adulthood\n18+y"
      )
    )
  )

epi_age_bin_counts <- perf_epi_sub_nogest %>%
  group_by(age_bin) %>%
  summarise(
    k = n_distinct(paste(cohort, timepoint)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("k=", k))

meta_results <- perf_epi_sub_nogest %>%
  filter(
    !is.na(log_wMAE),
    !is.na(log_wMAE_var),
    log_wMAE_var > 0
  ) %>%
  group_by(age_bin, model) %>%
  group_modify(~{
    
    dat <- .x
    
    if (nrow(dat) == 1) {
      return(tibble(
        pooled_log_mae = dat$log_wMAE[1],
        ci_lb_log = NA_real_,
        ci_ub_log = NA_real_,
        pooled_mae = exp(dat$log_wMAE[1]),
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        k = 1,
        meta_model = "single estimate"
      ))
    }
    
    res <- rma.mv(
      yi = log_wMAE,
      V = log_wMAE_var,
      random = ~ 1 | cohort/timepoint,
      data = dat,
      method = "REML"
    )
    
    tibble(
      pooled_log_mae = res$b[1, 1],
      ci_lb_log = res$ci.lb,
      ci_ub_log = res$ci.ub,
      pooled_mae = exp(res$b[1, 1]),
      ci_lb = exp(res$ci.lb),
      ci_ub = exp(res$ci.ub),
      k = res$k,
      meta_model = "random-effects"
    )
  }) %>%
  ungroup()

pd <- position_dodge(width = 0.7)

y_min <- min(
  perf_epi_sub_nogest$log_wMAE,
  meta_results$ci_lb_log,
  meta_results$pooled_log_mae,
  na.rm = TRUE
)

y_max <- max(
  perf_epi_sub_nogest$log_wMAE,
  meta_results$ci_ub_log,
  meta_results$pooled_log_mae,
  na.rm = TRUE
)

y_pad <- 0.08 * (y_max - y_min)

png(
  filename = 
    file.path(plots_dir, "plot_violin_perf_epi_age_brackets_morebrackets_k_nogest.png"),
  width = 2340,
  height = 1950,
  res = 300
)
ggplot(perf_epi_sub_nogest, aes(x = age_bin, y = log_wMAE)) +
  
  geom_hline(
    yintercept = 0,
    colour = "grey92",
    linewidth = 0.5
  ) +
  
  geom_violin(fill = "grey85", colour = NA, trim = FALSE) +
  
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = NA
  ) +
  
  geom_errorbar(
    data = meta_results,
    aes(
      x = age_bin,
      ymin = ci_lb_log,
      ymax = ci_ub_log,
      colour = model,
      group = model
    ),
    position = pd,
    width = 0.08,
    alpha = 0.45,
    na.rm = TRUE,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = meta_results,
    aes(
      x = age_bin,
      y = pooled_log_mae,
      colour = model,
      group = model
    ),
    position = pd,
    alpha = 0.6,
    size = 2.4,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(values = model_palette) +
  
  geom_text(
    data = epi_age_bin_counts,
    aes(
      x = age_bin,
      y = y_min - 0.35 * y_pad,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = 1,
    size = 3.5
  ) +
  
  coord_cartesian(
    ylim = c(
      y_min - y_pad,
      y_max + y_pad
    ),
    clip = "off"
  ) +
  
  labs(
    x = "Developmental age group",
    y = "log(weighted MAE)",
    title = "Epigenetic clock performance across age groups",
    colour = "Model"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  )
dev.off()

##### USED violin boxplot epi unweighted MAE more brackets -----
# includes meta estimates 

perf_epi_sub <- perf_epi_sub %>%
  mutate(
    age_bin = case_when(
      mean_age < 0.25  ~ "birth\n<3m",
      mean_age < 3  ~ "Infancy\n<3y",
      mean_age < 6  ~ "Early childhood\n3–5y",
      mean_age < 10 ~ "Middle childhood\n6–9y",
      mean_age < 15 ~ "Late childhood\n10–14y",
      mean_age < 18 ~ "Adolescence\n15–17y",
      TRUE ~ "Young adulthood\n18+y"
    ),
    age_bin = factor(
      age_bin,
      levels = c(
        "birth\n<3m",
        "Infancy\n<3y",
        "Early childhood\n3–5y",
        "Middle childhood\n6–9y",
        "Late childhood\n10–14y",
        "Adolescence\n15–17y",
        "Young adulthood\n18+y"
      )
    )
  )

epi_age_bin_counts <- perf_epi_sub %>%
  group_by(age_bin) %>%
  summarise(
    k = n_distinct(paste(cohort, timepoint)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("k=", k))

meta_results <- perf_epi_sub %>%
  filter(
    !is.na(log_MAE),
    !is.na(log_MAE_var),
    log_MAE_var > 0
  ) %>%
  group_by(age_bin, model) %>%
  group_modify(~{
    
    dat <- .x
    
    if (nrow(dat) == 1) {
      return(tibble(
        pooled_log_mae = dat$log_MAE[1],
        ci_lb_log = NA_real_,
        ci_ub_log = NA_real_,
        pooled_mae = exp(dat$log_MAE[1]),
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        k = 1,
        meta_model = "single estimate"
      ))
    }
    
    res <- rma.mv(
      yi = log_MAE,
      V = log_MAE_var,
      random = list(~ 1 | cohort/timepoint),
      data = dat,
      method = "REML",
      control = list(
        optimizer = "optim",
        optmethod = "BFGS",
        maxit = 10000
      )
    )
    
    tibble(
      pooled_log_mae = res$b[1, 1],
      ci_lb_log = res$ci.lb,
      ci_ub_log = res$ci.ub,
      pooled_mae = exp(res$b[1, 1]),
      ci_lb = exp(res$ci.lb),
      ci_ub = exp(res$ci.ub),
      k = res$k,
      meta_model = "random-effects"
    )
  }) %>%
  ungroup()

pd <- position_dodge(width = 0.7)

png(
  filename = 
    file.path(plots_dir, "plot_violin_perf_epi_age_unweightedMAE_brackets_morebrackets_k.png"),
  width = 2340,
  height = 1950,
  res = 300
)
ggplot(perf_epi_sub, aes(x = age_bin, y = log_MAE)) +
  geom_violin(fill = "grey85", colour = NA, trim = FALSE) +
  
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = NA
  ) +
  
  geom_errorbar(
    data = meta_results,
    aes(
      x = age_bin,
      ymin = ci_lb_log,
      ymax = ci_ub_log,
      colour = model,
      group = model
    ),
    position = pd,
    width = 0.08,
    alpha = 0.45,
    na.rm = TRUE,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = meta_results,
    aes(
      x = age_bin,
      y = pooled_log_mae,
      colour = model,
      group = model
    ),
    position = pd,
    alpha = 0.6,
    size = 2.4,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(values = model_palette) +
  
  geom_text(
    data = epi_age_bin_counts,
    aes(
      x = age_bin,
      y = min(perf_epi_sub$log_MAE, na.rm = TRUE),
      label = label
    ),
    inherit.aes = FALSE,
    vjust = 1.8,
    size = 3.5
  ) +
  
  coord_cartesian(
    ylim = c(
      min(perf_epi_sub$log_MAE, na.rm = TRUE) - 0.35,
      max(perf_epi_sub$log_MAE, na.rm = TRUE)
    ),
    clip = "off"
  ) +
  
  labs(
    x = "Developmental age group",
    y = "log(MAE)",
    title = "Epigenetic clock performance across age groups",
    colour = "Model"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  )
dev.off()

##### USED violin boxplot brain unweighted MAE more brackets -----
# includes meta estimates

model_palette <- c(
  "Pyment" = "#222222",
  "DevBrainAge" = "#f3c300",
  "Centile2" = "#875692",
  "DBN" = "#8db600",
  "Kaufmann" = "#a1caf1",
  "PyBrainAge" = "#be0032",
  "ENIGMA" = "#c2b280"
)

perf_brain_sub <- perf_brain_sub %>%
  mutate(
    age_bin = case_when(
      mean_age < 6  ~ "Early childhood\n4–5y",
      mean_age < 10 ~ "Middle childhood\n6–9y",
      mean_age < 15 ~ "Late childhood\n10–14y",
      mean_age < 18 ~ "Adolescence\n15–17y",
      TRUE ~ "Young adulthood\n18+y"
    ),
    age_bin = factor(
      age_bin,
      levels = c(
        "Early childhood\n4–5y",
        "Middle childhood\n6–9y",
        "Late childhood\n10–14y",
        "Adolescence\n15–17y",
        "Young adulthood\n18+y"
      )
    )
  )

brain_age_bin_counts <- perf_brain_sub %>%
  group_by(age_bin) %>%
  summarise(
    k = n_distinct(paste(cohort, timepoint)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("k=", k))

meta_results <- perf_brain_sub %>%
  filter(
    !is.na(MAE),
    !is.na(MAE_SE_boot),
    MAE_SE_boot > 0
  ) %>%
  group_by(age_bin, model) %>%
  group_modify(~{
    
    dat <- .x
    
    if (nrow(dat) == 1) {
      return(tibble(
        pooled_mae = dat$MAE[1],
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        k = 1,
        meta_model = "single estimate"
      ))
    }
    
    res <- rma.mv(
      yi = MAE,
      V = MAE_SE_boot^2,
      random = list(
        ~ 1 | cohort/timepoint
      ),
      data = dat,
      method = "REML"
    )
    
    tibble(
      pooled_mae = res$b[1, 1],
      ci_lb = res$ci.lb,
      ci_ub = res$ci.ub,
      k = res$k,
      meta_model = "random-effects"
    )
  }) %>%
  ungroup()

pd <- position_dodge(width = 0.7)

# plot went into axis labels originally so had to adapt coordinates based on CIs
y_min <- min(
  perf_brain_sub$MAE,
  meta_results$ci_lb,
  meta_results$pooled_mae,
  na.rm = TRUE
)

y_max <- max(
  perf_brain_sub$MAE,
  meta_results$ci_ub,
  meta_results$pooled_mae,
  na.rm = TRUE
)

y_pad <- 0.08 * (y_max - y_min)

png(
  filename = file.path(
    plots_dir,
    "plot_violin_perf_brain_age_unweightedMAE_brackets_morebrackets_k.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot(perf_brain_sub, aes(x = age_bin, y = MAE)) +
  
  geom_hline(
    yintercept = 0,
    colour = "grey92",
    linewidth = 0.5
  ) +
  
  geom_violin(fill = "grey85", colour = NA, trim = FALSE) +
  
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = NA
  ) +
  
  geom_errorbar(
    data = meta_results,
    aes(
      x = age_bin,
      ymin = ci_lb,
      ymax = ci_ub,
      colour = model,
      group = model
    ),
    position = pd,
    width = 0.08,
    alpha = 0.45,
    na.rm = TRUE,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = meta_results,
    aes(
      x = age_bin,
      y = pooled_mae,
      colour = model,
      group = model
    ),
    position = pd,
    alpha = 0.6,
    size = 2.4,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(values = model_palette) +
  
  geom_text(
    data = brain_age_bin_counts,
    aes(
      x = age_bin,
      y = y_min - 0.35 * y_pad,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = 1,
    size = 3.5
  ) +
  
  coord_cartesian(
    ylim = c(
      y_min - y_pad,
      y_max + y_pad
    ),
    clip = "off"
  ) +

  
  labs(
    x = "Developmental age group",
    y = "MAE",
    title = "Brain age model performance across age groups",
    colour = "Model"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  )

dev.off()

##### USED violin boxplot brain wMAE more brackets (3C) -----
# includes meta estimates

model_palette <- c(
  "Pyment" = "#222222",
  "DevBrainAge" = "#f3c300",
  "Centile2" = "#875692",
  "DBN" = "#8db600",
  "Kaufmann" = "#a1caf1",
  "PyBrainAge" = "#be0032",
  "ENIGMA" = "#c2b280"
)

perf_brain_sub <- perf_brain_sub %>%
  mutate(
    age_bin = case_when(
      mean_age < 6  ~ "Early childhood\n4–5y",
      mean_age < 10 ~ "Middle childhood\n6–9y",
      mean_age < 15 ~ "Late childhood\n10–14y",
      mean_age < 18 ~ "Adolescence\n15–17y",
      TRUE ~ "Young adulthood\n18+y"
    ),
    age_bin = factor(
      age_bin,
      levels = c(
        "Early childhood\n4–5y",
        "Middle childhood\n6–9y",
        "Late childhood\n10–14y",
        "Adolescence\n15–17y",
        "Young adulthood\n18+y"
      )
    )
  )

brain_age_bin_counts <- perf_brain_sub %>%
  group_by(age_bin) %>%
  summarise(
    k = n_distinct(paste(cohort, timepoint)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("k=", k))

meta_results <- perf_brain_sub %>%
  filter(
    !is.na(wMAE_test),
    !is.na(wMAE_SE_boot),
    wMAE_SE_boot > 0
  ) %>%
  group_by(age_bin, model) %>%
  group_modify(~{
    
    dat <- .x
    
    if (nrow(dat) == 1) {
      return(tibble(
        pooled_mae = dat$wMAE_test[1],
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        k = 1,
        meta_model = "single estimate"
      ))
    }
    
    res <- rma.mv(
      yi = wMAE_test,
      V = wMAE_SE_boot^2,
      random = list(
        ~ 1 | cohort/timepoint
      ),
      data = dat,
      method = "REML"
    )
    
    tibble(
      pooled_mae = res$b[1, 1],
      ci_lb = res$ci.lb,
      ci_ub = res$ci.ub,
      k = res$k,
      meta_model = "random-effects"
    )
  }) %>%
  ungroup()

pd <- position_dodge(width = 0.7)

# plot went into axis labels originally so had to adapt coordinates based on CIs
y_min <- min(
  perf_brain_sub$wMAE_test,
  meta_results$ci_lb,
  meta_results$pooled_mae,
  na.rm = TRUE
)

y_max <- max(
  perf_brain_sub$wMAE_test,
  meta_results$ci_ub,
  meta_results$pooled_mae,
  na.rm = TRUE
)

y_pad <- 0.08 * (y_max - y_min)

png(
  filename = file.path(
    plots_dir,
    "plot_violin_perf_brain_age_wMAE_brackets_morebrackets_k.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot(perf_brain_sub, aes(x = age_bin, y = wMAE_test)) +
  
  geom_hline(
    yintercept = 0,
    colour = "grey92",
    linewidth = 0.5
  ) +
  
  geom_violin(fill = "grey85", colour = NA, trim = FALSE) +
  
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = NA
  ) +
  
  geom_errorbar(
    data = meta_results,
    aes(
      x = age_bin,
      ymin = ci_lb,
      ymax = ci_ub,
      colour = model,
      group = model
    ),
    position = pd,
    width = 0.08,
    alpha = 0.45,
    na.rm = TRUE,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = meta_results,
    aes(
      x = age_bin,
      y = pooled_mae,
      colour = model,
      group = model
    ),
    position = pd,
    alpha = 0.6,
    size = 2.4,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(values = model_palette) +
  
  geom_text(
    data = brain_age_bin_counts,
    aes(
      x = age_bin,
      y = y_min - 0.35 * y_pad,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = 1,
    size = 3.5
  ) +
  
  coord_cartesian(
    ylim = c(
      y_min - y_pad,
      y_max + y_pad
    ),
    clip = "off"
  ) +
  
  
  labs(
    x = "Developmental age group",
    y = "weighted MAE",
    title = "Brain age model performance across age groups",
    colour = "Model"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  )

dev.off()

#### performance with age info ----

##### main brain results -----

# load meta results 
brain_wMAE_mv_global <- readRDS(file.path(results_dir, "perf_brain_wMAE_noDunedin_global_final/model_brain_wMAE_mv_global.rds"))
brain_wMAE_mv_modelwise <- readRDS(file.path(results_dir, "perf_brain_wMAE_noDunedin_modelwise_final/model_brain_wMAE_mv_modelwise.rds"))

perf_brain_sub_plotting <- perf_brain_sub %>%
  mutate(
    timepoint_ord = mean_age
  )

cohort_levels <- unique(perf_brain_sub_plotting$cohort)

cohort_shapes <- setNames(
  rep(c(0,1,2,3,4,5,6,15, 16, 17, 18, 7, 8, 9,10, 14), length.out = length(cohort_levels)),
  cohort_levels
)

# extracting meta mean 
overall_est <- data.frame(
  model = "Pooled",
  wMAE  = as.numeric(coef(brain_wMAE_mv_global)),
  ci.lb = as.numeric(brain_wMAE_mv_global$ci.lb),
  ci.ub = as.numeric(brain_wMAE_mv_global$ci.ub)
)

# extracting modelwise 
brain_wMAE_est <- tibble(
  term  = rownames(brain_wMAE_mv_modelwise$b),
  wMAE     = as.numeric(brain_wMAE_mv_modelwise$b),
  se    = brain_wMAE_mv_modelwise$se,
  ci.lb = brain_wMAE_mv_modelwise$ci.lb,
  ci.ub = brain_wMAE_mv_modelwise$ci.ub,
  pval  = brain_wMAE_mv_modelwise$pval
) %>%
  mutate(model = str_remove(term, "^model")) %>%
  select(-term) %>%
  left_join(perf_brain_sub %>% count(model, name = "k"), by = "model") %>%
  arrange(wMAE)

brain_wMAE_est <- brain_wMAE_est %>%
  mutate(label = paste0(model, " (k=", k, ")"))

# --- combine modelwise + overall ---

plot_est <- bind_rows(brain_wMAE_est, overall_est)

# --- ordering ---
# models ordered by generation and wMAE,
# overall forced to bottom
# sort models from best to worst (lowest wMAE = best)

model_order <- brain_wMAE_est %>%
  arrange(wMAE) %>%
  pull(model)

# we want worst just above Overall
# so reverse model order
model_order <- rev(model_order)

# Overall must be FIRST level to appear at bottom
model_order <- c("Pooled", model_order)

# apply factor levels
perf_brain_sub_plotting$model <- factor(
  as.character(perf_brain_sub_plotting$model),
  levels = model_order
)

plot_est$model <- factor(
  as.character(plot_est$model),
  levels = model_order
)

# divider above overall row
overall_divider_position <- 1.5

cohort_levels <- unique(perf_brain_sub_plotting$cohort)


png(
  filename = file.path(
    plots_dir,
    "plot_perf_brain_modelwise_plus_overall_meta_raw_combined.png"
  ),
  width = 2340,
  height = 2340,
  res = 300
)

ggplot() +
  
  geom_point(
    data = perf_brain_sub_plotting,
    aes(
      x = wMAE_test,
      y = model,
      shape = cohort,
      colour = timepoint_ord
    ),
    alpha = .65,
    size = 2.8,
    position = position_jitter(height = .12, width = 0)
  ) +
  
  geom_segment(
    data = plot_est,
    aes(x = ci.lb, xend = ci.ub, y = model, yend = model),
    colour = "black",
    linewidth = 1.1
  ) +
  
  geom_point(
    data = plot_est,
    aes(x = wMAE, y = model),
    colour = "black",
    fill = "black",
    shape = 23,
    size = 4.2
  ) + 
  
  geom_hline(
    yintercept = overall_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  scale_y_discrete(limits = levels(plot_est$model)) +
  
  scale_colour_viridis_c(
    option = "magma",
    direction = 1,
    name = "Mean age\nat time point"
  ) +
  
  scale_shape_manual(
    values = cohort_shapes,
    name = "Cohort"
  ) +
  
  labs(
    title = "Brain Age Model Performance",
    x = "Weighted MAE (MAE ÷ age range)",
    y = NULL
  ) +
  
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.background = element_rect(
      fill = NA,
      colour = NA
    ),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust=0.5),
    legend.key.height = unit(0.45, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank()
  ) +
  
  guides(
    shape = guide_legend(
      ncol = 1,
      override.aes = list(
        fill = "grey75",
        colour = "grey25",
        alpha = 1,
        size = 3
      )
    ),
    color = guide_colourbar(
      barheight = unit(5, "cm"),
      title.position = "top", 
      reverse = TRUE
    )
  ) +
  scale_x_continuous(
    limits = c(0, 15),
    oob = scales::squish,
    breaks = seq(0, 15, 5),
    labels = function(x) ifelse(x == 15, "15+", x)
  )
dev.off()

##### main brain results alt -----

# load meta results 
brain_wMAE_mv_global <- readRDS(file.path(results_dir, "perf_brain_wMAE_noDunedin_global_final/model_brain_wMAE_mv_global.rds"))
brain_wMAE_mv_modelwise <- readRDS(file.path(results_dir, "perf_brain_wMAE_noDunedin_modelwise_final/model_brain_wMAE_mv_modelwise.rds"))

perf_brain_sub_plotting <- perf_brain_sub %>%
  mutate(
    timepoint_ord = mean_age
  )

cohort_levels <- unique(perf_brain_sub_plotting$cohort)

cohort_shapes <- setNames(
  rep(c(15, 16, 17, 18, 3, 4, 7, 8, 9, 10, 12, 13, 14), length.out = length(cohort_levels)),
  cohort_levels
)

# extracting meta mean 
overall_est <- data.frame(
  model = "Pooled",
  wMAE  = as.numeric(coef(brain_wMAE_mv_global)),
  ci.lb = as.numeric(brain_wMAE_mv_global$ci.lb),
  ci.ub = as.numeric(brain_wMAE_mv_global$ci.ub)
)

# extracting modelwise 
brain_wMAE_est <- tibble(
  term  = rownames(brain_wMAE_mv_modelwise$b),
  wMAE     = as.numeric(brain_wMAE_mv_modelwise$b),
  se    = brain_wMAE_mv_modelwise$se,
  ci.lb = brain_wMAE_mv_modelwise$ci.lb,
  ci.ub = brain_wMAE_mv_modelwise$ci.ub,
  pval  = brain_wMAE_mv_modelwise$pval
) %>%
  mutate(model = str_remove(term, "^model")) %>%
  select(-term) %>%
  left_join(perf_brain_sub %>% count(model, name = "k"), by = "model") %>%
  arrange(wMAE)

brain_wMAE_est <- brain_wMAE_est %>%
  mutate(label = paste0(model, " (k=", k, ")"))

# --- combine modelwise + overall ---

plot_est <- bind_rows(brain_wMAE_est, overall_est)

# --- ordering ---
# models ordered by generation and wMAE,
# overall forced to bottom
# sort models from best to worst (lowest wMAE = best)

model_order <- brain_wMAE_est %>%
  arrange(wMAE) %>%
  pull(model)

# we want worst just above Overall
# so reverse model order
model_order <- rev(model_order)

# Overall must be FIRST level to appear at bottom
model_order <- c("Pooled", model_order)

# apply factor levels
perf_brain_sub_plotting$model <- factor(
  as.character(perf_brain_sub_plotting$model),
  levels = model_order
)

plot_est$model <- factor(
  as.character(plot_est$model),
  levels = model_order
)

# divider above overall row
overall_divider_position <- 1.5

cohort_levels <- unique(perf_brain_sub_plotting$cohort)

png(
  filename = file.path(
    plots_dir,
    "plot_perf_brain_modelwise_plus_overall_meta_raw_combined_alt.png"
  ),
  width = 2340,
  height = 2340,
  res = 300
)

ggplot() +
  
  geom_point(
    data = perf_brain_sub_plotting,
    aes(
      x = wMAE_test,
      y = model,
      shape = cohort,
      color = timepoint_ord
    ),
    alpha = .65,
    size = 2.8,
    position = position_jitter(height = .12, width = 0)
  ) +
  
  geom_segment(
    data = plot_est,
    aes(x = ci.lb, xend = ci.ub, y = model, yend = model),
    colour = "black",
    linewidth = 1.1
  ) +
  
  geom_point(
    data = plot_est,
    aes(x = wMAE, y = model),
    colour = "black",
    fill = "black",
    shape = 23,
    size = 4.2
  ) + 
  
  geom_hline(
    yintercept = overall_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  scale_y_discrete(limits = levels(plot_est$model)) +
  
  scale_color_viridis_c(
    option = "magma",
    direction = 1,
    name = "Timepoint\nyounger = darker"
  ) +
  
  scale_shape_manual(
    values = cohort_shapes,
    name = "Cohort"
  ) +
  
  labs(
    title = "Brain Age Model Performance",
    x = "Weighted MAE (MAE ÷ age range)",
    y = NULL
  ) +
  
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.background = element_rect(
      fill = NA,
      colour = NA
    ),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11),
    legend.key.height = unit(0.45, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank()
  ) +
  
  guides(
    shape = guide_legend(
      ncol = 1,
      override.aes = list(
        fill = "grey75",
        colour = "grey25",
        alpha = 1,
        size = 3
      )
    ),
    color = guide_colourbar(
      barheight = unit(5, "cm"),
      title.position = "top",
      reverse = TRUE
    )
   ) + 
  scale_x_continuous(
    trans = scales::log_trans(base = exp(1)),
    breaks = c(1, 2, 5, 10, 30),
    labels = scales::label_number()
  )
#+ 
  # scale_x_continuous(
  #   limits = c(0, 15),
  #   oob = scales::squish,
  #   breaks = seq(0, 15, 5),
  #   labels = function(x) ifelse(x == 15, "15+", x)
  # )

dev.off()

#### performance with just cohort info ----

##### USED brain perf with just cohort info (3A) -----

# load meta results 
brain_wMAE_mv_global <- readRDS(file.path(results_dir, "perf_brain_wMAE_noDunedin_global_final/model_brain_wMAE_mv_global.rds"))
brain_wMAE_mv_modelwise <- readRDS(file.path(results_dir, "perf_brain_wMAE_noDunedin_modelwise_final/model_brain_wMAE_mv_modelwise.rds"))


perf_brain_sub_plotting <- perf_brain_sub

cohort_levels <- unique(perf_brain_sub_plotting$cohort)

cohort_palette <- c(
  "ALSPAC" = "#222222",
  "BHRC" = "#f3c300",
  "CannTeen" = "#875692",
  "DCHS" = "#f38400",
  "FinnBrain" = "#a1caf1",
  "GenR" = "#be0032",
  "GUSTO" = "#c2b280",
  "K2H Childhood" = "#848482",
  "K2H Infancy" = "#008856",
  "MTwiNS" = "#e68fac",
  "NICAP" = "#0067a5",
  "Oregon ADHD-1000" = "#f99379",
  "FFCWS" = "#604e97",
  "TAG" = "#f6a600",
  "UCI Echo" = "#b3446c"
)

# extracting meta mean 
overall_est <- data.frame(
  model = "Pooled",
  wMAE  = as.numeric(coef(brain_wMAE_mv_global)),
  ci.lb = as.numeric(brain_wMAE_mv_global$ci.lb),
  ci.ub = as.numeric(brain_wMAE_mv_global$ci.ub)
)

# extracting modelwise 
brain_wMAE_est <- tibble(
  term  = rownames(brain_wMAE_mv_modelwise$b),
  wMAE     = as.numeric(brain_wMAE_mv_modelwise$b),
  se    = brain_wMAE_mv_modelwise$se,
  ci.lb = brain_wMAE_mv_modelwise$ci.lb,
  ci.ub = brain_wMAE_mv_modelwise$ci.ub,
  pval  = brain_wMAE_mv_modelwise$pval
) %>%
  mutate(model = str_remove(term, "^model")) %>%
  select(-term) %>%
  left_join(perf_brain_sub %>% count(model, name = "k"), by = "model") %>%
  arrange(wMAE)

brain_wMAE_est <- brain_wMAE_est %>%
  mutate(label = paste0(model, " (k=", k, ")"))

# --- combine modelwise + overall ---
plot_est <- bind_rows(brain_wMAE_est, overall_est)

# ---  ordering ---
# models ordered by generation and wMAE,
# overall forced to bottom
# sort models from best to worst (lowest wMAE = best)

model_order <- brain_wMAE_est %>%
  arrange(wMAE) %>%
  pull(model)

# we want worst just above Overall
# so reverse model order
model_order <- rev(model_order)

# Overall must be FIRST level to appear at bottom
model_order <- c("Pooled", model_order)

# apply factor levels
perf_brain_sub_plotting$model <- factor(
  as.character(perf_brain_sub_plotting$model),
  levels = model_order
)

plot_est$model <- factor(
  as.character(plot_est$model),
  levels = model_order
)

# divider above overall row
overall_divider_position <- 1.5

cohort_levels <- unique(perf_brain_sub_plotting$cohort)

png(
  filename = file.path(
    plots_dir,
    "plot_perf_brain_modelwise_plus_overall_meta_raw_combined_cohortcolors.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot() +
  
  # raw cohort-level results
  geom_point(
    data = perf_brain_sub_plotting,
    aes(
      x = wMAE_test,
      y = model,
      colour = cohort
    ),
    alpha = .45,
    size = 2,
    position = position_jitter(height = .12)
  ) +
  
  # meta CIs
  geom_segment(
    data = plot_est,
    aes(
      x = ci.lb,
      xend = ci.ub,
      y = model,
      yend = model,
      colour = "Overall pooled"
    ),
    linewidth = 0.9
  ) +
  
  # meta pooled estimate
  geom_point(
    data = plot_est,
    aes(
      x = wMAE,
      y = model,
      colour = "Overall pooled"
    ),
    shape = 18,
    size = 5
  ) +
  
  
  
  # divider above Overall
  geom_hline(
    yintercept = overall_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  scale_y_discrete(limits = levels(plot_est$model)) +
  
  scale_colour_manual(
    values = c(
      cohort_palette,
      "Overall pooled" = "black"
    ),
    breaks = c(names(cohort_palette), "Overall pooled")
  ) +
  
  labs(
    title = "Brain age model performance (weighted MAE)",
    x = "Weighted MAE (MAE ÷ age range)",
    y = NULL,
    colour = "Source"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    # legend.background = element_rect(
    #   fill = "white",
    #   colour = "grey70"
    # ),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  ) +
  guides(
    colour = guide_legend(ncol = 1)
  )
dev.off()

##### USED brain perf unweighted MAE with just cohort info -----

# load meta results 
brain_MAE_mv_global <- readRDS(file.path(results_dir, "perf_brain_MAE_noDunedin_global_final/model_brain_MAE_mv_global.rds"))
brain_MAE_mv_modelwise <- readRDS(file.path(results_dir, "perf_brain_MAE_noDunedin_modelwise_final/model_brain_MAE_mv_modelwise.rds"))


perf_brain_sub_plotting <- perf_brain_sub

cohort_levels <- unique(perf_brain_sub_plotting$cohort)

cohort_palette <- c(
  "ALSPAC" = "#222222",
  "BHRC" = "#f3c300",
  "CannTeen" = "#875692",
  "DCHS" = "#f38400",
  "FinnBrain" = "#a1caf1",
  "GenR" = "#be0032",
  "GUSTO" = "#c2b280",
  "K2H Childhood" = "#848482",
  "K2H Infancy" = "#008856",
  "MTwiNS" = "#e68fac",
  "NICAP" = "#0067a5",
  "Oregon ADHD-1000" = "#f99379",
  "FFCWS" = "#604e97",
  "TAG" = "#f6a600",
  "UCI Echo" = "#b3446c"
)

# extracting meta mean 
overall_est <- data.frame(
  model = "Pooled",
  MAE  = as.numeric(coef(brain_MAE_mv_global)),
  ci.lb = as.numeric(brain_MAE_mv_global$ci.lb),
  ci.ub = as.numeric(brain_MAE_mv_global$ci.ub)
)

# extracting modelwise 
brain_MAE_est <- tibble(
  term  = rownames(brain_MAE_mv_modelwise$b),
  MAE     = as.numeric(brain_MAE_mv_modelwise$b),
  se    = brain_MAE_mv_modelwise$se,
  ci.lb = brain_MAE_mv_modelwise$ci.lb,
  ci.ub = brain_MAE_mv_modelwise$ci.ub,
  pval  = brain_MAE_mv_modelwise$pval
) %>%
  mutate(model = str_remove(term, "^model")) %>%
  select(-term) %>%
  left_join(perf_brain_sub %>% count(model, name = "k"), by = "model") %>%
  arrange(MAE)

brain_MAE_est <- brain_MAE_est %>%
  mutate(label = paste0(model, " (k=", k, ")"))

# --- combine modelwise + overall ---
plot_est <- bind_rows(brain_MAE_est, overall_est)

# ---  ordering ---
# models ordered by generation and MAE,
# overall forced to bottom
# sort models from best to worst (lowest MAE = best)

model_order <- brain_MAE_est %>%
  arrange(MAE) %>%
  pull(model)

# we want worst just above Overall
# so reverse model order
model_order <- rev(model_order)

# Overall must be FIRST level to appear at bottom
model_order <- c("Pooled", model_order)

# apply factor levels
perf_brain_sub_plotting$model <- factor(
  as.character(perf_brain_sub_plotting$model),
  levels = model_order
)

plot_est$model <- factor(
  as.character(plot_est$model),
  levels = model_order
)

# divider above overall row
overall_divider_position <- 1.5

cohort_levels <- unique(perf_brain_sub_plotting$cohort)

png(
  filename = file.path(
    plots_dir,
    "plot_perf_brain_modelwise_plus_overall_meta_raw_combined_cohortcolors_unweightedMAE.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot() +
  
  # raw cohort-level results
  geom_point(
    data = perf_brain_sub_plotting,
    aes(
      x = MAE,
      y = model,
      colour = cohort
    ),
    alpha = .45,
    size = 2,
    position = position_jitter(height = .12)
  ) +
  
  # meta CIs
  geom_segment(
    data = plot_est,
    aes(
      x = ci.lb,
      xend = ci.ub,
      y = model,
      yend = model,
      colour = "Meta-analysis"
    ),
    linewidth = 0.9
  ) +
  
  # meta pooled estimate
  geom_point(
    data = plot_est,
    aes(
      x = MAE,
      y = model,
      colour = "Meta-analysis"
    ),
    shape = 18,
    size = 6
  ) +
  
  
  
  # divider above Overall
  geom_hline(
    yintercept = overall_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  scale_y_discrete(limits = levels(plot_est$model)) +
  
  scale_colour_manual(
    values = c(
      cohort_palette,
      "Meta-analysis" = "black"
    ),
    breaks = c(names(cohort_palette), "Meta-analysis")
  ) +
  
  labs(
    title = "Brain Age Model Performance",
    x = "MAE",
    y = NULL,
    colour = "Source"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    # legend.background = element_rect(
    #   fill = "white",
    #   colour = "grey70"
    # ),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank()
  ) +
  guides(
    colour = guide_legend(ncol = 1)
  )
dev.off()


##### USED epi perf modelwise performance (2A) ----
epi_wMAE_mv_global <- readRDS(file.path(results_dir, "perf_epi_wMAE_noDunedinDNAmTL_global_final/model_epi_wMAE_mv_global.rds"))
epi_wMAE_mv_modelwise <- readRDS(file.path(results_dir, "perf_epi_wMAE_noDunedinDNAmTL_modelwise_final/model_epi_wMAE_mv_modelwise.rds"))

perf_epi_plotting <- perf_epi_sub

cohort_palette <- c(
  "ALSPAC" = "#222222",
  "BHRC" = "#f3c300",
  "CannTeen" = "#875692",
  "DCHS" = "#f38400",
  "FinnBrain" = "#a1caf1",
  "GenR" = "#be0032",
  "GUSTO" = "#c2b280",
  "K2H Childhood" = "#848482",
  "K2H Infancy" = "#008856",
  "MTwiNS" = "#e68fac",
  "NICAP" = "#0067a5",
  "Oregon ADHD-1000" = "#f99379",
  "FFCWS" = "#604e97",
  "TAG" = "#f6a600",
  "UCI Echo" = "#b3446c"
)

# extract global pooled estimate 
overall_est <- data.frame(
  model      = "Pooled",
  wMAE       = exp(epi_wMAE_mv_global$beta[1, 1]),
  ci.lb      = exp(epi_wMAE_mv_global$ci.lb[1]),
  ci.ub      = exp(epi_wMAE_mv_global$ci.ub[1]),
  generation = "Pooled"
)

# extract modelwise 
epi_wMAE_est <- data.frame(
  model = gsub("^model", "", rownames(epi_wMAE_mv_modelwise$beta)),
  wMAE  = exp(epi_wMAE_mv_modelwise$beta[, 1]), # geometric mean
  ci.lb = exp(epi_wMAE_mv_modelwise$ci.lb),
  ci.ub = exp(epi_wMAE_mv_modelwise$ci.ub)
)


# add generation var 
epi_wMAE_est <- epi_wMAE_est %>%
  mutate(
    generation = case_when(
      model %in% c(
        "Bohlin",
        "Knight",
        "EPIC"
      ) ~ "Gen1-gest",
      model %in% c(
        "ZhangEN",
        "ZhangBLUP",
        "Wu",
        "skinHorvath",
        "PedBE",
        "PCBrainAge",
        "Horvath2013",
        "Hannum",
        "CorticalClock",
        "cAge",
        "AltumAge"
      ) ~ "Gen1",
      
      model %in% c(
        "PCGrimAge",
        "PhenoAge",
        "DamAge",
        "AdaptAge"
      ) ~ "Gen2-4"
    )
  )


# order models within generation by performance 
gen1_gest_order <- epi_wMAE_est %>%
  filter(generation == "Gen1-gest") %>%
  arrange(wMAE) %>%
  pull(model)

gen1_order <- epi_wMAE_est %>%
  filter(generation == "Gen1") %>%
  arrange(wMAE) %>%
  pull(model)

gen2_order <- epi_wMAE_est %>%
  filter(generation == "Gen2-4") %>%
  arrange(wMAE) %>%
  pull(model)

# create final y axis order 
model_order <- c(
  "Pooled",
  rev(gen2_order),
  rev(gen1_order),
  rev(gen1_gest_order)
)

# build final ploitting table (combination of modelwise + pooled) 
plot_est <- bind_rows(
  epi_wMAE_est,
  overall_est
)

# apply factor ordering 
perf_epi_plotting$model <- factor(
  as.character(perf_epi_plotting$model),
  levels = model_order
)

plot_est$model <- factor(
  as.character(plot_est$model),
  levels = model_order
)

# divider positions 
overall_divider_position <- 1.5

# divider between Gen2-4 and Gen1
divider_position_gen2_gen1 <- 1.5 + length(gen2_order)

# divider between Gen1 and Gen1-gest
divider_position_gen1_gest <- 1.5 + length(gen2_order) + length(gen1_order)

cohort_levels <- unique(perf_epi_plotting$cohort)

# create plot
png(
  filename = file.path(
    plots_dir,
    "plot_perf_epi_modelwise_plus_overall_meta_raw_combined_cohortcolors.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot() +
  
  geom_point(
    data = perf_epi_plotting,
    aes(
      x = wMAE_test,
      y = model,
      colour = cohort
    ),
    alpha = .45,
    size = 2,
    position = position_jitter(height = .12)
  ) +
  
  geom_segment(
    data = plot_est,
    aes(
      x = ci.lb,
      xend = ci.ub,
      y = model,
      yend = model,
      colour = "Meta-analysis"
    ),
    linewidth = 0.9
  ) +
  
  geom_point(
    data = plot_est,
    aes(
      x = wMAE,
      y = model,
      colour = "Meta-analysis"
    ),
    shape = 18,
    size = 5
  ) +
  
  geom_hline(
    yintercept = divider_position_gen2_gen1,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = divider_position_gen1_gest,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = overall_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  scale_y_discrete(limits = levels(plot_est$model)) +
  
  scale_colour_manual(
    values = c(
      cohort_palette,
      "Meta-analysis" = "black"
    ),
    breaks = c(names(cohort_palette), "Meta-analysis")
  ) +
  
  labs(
    title = "Epigenetic clock performance (weighted MAE)",
    x = "Weighted MAE (MAE ÷ age range)",
    y = NULL,
    colour = "Source"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
  )) +
  guides(
    colour = guide_legend(ncol = 1)
  )

dev.off()

##### epi perf modelwise log -----
epi_wMAE_mv_global <- readRDS(file.path(results_dir, "perf_epi_wMAE_noDunedinDNAmTL_global_final/model_epi_wMAE_mv_global.rds"))
epi_wMAE_mv_modelwise <- readRDS(file.path(results_dir, "perf_epi_wMAE_noDunedinDNAmTL_modelwise_final/model_epi_wMAE_mv_modelwise.rds"))

perf_epi_plotting <- perf_epi_sub

cohort_palette <- c(
  "ALSPAC" = "#222222",
  "BHRC" = "#f3c300",
  "CannTeen" = "#875692",
  "DCHS" = "#f38400",
  "FinnBrain" = "#a1caf1",
  "GenR" = "#be0032",
  "GUSTO" = "#c2b280",
  "K2H Childhood" = "#848482",
  "K2H Infancy" = "#008856",
  "MTwiNS" = "#e68fac",
  "NICAP" = "#0067a5",
  "Oregon ADHD-1000" = "#f99379",
  "FFCWS" = "#604e97",
  "TAG" = "#f6a600",
  "UCI Echo" = "#b3446c"
)

# extract global pooled estimate 
overall_est <- data.frame(
  model      = "Pooled",
  wMAE       = exp(epi_wMAE_mv_global$beta[1, 1]),
  ci.lb      = exp(epi_wMAE_mv_global$ci.lb[1]),
  ci.ub      = exp(epi_wMAE_mv_global$ci.ub[1]),
  generation = "Pooled"
)

# extract modelwise 
epi_wMAE_est <- data.frame(
  model = gsub("^model", "", rownames(epi_wMAE_mv_modelwise$beta)),
  wMAE  = exp(epi_wMAE_mv_modelwise$beta[, 1]), # geometric mean
  ci.lb = exp(epi_wMAE_mv_modelwise$ci.lb),
  ci.ub = exp(epi_wMAE_mv_modelwise$ci.ub)
)


# add generation var 
epi_wMAE_est <- epi_wMAE_est %>%
  mutate(
    generation = case_when(
      model %in% c(
        "Bohlin",
        "Knight",
        "EPIC"
      ) ~ "Gen1-gest",
      model %in% c(
        "ZhangEN",
        "ZhangBLUP",
        "Wu",
        "skinHorvath",
        "PedBE",
        "PCBrainAge",
        "Horvath2013",
        "Hannum",
        "CorticalClock",
        "cAge",
        "AltumAge"
      ) ~ "Gen1",
      
      model %in% c(
        "PCGrimAge",
        "PhenoAge",
        "DamAge",
        "AdaptAge"
      ) ~ "Gen2-4"
    )
  )


# order models within generation by performance 
gen1_gest_order <- epi_wMAE_est %>%
  filter(generation == "Gen1-gest") %>%
  arrange(wMAE) %>%
  pull(model)

gen1_order <- epi_wMAE_est %>%
  filter(generation == "Gen1") %>%
  arrange(wMAE) %>%
  pull(model)

gen2_order <- epi_wMAE_est %>%
  filter(generation == "Gen2-4") %>%
  arrange(wMAE) %>%
  pull(model)

# create final y axis order 
model_order <- c(
  "Pooled",
  rev(gen2_order),
  rev(gen1_order),
  rev(gen1_gest_order)
)

# build final ploitting table (combination of modelwise + pooled) 
plot_est <- bind_rows(
  epi_wMAE_est,
  overall_est
)

# apply factor ordering 
perf_epi_plotting$model <- factor(
  as.character(perf_epi_plotting$model),
  levels = model_order
)

plot_est$model <- factor(
  as.character(plot_est$model),
  levels = model_order
)

# divider positions 
overall_divider_position <- 1.5

# divider between Gen2-4 and Gen1
divider_position_gen2_gen1 <- 1.5 + length(gen2_order)

# divider between Gen1 and Gen1-gest
divider_position_gen1_gest <- 1.5 + length(gen2_order) + length(gen1_order)

cohort_levels <- unique(perf_epi_plotting$cohort)

# create plot
png(
  filename = file.path(
    plots_dir,
    "plot_perf_epi_modelwise_plus_overall_meta_raw_combined_cohortcolors_log.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot() +
  
  geom_point(
    data = perf_epi_plotting,
    aes(
      x = wMAE_test,
      y = model,
      colour = cohort
    ),
    alpha = .45,
    size = 2,
    position = position_jitter(height = .12)
  ) +
  
  geom_segment(
    data = plot_est,
    aes(
      x = ci.lb,
      xend = ci.ub,
      y = model,
      yend = model,
      colour = "Meta-analysis"
    ),
    linewidth = 0.9
  ) +
  
  geom_point(
    data = plot_est,
    aes(
      x = wMAE,
      y = model,
      colour = "Meta-analysis"
    ),
    shape = 18,
    size = 6
  ) +
  
  geom_hline(
    yintercept = divider_position_gen2_gen1,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = divider_position_gen1_gest,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = overall_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  scale_y_discrete(limits = levels(plot_est$model)) +
  
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 20, 50, 100),
    labels = c("1", "2", "5", "10", "20", "50", "100")
  ) +
  
  scale_colour_manual(
    values = c(
      cohort_palette,
      "Meta-analysis" = "black"
    ),
    breaks = c(names(cohort_palette), "Meta-analysis")
  ) +
  
  labs(
    title = "Epigenetic Clock Performance",
    x = "Weighted MAE (MAE ÷ age range), log scale",
    y = NULL,
    colour = "Source"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major.y = element_blank()
  ) +
  guides(
    colour = guide_legend(ncol = 1)
  )

dev.off()

#####  USED epi unweighted MAEmodelwise performance ----
epi_MAE_mv_global <- readRDS(file.path(results_dir, "perf_epi_MAE_noDunedinDNAmTL_global_BFGSoptimiser_final/model_epi_MAE_mv_global.rds"))
epi_MAE_mv_modelwise <- readRDS(file.path(results_dir, "perf_epi_MAE_noDunedinDNAmTL_modelwise_BFGSoptimiser_final/model_epi_MAE_mv_modelwise.rds"))

perf_epi_plotting <- perf_epi_sub

cohort_palette <- c(
  "ALSPAC" = "#222222",
  "BHRC" = "#f3c300",
  "CannTeen" = "#875692",
  "DCHS" = "#f38400",
  "FinnBrain" = "#a1caf1",
  "GenR" = "#be0032",
  "GUSTO" = "#c2b280",
  "K2H Childhood" = "#848482",
  "K2H Infancy" = "#008856",
  "MTwiNS" = "#e68fac",
  "NICAP" = "#0067a5",
  "Oregon ADHD-1000" = "#f99379",
  "FFCWS" = "#604e97",
  "TAG" = "#f6a600",
  "UCI Echo" = "#b3446c"
)

# extract global pooled estimate 
overall_est <- data.frame(
  model      = "Pooled",
  MAE       = exp(epi_MAE_mv_global$beta[1, 1]),
  ci.lb      = exp(epi_MAE_mv_global$ci.lb[1]),
  ci.ub      = exp(epi_MAE_mv_global$ci.ub[1]),
  generation = "Pooled"
)

# extract modelwise 
epi_MAE_est <- data.frame(
  model = gsub("^model", "", rownames(epi_MAE_mv_modelwise$beta)),
  MAE  = exp(epi_MAE_mv_modelwise$beta[, 1]), # geometric mean
  ci.lb = exp(epi_MAE_mv_modelwise$ci.lb),
  ci.ub = exp(epi_MAE_mv_modelwise$ci.ub)
)


# add generation var 
epi_MAE_est <- epi_MAE_est %>%
  mutate(
    generation = case_when(
      model %in% c(
        "Bohlin",
        "Knight",
        "EPIC"
      ) ~ "Gen1-gest",
      model %in% c(
        "ZhangEN",
        "ZhangBLUP",
        "Wu",
        "skinHorvath",
        "PedBE",
        "PCBrainAge",
        "Horvath2013",
        "Hannum",
        "CorticalClock",
        "cAge",
        "AltumAge"
      ) ~ "Gen1",
      
      model %in% c(
        "PCGrimAge",
        "PhenoAge",
        "DamAge",
        "AdaptAge"
      ) ~ "Gen2-4"
    )
  )




# order models within generation by performance 
# order models within generation by performance 
gen1_gest_order <- epi_MAE_est %>%
  filter(generation == "Gen1-gest") %>%
  arrange(MAE) %>%
  pull(model)

gen1_order <- epi_MAE_est %>%
  filter(generation == "Gen1") %>%
  arrange(MAE) %>%
  pull(model)

gen2_order <- epi_MAE_est %>%
  filter(generation == "Gen2-4") %>%
  arrange(MAE) %>%
  pull(model)

# create final y axis order 
model_order <- c(
  "Pooled",
  rev(gen2_order),
  rev(gen1_order),
  rev(gen1_gest_order)
)

# build final ploitting table (combination of modelwise + pooled) 
plot_est <- bind_rows(
  epi_MAE_est,
  overall_est
)

# apply factor ordering 
perf_epi_plotting$model <- factor(
  as.character(perf_epi_plotting$model),
  levels = model_order
)

plot_est$model <- factor(
  as.character(plot_est$model),
  levels = model_order
)

# divider positions 
overall_divider_position <- 1.5

# divider between Gen2-4 and Gen1
divider_position_gen2_gen1 <- 1.5 + length(gen2_order)

# divider between Gen1 and Gen1-gest
divider_position_gen1_gest <- 1.5 + length(gen2_order) + length(gen1_order)

cohort_levels <- unique(perf_epi_plotting$cohort)

# create plot
png(
  filename = file.path(
    plots_dir,
    "plot_perf_epi_unweightedMAE_modelwise_plus_overall_meta_raw_combined_cohortcolors.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot() +
  
  geom_point(
    data = perf_epi_plotting,
    aes(
      x = MAE,
      y = model,
      colour = cohort
    ),
    alpha = .45,
    size = 2,
    position = position_jitter(height = .12)
  ) +
  
  geom_segment(
    data = plot_est,
    aes(
      x = ci.lb,
      xend = ci.ub,
      y = model,
      yend = model,
      colour = "Meta-analysis"
    ),
    linewidth = 0.9
  ) +
  
  geom_point(
    data = plot_est,
    aes(
      x = MAE,
      y = model,
      colour = "Meta-analysis"
    ),
    shape = 18,
    size = 6
  ) +
  
  geom_hline(
    yintercept = divider_position_gen2_gen1,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = divider_position_gen1_gest,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = overall_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  scale_y_discrete(limits = levels(plot_est$model)) +
  
  scale_colour_manual(
    values = c(
      cohort_palette,
      "Meta-analysis" = "black"
    ),
    breaks = c(names(cohort_palette), "Meta-analysis")
  ) +
  
  labs(
    title = "Epigenetic Clock Performance",
    x = "MAE",
    y = NULL,
    colour = "Source"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major.y = element_blank()
  ) +
  guides(
    colour = guide_legend(ncol = 1)
  )

dev.off()

#### correlation ----

##### USED violin boxplot epi Pearson r by age bin (SM Figure S1) -----
# includes model-specific meta estimates with robust fallback


# convergence issues in the subgroup meta-analyses so asked chatchpt to code a function for fallback options - first option isthe normal model, then we go to a different optimiser (the one also used for unweighted MAE for epi), then a RE model without multilevel structure, then a fixed effect model --> it's not important for the plot that all fo this is perfect 
fit_meta_cor <- function(dat) {
  
  # preferred: same multilevel model as main analysis
  res <- tryCatch(
    rma.mv(
      yi = pearson_z,
      V = pearson_var,
      random = list(~ 1 | cohort/timepoint),
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma.mv REML"))
  }
  
  # same model, alternative optimizer
  res <- tryCatch(
    rma.mv(
      yi = pearson_z,
      V = pearson_var,
      random = list(~ 1 | cohort/timepoint),
      data = dat,
      method = "REML",
      control = list(
        optimizer = "optim",
        optmethod = "BFGS",
        maxit = 10000
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma.mv REML BFGS"))
  }
  
  # simpler random-effects fallback
  res <- tryCatch(
    rma(
      yi = pearson_z,
      vi = pearson_var,
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma REML fallback"))
  }
  
  # final fixed-effect fallback
  res <- rma(
    yi = pearson_z,
    vi = pearson_var,
    data = dat,
    method = "FE"
  )
  
  list(fit = res, meta_model = "fixed-effect fallback")
}

model_palette <- c(
  "PCBrainAge" = "#222222",
  "PCGrimAge" = "#f3c300",
  "CorticalClock" = "#875692",
  "cAge" = "#8db600",
  "AltumAge" = "#a1caf1",
  "Horvath2013" = "#be0032",
  "DamAge" = "#c2b280",
  "AdaptAge" = "#848482",
  "Hannum" = "#008856",
  "PhenoAge" = "#e68fac",
  "skinHorvath" = "#0067a5",
  "PedBE" = "#f99379",
  "Wu" = "#604e97",
  "ZhangBLUP" = "#f6a600",
  "ZhangEN" = "#b3446c",
  "Bohlin" = "#dcd300",
  "Knight" = "#654522",
  "EPIC" = "#2b3d26",
  "DunedinPACE" = "#e25822"
)


perf_epi <- perf_epi %>%
  mutate(
    age_bin = case_when(
      mean_age < 0.25 ~ "birth\n<3m",
      mean_age < 3    ~ "Infancy\n<3y",
      mean_age < 6    ~ "Early childhood\n3–5y",
      mean_age < 10   ~ "Middle childhood\n6–9y",
      mean_age < 15   ~ "Late childhood\n10–14y",
      mean_age < 18   ~ "Adolescence\n15–17y",
      TRUE            ~ "Young adulthood\n18+y"
    ),
    age_bin = factor(
      age_bin,
      levels = c(
        "birth\n<3m",
        "Infancy\n<3y",
        "Early childhood\n3–5y",
        "Middle childhood\n6–9y",
        "Late childhood\n10–14y",
        "Adolescence\n15–17y",
        "Young adulthood\n18+y"
      )
    )
  )

epi_age_bin_counts <- perf_epi %>%
  group_by(age_bin) %>%
  summarise(
    k = n_distinct(paste(cohort, timepoint)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("k=", k))


meta_results_r <- perf_epi %>%
  filter(
    !is.na(pearson_z),
    !is.na(pearson_var),
    pearson_var > 0
  ) %>%
  group_by(age_bin, model) %>%
  group_modify(~{
    
    dat <- .x
    
    if (nrow(dat) == 1) {
      return(tibble(
        pooled_z = dat$pearson_z[1],
        ci_lb_z = NA_real_,
        ci_ub_z = NA_real_,
        pooled_r = transf.ztor(dat$pearson_z[1]),
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        k = 1,
        meta_model = "single estimate"
      ))
    }
    
    fit_out <- fit_meta_cor(dat)
    res <- fit_out$fit
    
    tibble(
      pooled_z = as.numeric(coef(res)[1]),
      ci_lb_z = res$ci.lb[1],
      ci_ub_z = res$ci.ub[1],
      pooled_r = transf.ztor(as.numeric(coef(res)[1])),
      ci_lb = transf.ztor(res$ci.lb[1]),
      ci_ub = transf.ztor(res$ci.ub[1]),
      k = res$k,
      meta_model = fit_out$meta_model
    )
  }) %>%
  ungroup()

# check fallback use
meta_results_r %>%
  count(meta_model)


pd <- position_dodge(width = 0.7)

y_min <- min(
  perf_epi$Pearson,
  meta_results_r$ci_lb,
  meta_results_r$pooled_r,
  na.rm = TRUE
)

y_max <- max(
  perf_epi$Pearson,
  meta_results_r$ci_ub,
  meta_results_r$pooled_r,
  na.rm = TRUE
)

y_pad <- 0.08 * (y_max - y_min)


png(
  filename = file.path(
    plots_dir,
    "plot_violin_perf_epi_age_pearson_brackets.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot(perf_epi, aes(x = age_bin, y = Pearson)) +
  
  geom_hline(
    yintercept = 0,
    colour = "grey92",
    linewidth = 0.5
  ) +
  
  geom_violin(
    fill = "grey85",
    colour = NA,
    trim = FALSE
  ) +
  
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = NA
  ) +
  
  geom_errorbar(
    data = meta_results_r,
    aes(
      x = age_bin,
      ymin = ci_lb,
      ymax = ci_ub,
      colour = model,
      group = model
    ),
    position = pd,
    width = 0.08,
    alpha = 0.45,
    na.rm = TRUE,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = meta_results_r,
    aes(
      x = age_bin,
      y = pooled_r,
      colour = model,
      group = model
    ),
    position = pd,
    alpha = 0.6,
    size = 2.4,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(values = model_palette) +
  
  geom_text(
    data = epi_age_bin_counts,
    aes(
      x = age_bin,
      y = y_min - 0.35 * y_pad,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = 1,
    size = 3.5
  ) +
  
  coord_cartesian(
    ylim = c(
      y_min - y_pad,
      y_max + y_pad
    ),
    clip = "off"
  ) +
  
  labs(
    x = "Developmental age group",
    y = "Pearson r",
    title = "Epigenetic clock correlations across age groups",
    colour = "Model"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  )

dev.off()

##### USED epi Pearson modelwise performance (2B) ----

epi_pearson_mv_global<- readRDS(
  file.path(
    results_dir,
    "perf_epi_pearson_global_final/model_epi_pearson_mv_global.rds"
  )
)

epi_pearson_mv_modelwise <- readRDS(
  file.path(
    results_dir,
    "perf_epi_pearson_modelwise_final/model_epi_pearson_mv_modelwise.rds"
  )
)

perf_epi_plotting <- perf_epi

cohort_palette <- c(
  "ALSPAC" = "#222222",
  "BHRC" = "#f3c300",
  "CannTeen" = "#875692",
  "DCHS" = "#f38400",
  "FinnBrain" = "#a1caf1",
  "GenR" = "#be0032",
  "GUSTO" = "#c2b280",
  "K2H Childhood" = "#848482",
  "K2H Infancy" = "#008856",
  "MTwiNS" = "#e68fac",
  "NICAP" = "#0067a5",
  "Oregon ADHD-1000" = "#f99379",
  "FFCWS" = "#604e97",
  "TAG" = "#f6a600",
  "UCI Echo" = "#b3446c"
)

# --- extract global pooled estimate ---

overall_est <- data.frame(
  model = "Pooled",
  pearson_z = epi_pearson_mv_global$beta[1, 1],
  pearson_r = transf.ztor(epi_pearson_mv_global$beta[1, 1]),
  ci.lb = transf.ztor(epi_pearson_mv_global$ci.lb[1]),
  ci.ub = transf.ztor(epi_pearson_mv_global$ci.ub[1]),
  generation = "Pooled"
)

# --- extract modelwise estimates ---

epi_pearson_est <- data.frame(
  model = gsub("^model", "", rownames(epi_pearson_mv_modelwise$beta)),
  pearson_z = epi_pearson_mv_modelwise$beta[, 1],
  pearson_r = transf.ztor(epi_pearson_mv_modelwise$beta[, 1]),
  ci.lb = transf.ztor(epi_pearson_mv_modelwise$ci.lb),
  ci.ub = transf.ztor(epi_pearson_mv_modelwise$ci.ub)
)

# --- add generation variable ---

epi_pearson_est <- epi_pearson_est %>%
  mutate(
    generation = case_when(
      model %in% c(
        "Bohlin",
        "Knight",
        "EPIC"
      ) ~ "Gen1-gest",
      
      model %in% c(
        "ZhangEN",
        "ZhangBLUP",
        "Wu",
        "skinHorvath",
        "PedBE",
        "PCBrainAge",
        "Horvath2013",
        "Hannum",
        "CorticalClock",
        "cAge",
        "AltumAge"
      ) ~ "Gen1",
      
      model %in% c(
        "PCGrimAge",
        "PhenoAge",
        "DamAge",
        "AdaptAge",
        "DunedinPACE",
        "DNAmTL"
      ) ~ "Gen2-4",
    )
  )

# --- order models within generation by correlation ---
# higher r = better performance, so order by pearson_r

gen1_gest_order <- epi_pearson_est %>%
  filter(generation == "Gen1-gest") %>%
  arrange(pearson_r) %>%
  pull(model)

gen1_order <- epi_pearson_est %>%
  filter(generation == "Gen1") %>%
  arrange(pearson_r) %>%
  pull(model)

gen2_order <- epi_pearson_est %>%
  filter(generation == "Gen2-4") %>%
  arrange(pearson_r) %>%
  pull(model)

other_order <- epi_pearson_est %>%
  filter(generation == "Other") %>%
  arrange(pearson_r) %>%
  pull(model)

model_order <- c(
  "Pooled",
  other_order,
  gen2_order,
  gen1_order,
  gen1_gest_order
)

# --- build plotting table ---

plot_est <- bind_rows(
  epi_pearson_est,
  overall_est
)

perf_epi_plotting <- perf_epi_plotting %>%
  mutate(
    model = factor(as.character(model), levels = model_order)
  )

plot_est <- plot_est %>%
  mutate(
    model = factor(as.character(model), levels = model_order)
  )

# --- divider positions ---

overall_divider_position <- 1.5
divider_position_other_gen2 <- 1.5 + length(other_order)
divider_position_gen2_gen1 <- 1.5 + length(other_order) + length(gen2_order)
divider_position_gen1_gest <- 1.5 + length(other_order) + length(gen2_order) + length(gen1_order)

# --- plot ---

png(
  filename = file.path(
    plots_dir,
    "plot_perf_epi_modelwise_plus_overall_pearson_combined_cohortcolors.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot() +
  
  geom_vline(
    xintercept = 0,
    colour = "grey85",
    linewidth = 0.5
  ) +
  
  geom_point(
    data = perf_epi_plotting,
    aes(
      x = Pearson,
      y = model,
      colour = cohort
    ),
    alpha = .45,
    size = 2,
    position = position_jitter(height = .12)
  ) +
  
  geom_segment(
    data = plot_est,
    aes(
      x = ci.lb,
      xend = ci.ub,
      y = model,
      yend = model,
      colour = "Overall pooled"
    ),
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = plot_est,
    aes(
      x = pearson_r,
      y = model,
      colour = "Overall pooled"
    ),
    shape = 18,
    size = 5,
    inherit.aes = FALSE
  ) +
  
  geom_hline(
    yintercept = overall_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  {
    if (length(other_order) > 0) {
      geom_hline(
        yintercept = divider_position_other_gen2,
        linetype = "dashed",
        linewidth = 0.4,
        colour = "grey40"
      )
    }
  } +
  
  geom_hline(
    yintercept = divider_position_gen2_gen1,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = divider_position_gen1_gest,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  scale_y_discrete(limits = model_order) +
  
  scale_colour_manual(
    values = c(
      cohort_palette,
      "Overall pooled" = "black"
    ),
    breaks = c(names(cohort_palette), "Overall pooled")
  ) +
  
  labs(
    title = "Epigenetic clock performance (correlation)",
    x = "Pearson correlation between epigenetic age and chronological age",
    y = NULL,
    colour = "Source"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  ) +
  
  guides(
    colour = guide_legend(ncol = 1)
  )

dev.off()

##### USED violin boxplot brain Pearson r by age bin -----
# includes model-specific meta estimates with robust fallback

fit_meta_cor <- function(dat) {
  
  res <- tryCatch(
    rma.mv(
      yi = pearson_z,
      V = pearson_var,
      random = list(~ 1 | cohort/timepoint),
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma.mv REML"))
  }
  
  res <- tryCatch(
    rma.mv(
      yi = pearson_z,
      V = pearson_var,
      random = list(~ 1 | cohort/timepoint),
      data = dat,
      method = "REML",
      control = list(
        optimizer = "optim",
        optmethod = "BFGS",
        maxit = 10000
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma.mv REML BFGS"))
  }
  
  res <- tryCatch(
    rma(
      yi = pearson_z,
      vi = pearson_var,
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma REML fallback"))
  }
  
  res <- rma(
    yi = pearson_z,
    vi = pearson_var,
    data = dat,
    method = "FE"
  )
  
  list(fit = res, meta_model = "fixed-effect fallback")
}


model_palette <- c(
  "Pyment" = "#222222",
  "DevBrainAge" = "#f3c300",
  "Centile2" = "#875692",
  "DBN" = "#8db600",
  "Kaufmann" = "#a1caf1",
  "PyBrainAge" = "#be0032",
  "ENIGMA" = "#c2b280",
  "DunedinPACNI" = "#e25822"
)

perf_brain <- perf_brain %>%
  mutate(
    age_bin = case_when(
      mean_age < 6  ~ "Early childhood\n3–5y",
      mean_age < 10 ~ "Middle childhood\n6–9y",
      mean_age < 15 ~ "Late childhood\n10–14y",
      mean_age < 18 ~ "Adolescence\n15–17y",
      TRUE          ~ "Young adulthood\n18+y"
    ),
    age_bin = factor(
      age_bin,
      levels = c(
        "Early childhood\n3–5y",
        "Middle childhood\n6–9y",
        "Late childhood\n10–14y",
        "Adolescence\n15–17y",
        "Young adulthood\n18+y"
      )
    )
  )

brain_age_bin_counts <- perf_brain %>%
  group_by(age_bin) %>%
  summarise(
    k = n_distinct(paste(cohort, timepoint)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("k=", k))

meta_results_r <- perf_brain %>%
  filter(
    !is.na(pearson_z),
    !is.na(pearson_var),
    pearson_var > 0
  ) %>%
  group_by(age_bin, model) %>%
  group_modify(~{
    
    dat <- .x
    
    if (nrow(dat) == 1) {
      return(tibble(
        pooled_z = dat$pearson_z[1],
        ci_lb_z = NA_real_,
        ci_ub_z = NA_real_,
        pooled_r = transf.ztor(dat$pearson_z[1]),
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        k = 1,
        meta_model = "single estimate"
      ))
    }
    
    fit_out <- fit_meta_cor(dat)
    res <- fit_out$fit
    
    tibble(
      pooled_z = as.numeric(coef(res)[1]),
      ci_lb_z = res$ci.lb[1],
      ci_ub_z = res$ci.ub[1],
      pooled_r = transf.ztor(as.numeric(coef(res)[1])),
      ci_lb = transf.ztor(res$ci.lb[1]),
      ci_ub = transf.ztor(res$ci.ub[1]),
      k = res$k,
      meta_model = fit_out$meta_model
    )
  }) %>%
  ungroup()

meta_results_r %>%
  count(meta_model)

pd <- position_dodge(width = 0.7)

y_min <- min(
  perf_brain$Pearson,
  meta_results_r$ci_lb,
  meta_results_r$pooled_r,
  na.rm = TRUE
)

y_max <- max(
  perf_brain$Pearson,
  meta_results_r$ci_ub,
  meta_results_r$pooled_r,
  na.rm = TRUE
)

y_pad <- 0.08 * (y_max - y_min)

png(
  filename = file.path(
    plots_dir,
    "plot_violin_perf_brain_age_pearson_brackets.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot(perf_brain, aes(x = age_bin, y = Pearson)) +
  
  geom_hline(
    yintercept = 0,
    colour = "grey92",
    linewidth = 0.5
  ) +
  
  geom_violin(
    fill = "grey85",
    colour = NA,
    trim = FALSE
  ) +
  
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = NA
  ) +
  
  geom_errorbar(
    data = meta_results_r,
    aes(
      x = age_bin,
      ymin = ci_lb,
      ymax = ci_ub,
      colour = model,
      group = model
    ),
    position = pd,
    width = 0.08,
    alpha = 0.45,
    na.rm = TRUE,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = meta_results_r,
    aes(
      x = age_bin,
      y = pooled_r,
      colour = model,
      group = model
    ),
    position = pd,
    alpha = 0.6,
    size = 2.4,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(values = model_palette) +
  
  geom_text(
    data = brain_age_bin_counts,
    aes(
      x = age_bin,
      y = y_min - 0.35 * y_pad,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = 1,
    size = 3.5
  ) +
  
  coord_cartesian(
    ylim = c(
      y_min - y_pad,
      y_max + y_pad
    ),
    clip = "off"
  ) +
  
  labs(
    x = "Developmental age group",
    y = "Pearson r",
    title = "Brain age model correlations across age groups",
    colour = "Model"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  )

dev.off()

##### USED brain Pearson modelwise performance (3B) ----

brain_pearson_mv_global <- readRDS(
  file.path(
    results_dir,
    "perf_brain_pearson_global_final/model_brain_pearson_mv_global.rds"
  )
)

brain_pearson_mv_modelwise <- readRDS(
  file.path(
    results_dir,
    "perf_brain_pearson_modelwise_final/model_brain_pearson_mv_modelwise.rds"
  )
)

perf_brain_plotting <- perf_brain

cohort_palette <- c(
  "ALSPAC" = "#222222",
  "BHRC" = "#f3c300",
  "CannTeen" = "#875692",
  "DCHS" = "#f38400",
  "FinnBrain" = "#a1caf1",
  "GenR" = "#be0032",
  "GUSTO" = "#c2b280",
  "K2H Childhood" = "#848482",
  "K2H Infancy" = "#008856",
  "MTwiNS" = "#e68fac",
  "NICAP" = "#0067a5",
  "Oregon ADHD-1000" = "#f99379",
  "FFCWS" = "#604e97",
  "TAG" = "#f6a600",
  "UCI Echo" = "#b3446c"
)

overall_est <- data.frame(
  model = "Pooled",
  pearson_z = brain_pearson_mv_global$beta[1, 1],
  pearson_r = transf.ztor(brain_pearson_mv_global$beta[1, 1]),
  ci.lb = transf.ztor(brain_pearson_mv_global$ci.lb[1]),
  ci.ub = transf.ztor(brain_pearson_mv_global$ci.ub[1]),
  generation = "Pooled"
)

brain_pearson_est <- data.frame(
  model = gsub("^model", "", rownames(brain_pearson_mv_modelwise$beta)),
  pearson_z = brain_pearson_mv_modelwise$beta[, 1],
  pearson_r = transf.ztor(brain_pearson_mv_modelwise$beta[, 1]),
  ci.lb = transf.ztor(brain_pearson_mv_modelwise$ci.lb),
  ci.ub = transf.ztor(brain_pearson_mv_modelwise$ci.ub)
)

brain_pearson_est <- brain_pearson_est %>%
  mutate(
    generation = case_when(
      model %in% c(
        "Pyment",
        "DevBrainAge",
        "Centile2",
        "DBN",
        "Kaufmann",
        "PyBrainAge",
        "ENIGMA"
      ) ~ "Brain age",
      
      model %in% c("DunedinPACNI") ~ "Next generation"
    )
  )

brain_age_order <- brain_pearson_est %>%
  filter(generation == "Brain age") %>%
  arrange(pearson_r) %>%
  pull(model)

nextgen_order <- brain_pearson_est %>%
  filter(generation == "Next generation") %>%
  arrange(pearson_r) %>%
  pull(model)

model_order <- c(
  "Pooled",
  nextgen_order,
  brain_age_order
)

plot_est <- bind_rows(
  brain_pearson_est,
  overall_est
)

perf_brain_plotting <- perf_brain_plotting %>%
  mutate(
    model = factor(as.character(model), levels = model_order)
  )

plot_est <- plot_est %>%
  mutate(
    model = factor(as.character(model), levels = model_order)
  )

# sanity checks
setdiff(unique(as.character(perf_brain$model)), model_order)
setdiff(unique(as.character(plot_est$model)), model_order)

overall_divider_position <- 1.5
divider_position_nextgen_brainage <- 1.5 + length(nextgen_order)

png(
  filename = file.path(
    plots_dir,
    "plot_perf_brain_modelwise_plus_overall_pearson_combined_cohortcolors.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot() +
  
  geom_vline(
    xintercept = 0,
    colour = "grey85",
    linewidth = 0.5
  ) +
  
  geom_point(
    data = perf_brain_plotting,
    aes(
      x = Pearson,
      y = model,
      colour = cohort
    ),
    alpha = .45,
    size = 2,
    position = position_jitter(height = .12)
  ) +
  
  geom_segment(
    data = plot_est,
    aes(
      x = ci.lb,
      xend = ci.ub,
      y = model,
      yend = model,
      colour = "Overall pooled"
    ),
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = plot_est,
    aes(
      x = pearson_r,
      y = model,
      colour = "Overall pooled"
    ),
    shape = 18,
    size = 5,
    inherit.aes = FALSE
  ) +
  
  geom_hline(
    yintercept = overall_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = divider_position_nextgen_brainage,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  scale_y_discrete(limits = model_order) +
  
  scale_colour_manual(
    values = c(
      cohort_palette,
      "Overall pooled" = "black"
    ),
    breaks = c(names(cohort_palette), "Overall pooled")
  ) +
  
  labs(
    title = "Brain age model performance (correlation)",
    x = "Pearson correlation between brain age and chronological age",
    y = NULL,
    colour = "Source"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  ) +
  
  guides(
    colour = guide_legend(ncol = 1)
  )

dev.off()

#### association ----
##### USED violin boxplot association beta by brain age bin (4C) -----
# includes pooled meta estimates by age bin

# --- age bins ---

mod2_eff2_age_brain_bins <- mod2_eff2 %>%
  mutate(
    age_bin = case_when(
      mean_age < 6  ~ "Early childhood\n4–5y",
      mean_age < 10 ~ "Middle childhood\n6–9y",
      mean_age < 15 ~ "Late childhood\n10–14y",
      mean_age < 18 ~ "Adolescence\n15–17y",
      TRUE ~ "Young adulthood\n18+y"
    ),
    age_bin = factor(
      age_bin,
      levels = c(
        "Early childhood\n4–5y",
        "Middle childhood\n6–9y",
        "Late childhood\n10–14y",
        "Adolescence\n15–17y",
        "Young adulthood\n18+y"
      )
    )
  )

# --- age-bin counts ---

brain_assoc_age_bin_counts <- mod2_eff2_age_brain_bins %>%
  group_by(age_bin) %>%
  summarise(
    k = n_distinct(paste(cohort, timepoint)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("k=", k))

# We again have the problem that the model doesn't converge for all cells , so need fallback options 

# --- helper: robust meta-analysis fit for age_bin x model_combi cells ---

fit_meta_assoc <- function(dat) {
  
  # 1. preferred model: multilevel REML with default optimiser
  res <- tryCatch(
    rma.mv(
      yi = RLM_Estimate_scaled,
      V  = RLM_SE_coeftest_HC3_scaled^2,
      random = list(
        ~ 1 | cohort/timepoint
      ),
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma.mv REML"))
  }
  
  # 2. same multilevel model, BFGS optimiser
  res <- tryCatch(
    rma.mv(
      yi = RLM_Estimate_scaled,
      V  = RLM_SE_coeftest_HC3_scaled^2,
      random = list(
        ~ 1 | cohort/timepoint
      ),
      data = dat,
      method = "REML",
      control = list(
        optimizer = "optim",
        optmethod = "BFGS",
        maxit = 10000
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma.mv REML BFGS"))
  }
  
  # 3. simpler random-effects fallback
  res <- tryCatch(
    rma(
      yi = RLM_Estimate_scaled,
      vi = RLM_SE_coeftest_HC3_scaled^2,
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma REML fallback"))
  }
  
  # 4. final fixed-effect fallback
  res <- rma(
    yi = RLM_Estimate_scaled,
    vi = RLM_SE_coeftest_HC3_scaled^2,
    data = dat,
    method = "FE"
  )
  
  list(fit = res, meta_model = "fixed-effect fallback")
}

# --- pooled meta-analysis by age bin ---

assoc_agebin_meta_brain <- mod2_eff2_age_brain_bins %>%
  filter(
    !is.na(age_bin),
    !is.na(model_combi),
    !is.na(RLM_Estimate_scaled),
    !is.na(RLM_SE_coeftest_HC3_scaled),
    RLM_SE_coeftest_HC3_scaled > 0
  ) %>%
  group_by(age_bin, model_combi) %>%
  group_modify(~{
    
    dat <- .x
    
    if (nrow(dat) == 1) {
      return(tibble(
        pooled_beta = dat$RLM_Estimate_scaled[1],
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        k = 1,
        meta_model = "single estimate"
      ))
    }
    
    fit_out <- fit_meta_assoc(dat)
    res <- fit_out$fit
    
    tibble(
      pooled_beta = as.numeric(coef(res)[1]),
      ci_lb = res$ci.lb[1],
      ci_ub = res$ci.ub[1],
      k = res$k,
      meta_model = fit_out$meta_model
    )
  }) %>%
  ungroup()

# what got used? 
assoc_agebin_meta_brain %>%
  count(meta_model)

assoc_agebin_meta_brain %>%
  filter(meta_model != "rma.mv REML") %>%
  arrange(age_bin, model_combi) %>%
  print(n = Inf)

# --- y-axis padding based on plotted values + meta CIs ---

pd <- position_dodge(width = 0.7)

y_min <- min(
  mod2_eff2_age_brain_bins$RLM_Estimate_scaled,
  assoc_agebin_meta_brain$ci_lb,
  assoc_agebin_meta_brain$pooled_beta,
  na.rm = TRUE
)

y_max <- max(
  mod2_eff2_age_brain_bins$RLM_Estimate_scaled,
  assoc_agebin_meta_brain$ci_ub,
  assoc_agebin_meta_brain$pooled_beta,
  na.rm = TRUE
)

y_pad <- 0.005 * (y_max - y_min)


model_combis <- sort(unique(assoc_agebin_meta_brain$model_combi))

library(Polychrome)

kelly_base <- kelly.colors(22)

model_palette <- setNames(
  rep(kelly_base, length.out = length(model_combis)),
  model_combis
)

# library(colorspace)
# 
# 
# model_palette <- setNames(
#   qualitative_hcl(
#     n = length(model_combis),
#     palette = "Harmonic"
#   ),
#   model_combis
# )

# numeric x positions
mod2_eff2_age_brain_bins <- mod2_eff2_age_brain_bins %>%
  mutate(age_x = as.numeric(age_bin) - 0.18)

assoc_agebin_meta_brain <- assoc_agebin_meta_brain %>%
  mutate(age_x = as.numeric(age_bin) + 0.18)

pd_meta <- position_dodge(width = 0.45)

png(
  filename = file.path(
    plots_dir,
    "plot_violin_assoc_brain_mean_age_brackets_modelcombi_meta_new.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

ggplot() +
  
  geom_hline(yintercept = 0, colour = "grey92", linewidth = 0.5) +
  
  geom_violin(
    data = mod2_eff2_age_brain_bins,
    aes(x = age_x, y = RLM_Estimate_scaled, group = age_bin),
    fill = "grey85",
    colour = NA,
    trim = FALSE,
    width = 0.28
  ) +
  
  geom_boxplot(
    data = mod2_eff2_age_brain_bins,
    aes(x = age_x, y = RLM_Estimate_scaled, group = age_bin),
    width = 0.08,
    outlier.shape = NA,
    fill = NA
  ) +
  
  geom_errorbar(
    data = assoc_agebin_meta_brain,
    aes(
      x = age_x,
      ymin = ci_lb,
      ymax = ci_ub,
      colour = model_combi,
      group = model_combi
    ),
    position = pd_meta,
    width = 0.03,
    alpha = 0.3,
    linewidth = 0.2,
    na.rm = TRUE
  ) +
  
  geom_point(
    data = assoc_agebin_meta_brain,
    aes(
      x = age_x,
      y = pooled_beta,
      colour = model_combi,
      group = model_combi
    ),
    position = pd_meta,
    alpha = 0.5,
    size = 2.2,
    shape = 16,
    stroke = 0
  ) +
  
  scale_x_continuous(
    breaks = seq_along(levels(mod2_eff2_age_brain_bins$age_bin)),
    labels = levels(mod2_eff2_age_brain_bins$age_bin)
  ) +
  
  scale_colour_manual(values = model_palette) +
  
  geom_text(
    data = brain_assoc_age_bin_counts %>%
      mutate(age_x = as.numeric(age_bin) - 0.18),
    aes(
      x = age_x,
      y = y_min - 0.35 * y_pad,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = 1,
    size = 3.5
  ) +
  
  coord_cartesian(
    ylim = c(y_min - y_pad, y_max + y_pad),
    clip = "off"
  ) +
  
  labs(
    x = "Developmental age group",
    y = "Robust standardized beta",
    title = "Brain-PAR - Epigenetic-PAR associations across age groups",
    colour = "Model combination",
    caption = "colours reflect different brain - epigenetic model pairs"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(
      size = 10,
      hjust = 0.5
    ),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

dev.off()

# checking highest and lowest 

largest_assoc_by_agebin <- assoc_agebin_meta_brain %>%
  filter(!is.na(pooled_beta)) %>%
  group_by(age_bin) %>%
  summarise(
    most_positive_model = model_combi[which.max(pooled_beta)],
    most_positive_beta  = max(pooled_beta, na.rm = TRUE),
    most_negative_model = model_combi[which.min(pooled_beta)],
    most_negative_beta  = min(pooled_beta, na.rm = TRUE),
    .groups = "drop"
  )

print(largest_assoc_by_agebin, width = Inf)

# count contributing effect sizes per cell
effect_size_counts <- mod2_eff2_age_brain_bins %>%
  group_by(age_bin, model_combi) %>%
  summarise(
    n_effect_sizes = n(),
    .groups = "drop"
  )

# merge onto meta results
assoc_agebin_meta_brain2 <- assoc_agebin_meta_brain %>%
  left_join(
    effect_size_counts,
    by = c("age_bin", "model_combi")
  )

# most positive
largest_positive_by_agebin <- assoc_agebin_meta_brain2 %>%
  filter(!is.na(pooled_beta)) %>%
  group_by(age_bin) %>%
  slice_max(pooled_beta, n = 1, with_ties = FALSE) %>%
  mutate(direction = "most positive") %>%
  ungroup()

# most negative
largest_negative_by_agebin <- assoc_agebin_meta_brain2 %>%
  filter(!is.na(pooled_beta)) %>%
  group_by(age_bin) %>%
  slice_min(pooled_beta, n = 1, with_ties = FALSE) %>%
  mutate(direction = "most negative") %>%
  ungroup()

# combine
largest_assoc_by_agebin <- bind_rows(
  largest_positive_by_agebin,
  largest_negative_by_agebin
) %>%
  arrange(age_bin, direction)

# inspect
largest_assoc_by_agebin %>%
  select(
    age_bin,
    direction,
    model_combi,
    pooled_beta,
    k,
    meta_model
  ) %>% print(width = Inf)

##### association modelwise - with legend outside ------

# load results needed 
model_assoc_mv_global<- readRDS(
  file.path(
    results_dir,
    "assoc_model2_PAR_global_final/model_assoc_mv_global.rds"
  )
)

model_assoc_mv_modelwise <- readRDS(
  file.path(
    results_dir,
    "assoc_model2_PAR_modelwise_final/model_assoc_mv_modelwise.rds"
  )
)



# 1. Pairwise pooled estimates

pair_est <- tibble(
  term  = rownames(model_assoc_mv_modelwise$b),
  b     = as.numeric(model_assoc_mv_modelwise$b),
  se    = model_assoc_mv_modelwise$se,
  ci.lb = model_assoc_mv_modelwise$ci.lb,
  ci.ub = model_assoc_mv_modelwise$ci.ub,
  pval  = model_assoc_mv_modelwise$pval
) %>%
  mutate(term = str_remove(term, "^model_combi")) %>%
  separate(term, into = c("brain_model", "epi_model"), sep = "\\.")

pair_raw <- mod2_eff2 %>%
  select(cohort, brain_model, epi_model, RLM_Estimate_scaled)

# 2. Define epi clock groups and desired order

gen1_models <- sort(c(
  "AltumAge", "CorticalClock", "Hannum", "Horvath2013",
  "PCBrainAge", "PedBE", "Wu", "ZhangBLUP", "ZhangEN",
  "cAge", "skinHorvath"
))

gen2plus_models <- sort(c(
  "AdaptAge", "DamAge", "DNAmTL", "DunedinPACE",
  "PCGrimAge", "PhenoAge"
))

# We want TOP -> BOTTOM:
# Gen1 alphabetical, then Gen2plus alphabetical, then Pooled
# For ggplot y-axis, levels must be bottom -> top:
epi_levels <- c(
  "Pooled",
  rev(gen2plus_models),
  rev(gen1_models)
)

# 3. Marginal pooled effects

# pooled across epi clocks for each brain model
mod_brain_margin <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ brain_model - 1,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | epi_model
  ),
  data = mod2_eff2,
  method = "REML"
)

brain_margin_est <- tibble(
  brain_model = rownames(mod_brain_margin$b),
  epi_model   = "Pooled",
  b           = as.numeric(mod_brain_margin$b),
  ci.lb       = mod_brain_margin$ci.lb,
  ci.ub       = mod_brain_margin$ci.ub
) %>%
  mutate(
    brain_model = str_remove(brain_model, "^brain_model")
  )

# pooled across brain models for each epi model
mod_epi_margin <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ epi_model - 1,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | brain_model
  ),
  data = mod2_eff2,
  method = "REML"
)

epi_margin_est <- tibble(
  epi_model   = rownames(mod_epi_margin$b),
  brain_model = "Pooled",
  b           = as.numeric(mod_epi_margin$b),
  ci.lb       = mod_epi_margin$ci.lb,
  ci.ub       = mod_epi_margin$ci.ub
) %>%
  mutate(
    epi_model = str_remove(epi_model, "^epi_model")
  )

# grand pooled effect
grand_est <- tibble(
  brain_model = "Pooled",
  epi_model   = "Pooled",
  b           = model_assoc_mv_global$b[1,1],
  ci.lb       = model_assoc_mv_global$ci.lb[1],
  ci.ub       = model_assoc_mv_global$ci.ub[1]
)

# 4. Combine all estimates

plot_est <- bind_rows(
  pair_est %>% select(brain_model, epi_model, b, ci.lb, ci.ub),
  brain_margin_est,
  epi_margin_est,
  grand_est
)

# 5. Brain model order

brain_levels <- c(
  "Centile2",
  "DBN",
  "DevBrainAge",
  "ENIGMA",
  "Kaufmann",
  "PyBrainAge",
  "Pyment",
  "DunedinPACNI",
  "Pooled"
)

# 6. Apply factor levels carefully

pair_raw <- pair_raw %>%
  mutate(
    brain_model = factor(as.character(brain_model), levels = brain_levels),
    epi_model   = factor(as.character(epi_model), levels = epi_levels)
  )

plot_est <- plot_est %>%
  mutate(
    brain_model = factor(as.character(brain_model), levels = brain_levels),
    epi_model   = factor(as.character(epi_model), levels = epi_levels)
  )

plot_est <- plot_est %>%
  mutate(
    brain_model_facet = as.character(brain_model),
    brain_model_facet = ifelse(brain_model_facet == "Pooled",
                               "Brain age models (pooled)",
                               brain_model_facet),
    brain_model_facet = factor(
      brain_model_facet,
      levels = c(
        "Centile2",
        "DBN",
        "DevBrainAge",
        "ENIGMA",
        "Kaufmann",
        "PyBrainAge",
        "Pyment",
        "DunedinPACNI",
        "Brain age models (pooled)"
      )
    )
  )

pair_raw <- pair_raw %>%
  mutate(
    brain_model_facet = as.character(brain_model),
    brain_model_facet = ifelse(brain_model_facet == "Pooled",
                               "Brain age models (pooled)",
                               brain_model_facet),
    brain_model_facet = factor(
      brain_model_facet,
      levels = c(
        "Centile2",
        "DBN",
        "DevBrainAge",
        "ENIGMA",
        "Kaufmann",
        "PyBrainAge",
        "Pyment",
        "DunedinPACNI",
        "Brain age models (pooled)"
      )
    )
  )

# 7. Consistent cohort colours

# cohort_palette <- setNames(
#   viridis(length(cohort_levels)),
#   cohort_levels
# )

cohort_levels <- sort(unique(pair_raw$cohort))

# cohort_palette <- setNames(
#   kelly.colors(22)[-1], # removing the whitish colour that's not very visible 
#   cohort_levels
# )


cohort_palette <- c(
  "ALSPAC" = "#222222",
  "BHRC" = "#f3c300",
  "CannTeen" = "#875692",
  "DCHS" = "#f38400",
  "FinnBrain" = "#a1caf1",
  "GenR" = "#be0032",
  "GUSTO" = "#c2b280",
  "K2H Childhood" = "#848482",
  "K2H Infancy" = "#008856",
  "MTwiNS" = "#e68fac",
  "NICAP" = "#0067a5",
  "Oregon ADHD-1000" = "#f99379",
  "FFCWS" = "#604e97",
  "TAG" = "#f6a600",
  "UCI Echo" = "#b3446c"
)

# dividers 
pooled_divider_position <- 1.5
gen_divider_position <- 1.5 + length(gen2plus_models)

# 8. Plot

png(
  file.path(plots_dir, "plot_assoc_pairwise_with_margins_fixed_5cols_wide.png"),
  # width = 3600,
  # height = 2800,
  # res = 300
  width = 5616, #2340 #4680
  height = 3042, #1950 #2535
  res = 300
)

ggplot() +
  
  # reference line
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey85") +
  
  geom_hline(
    yintercept = pooled_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = gen_divider_position,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  # raw cohort-level associations
  geom_point(
    data = pair_raw,
    aes(
      x = RLM_Estimate_scaled,
      y = epi_model,
      colour = cohort
    ),
    alpha = .45,
    size = 2,
    position = position_jitter(height = .12)
  ) +
  
  # pooled meta-analysis CIs
  geom_segment(
    data = plot_est,
    aes(
      x = ci.lb,
      xend = ci.ub,
      y = epi_model,
      yend = epi_model,
      colour = "Overall pooled"
    ),
    linewidth = 0.9
  ) +
  
# pooled meta-analysis point estimates
  geom_point(
    data = plot_est,
    aes(
      x = b,
      y = epi_model,
      colour = "Overall pooled"
    ),
    shape = 18,
    size = 4
  ) +
  
  # one facet per brain model, including pooled column
  facet_wrap(~ brain_model_facet, ncol = 5, drop = FALSE) +
  
  scale_y_discrete(limits = epi_levels,
                   labels = function(x) ifelse(x == "Pooled", "Pooled (epi clocks)", x),
                   drop = FALSE) +
  
  scale_colour_manual(
    values = c(
      cohort_palette,
      "Overall pooled" = "black"
    ),
    breaks = c(names(cohort_palette), "Overall pooled")
  ) +
  
  labs(
    # title = "brain age - epi age associations",
    x = "Robust standardized beta",
    y = NULL,
    colour = "Source",
    title = "Brain-PAR - Epigenetic-PAR associations across model combinations",
  ) +
  
  coord_cartesian(clip = "off") +
  
  # theme(
  #   panel.grid.major.y = element_blank(),
  #   strip.text = element_text(face = "bold"),
  #   legend.position = "right"
  # )
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, hjust = 0.5),
    
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold"),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  ) +
  guides(
    colour = guide_legend(ncol = 1)
  )

dev.off()

##### USED association modelwise - legend inside (4A) ------

# load results needed 
model_assoc_mv_global<- readRDS(
  file.path(
    results_dir,
    "assoc_model2_PAR_global_final/model_assoc_mv_global.rds" #EW: updated input without GenR dupl: assoc_model2_PAR_global_final; previoiusly: assoc_model2_PAR_global_final
  )
)

model_assoc_mv_modelwise <- readRDS(
  file.path(
    results_dir,
    "assoc_model2_PAR_modelwise_final/model_assoc_mv_modelwise.rds" #EW: updated input without GenR dupl: assoc_model2_PAR_modelwise_final; previously: assoc_model2_PAR_modelwise_final
  )
)



# 1. Pairwise pooled estimates

pair_est <- tibble(
  term  = rownames(model_assoc_mv_modelwise$b),
  b     = as.numeric(model_assoc_mv_modelwise$b),
  se    = model_assoc_mv_modelwise$se,
  ci.lb = model_assoc_mv_modelwise$ci.lb,
  ci.ub = model_assoc_mv_modelwise$ci.ub,
  pval  = model_assoc_mv_modelwise$pval
) %>%
  mutate(term = str_remove(term, "^model_combi")) %>%
  separate(term, into = c("brain_model", "epi_model"), sep = "\\.")

pair_raw <- mod2_eff2 %>%
  select(cohort, brain_model, epi_model, RLM_Estimate_scaled)

# 2. Define epi clock groups and desired order

gen1_models <- sort(c(
  "AltumAge", "CorticalClock", "Hannum", "Horvath2013",
  "PCBrainAge", "PedBE", "Wu", "ZhangBLUP", "ZhangEN",
  "cAge", "skinHorvath"
))

gen2plus_models <- sort(c(
  "AdaptAge", "DamAge", "DNAmTL", "DunedinPACE",
  "PCGrimAge", "PhenoAge"
))

# We want TOP -> BOTTOM:
# Gen1 alphabetical, then Gen2plus alphabetical, then Pooled
# For ggplot y-axis, levels must be bottom -> top:
epi_levels <- c(
  "Pooled",
  rev(gen2plus_models),
  rev(gen1_models)
)

# 3. Marginal pooled effects

# pooled across epi clocks for each brain model
mod_brain_margin <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ brain_model - 1,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | epi_model
  ),
  data = mod2_eff2,
  method = "REML"
)

brain_margin_est <- tibble(
  brain_model = rownames(mod_brain_margin$b),
  epi_model   = "Pooled",
  b           = as.numeric(mod_brain_margin$b),
  ci.lb       = mod_brain_margin$ci.lb,
  ci.ub       = mod_brain_margin$ci.ub
) %>%
  mutate(
    brain_model = str_remove(brain_model, "^brain_model")
  )

# pooled across brain models for each epi model
mod_epi_margin <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ epi_model - 1,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | brain_model
  ),
  data = mod2_eff2,
  method = "REML"
)

epi_margin_est <- tibble(
  epi_model   = rownames(mod_epi_margin$b),
  brain_model = "Pooled",
  b           = as.numeric(mod_epi_margin$b),
  ci.lb       = mod_epi_margin$ci.lb,
  ci.ub       = mod_epi_margin$ci.ub
) %>%
  mutate(
    epi_model = str_remove(epi_model, "^epi_model")
  )

# grand pooled effect
grand_est <- tibble(
  brain_model = "Pooled",
  epi_model   = "Pooled",
  b           = model_assoc_mv_global$b[1,1],
  ci.lb       = model_assoc_mv_global$ci.lb[1],
  ci.ub       = model_assoc_mv_global$ci.ub[1]
)

# 4. Combine all estimates

plot_est <- bind_rows(
  pair_est %>% select(brain_model, epi_model, b, ci.lb, ci.ub),
  brain_margin_est,
  epi_margin_est,
  grand_est
)

# 5. Brain model order

brain_levels <- c(
  "Centile2",
  "DBN",
  "DevBrainAge",
  "ENIGMA",
  "Kaufmann",
  "PyBrainAge",
  "Pyment",
  "DunedinPACNI",
  "Pooled"
)

# 6. Apply factor levels carefully

pair_raw <- pair_raw %>%
  mutate(
    brain_model = factor(as.character(brain_model), levels = brain_levels),
    epi_model   = factor(as.character(epi_model), levels = epi_levels)
  )

plot_est <- plot_est %>%
  mutate(
    brain_model = factor(as.character(brain_model), levels = brain_levels),
    epi_model   = factor(as.character(epi_model), levels = epi_levels)
  )

plot_est <- plot_est %>%
  mutate(
    brain_model_facet = as.character(brain_model),
    brain_model_facet = ifelse(brain_model_facet == "Pooled",
                               "Brain age models (pooled)",
                               brain_model_facet),
    brain_model_facet = factor(
      brain_model_facet,
      levels = c(
        "Centile2",
        "DBN",
        "DevBrainAge",
        "ENIGMA",
        "Kaufmann",
        "PyBrainAge",
        "Pyment",
        "DunedinPACNI",
        "Brain age models (pooled)"
      )
    )
  )

pair_raw <- pair_raw %>%
  mutate(
    brain_model_facet = as.character(brain_model),
    brain_model_facet = ifelse(brain_model_facet == "Pooled",
                               "Brain age models (pooled)",
                               brain_model_facet),
    brain_model_facet = factor(
      brain_model_facet,
      levels = c(
        "Centile2",
        "DBN",
        "DevBrainAge",
        "ENIGMA",
        "Kaufmann",
        "PyBrainAge",
        "Pyment",
        "DunedinPACNI",
        "Brain age models (pooled)"
      )
    )
  )

# 7. Consistent cohort colours

# cohort_palette <- setNames(
#   viridis(length(cohort_levels)),
#   cohort_levels
# )

cohort_levels <- sort(unique(pair_raw$cohort))

# cohort_palette <- setNames(
#   kelly.colors(22)[-1], # removing the whitish colour that's not very visible 
#   cohort_levels
# )


cohort_palette <- c(
  "ALSPAC" = "#222222",
  "BHRC" = "#f3c300",
  "CannTeen" = "#875692",
  "DCHS" = "#f38400",
  "FinnBrain" = "#a1caf1",
  "GenR" = "#be0032",
  "GUSTO" = "#c2b280",
  "K2H Childhood" = "#848482",
  "K2H Infancy" = "#008856",
  "MTwiNS" = "#e68fac",
  "NICAP" = "#0067a5",
  "Oregon ADHD-1000" = "#f99379",
  "FFCWS" = "#604e97",
  "TAG" = "#f6a600",
  "UCI Echo" = "#b3446c"
)

# dividers 
pooled_divider_position <- 1.5
gen_divider_position <- 1.5 + length(gen2plus_models)

# 8. Plot

png(
  file.path(plots_dir, "plot_assoc_pairwise_with_margins_fixed_5cols_wide.png"),
  # width = 3600,
  # height = 2800,
  # res = 300
  width = 5616, #2340 #4680
  height = 3042, #1950 #2535
  res = 300
)

p <- ggplot() +
  
  # reference line
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey85") +
  
  geom_hline(
    yintercept = pooled_divider_position,
    linetype = "solid",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = gen_divider_position,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  
  # raw cohort-level associations
  geom_point(
    data = pair_raw,
    aes(
      x = RLM_Estimate_scaled,
      y = epi_model,
      colour = cohort
    ),
    alpha = .45,
    size = 2,
    position = position_jitter(height = .12)
  ) +
  
  # pooled meta-analysis CIs
  geom_segment(
    data = plot_est,
    aes(
      x = ci.lb,
      xend = ci.ub,
      y = epi_model,
      yend = epi_model,
      colour = "Overall pooled"
    ),
    linewidth = 0.9
  ) +
  
  # pooled meta-analysis point estimates
  geom_point(
    data = plot_est,
    aes(
      x = b,
      y = epi_model,
      colour = "Overall pooled"
    ),
    shape = 18,
    size = 4
  ) +
  
  # one facet per brain model, including pooled column
  facet_wrap(~ brain_model_facet, ncol = 5, drop = FALSE) +
  
  scale_y_discrete(limits = epi_levels,
                   labels = function(x) ifelse(x == "Pooled", "Pooled (epi clocks)", x),
                   drop = FALSE) +
  
  scale_colour_manual(
    values = c(
      cohort_palette,
      "Overall pooled" = "black"
    ),
    breaks = c(names(cohort_palette), "Overall pooled")
  ) +
  
  labs(
    # title = "brain age - epi age associations",
    x = "Robust standardized beta",
    y = NULL,
    colour = "Source",
    title = "Brain-PAR - Epigenetic-PAR associations across model combinations",
  ) +
  
  coord_cartesian(clip = "off") +
  
  # theme(
  #   panel.grid.major.y = element_blank(),
  #   strip.text = element_text(face = "bold"),
  #   legend.position = "right"
  # )
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15, hjust = 0.5),
    
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold"),
    
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  ) +
  guides(
    colour = guide_legend(ncol = 1)
  )

legend_plot <- cowplot::get_legend(p)

p_no_legend <- p + theme(legend.position = "none")

final_plot <- p_no_legend + patchwork::inset_element(
  legend_plot,
  left = 0.81,
  bottom = 0.1,
  right = 0.99,
  top = 0.47,
  align_to = "full"
)

print(final_plot)

dev.off()



