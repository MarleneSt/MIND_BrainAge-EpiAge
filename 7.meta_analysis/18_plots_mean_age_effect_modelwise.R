
library(tidyverse)
library(metafor)
library(patchwork)
library(scales)

results_dir <- "results"
plots_dir   <- "results/plots"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)


#### mean age effect modelwise ----

##### epi wMAE (2D) -----

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

epi_wMAE_mv_global_ageint <- readRDS(
  file.path(
    results_dir,
    "perf_epi_wMAE_noDunedinDNAmTL_global_ageint_final/model_epi_wMAE_mv_global_ageint_logscale.rds"
  )
)

epi_wMAE_mv_modelwise_ageint <- readRDS(
  file.path(
    results_dir,
    "perf_epi_wMAE_noDunedinDNAmTL_modelwise_ageint_final/model_epi_wMAE_mv_modelwise_ageint.rds"
  )
)
# overall 
age_model_est <- tibble(
  term = names(coef(epi_wMAE_mv_global_ageint)),
  estimate_log = as.numeric(coef(epi_wMAE_mv_global_ageint)),
  se_log = epi_wMAE_mv_global_ageint$se,
  ci.lb_log = epi_wMAE_mv_global_ageint$ci.lb,
  ci.ub_log = epi_wMAE_mv_global_ageint$ci.ub,
  pval = epi_wMAE_mv_global_ageint$pval
) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

# modelwise 
modelwise_age_est <- tibble(
  term = rownames(epi_wMAE_mv_modelwise_ageint$b),
  estimate_log = as.numeric(epi_wMAE_mv_modelwise_ageint$b),
  se_log = epi_wMAE_mv_modelwise_ageint$se,
  ci.lb_log = epi_wMAE_mv_modelwise_ageint$ci.lb,
  ci.ub_log = epi_wMAE_mv_modelwise_ageint$ci.ub,
  pval = epi_wMAE_mv_modelwise_ageint$pval
) %>%
  mutate(
    type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
    model = term %>%
      str_remove("^model") %>%
      str_remove(":mean_age_c") %>%
      str_remove("\\.mean_age_c")
  ) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

# --- plotting model-specific age slopes + global mean age effect ---

# --- model-specific age slopes + pooled/global age slope row ---

pooled_age_slope <- age_model_est %>%
  filter(term == "mean_age_c") %>%
  transmute(
    model = "Pooled",
    generation = "Pooled",
    type = "age_slope",
    percent_change = percent_change,
    percent_ci.lb = percent_ci.lb,
    percent_ci.ub = percent_ci.ub
  )

age_slopes <- modelwise_age_est %>%
  filter(type == "age_slope") %>%
  transmute(
    model = as.character(model),
    type,
    percent_change,
    percent_ci.lb,
    percent_ci.ub
  ) %>%
  mutate(
    generation = case_when(
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
        "DunedinPACE"
      ) ~ "Gen2-4",
      
      TRUE ~ "Other"
    )
  )

gen1_order <- age_slopes %>%
  filter(generation == "Gen1") %>%
  arrange(percent_change) %>%
  pull(model)

gen2_order <- age_slopes %>%
  filter(generation == "Gen2-4") %>%
  arrange(percent_change) %>%
  pull(model)

other_order <- age_slopes %>%
  filter(generation == "Other") %>%
  arrange(percent_change) %>%
  pull(model)

model_order <- c(
  "Pooled",
  rev(other_order),
  rev(gen2_order),
  rev(gen1_order)
)

plot_age_slopes <- bind_rows(
  age_slopes,
  pooled_age_slope
) %>%
  mutate(
    model = factor(model, levels = model_order)
  )

overall_divider_position <- 1.5
divider_position_other_gen2 <- 1.5 + length(other_order)
divider_position_gen2_gen1 <- 1.5 + length(other_order) + length(gen2_order)

x_range <- range(
  c(plot_age_slopes$percent_ci.lb, plot_age_slopes$percent_ci.ub),
  na.rm = TRUE
)

x_pad <- diff(x_range) * 0.08


png(
  filename = file.path(
    plots_dir,
    "plot_epi_modelwise_age_slopes_plus_pooled_by_generation.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

p <- ggplot(plot_age_slopes, aes(x = percent_change, y = model)) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.5
  ) +
  
  geom_segment(
    data = filter(plot_age_slopes, model != "Pooled"),
    aes(
      x = percent_ci.lb,
      xend = percent_ci.ub,
      y = model,
      yend = model,
      colour = model
    ),
    linewidth = 0.9,
    alpha = 0.6
  ) +
  
  geom_point(
    data = filter(plot_age_slopes, model != "Pooled"),
    aes(colour = model),
    shape = 18,
    size = 4,
    alpha = 1
  ) +
  
  geom_segment(
    data = filter(plot_age_slopes, model == "Pooled"),
    aes(
      x = percent_ci.lb,
      xend = percent_ci.ub,
      y = model,
      yend = model
    ),
    colour = "black",
    linewidth = 0.9,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  
    geom_point(
      data = filter(plot_age_slopes, model == "Pooled"),
      aes(x = percent_change, y = model),
      colour = "black",
      shape = 18,
      size = 5,
      alpha = 1,
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
  
  scale_y_discrete(
    limits = model_order
  ) +
  
  scale_x_continuous(
    limits = c(x_range[1] - x_pad, x_range[2] + x_pad)
  ) +
  
  scale_colour_manual(
    values = c(model_palette, "Pooled" = "black"),
    guide = "none"
  ) +
  
  labs(
    title = "Model-specific age moderation",
    caption = "negative = improved performance with increasing age\n positive = worse performance with increasing age",
    x = "% change in geometric mean wMAE per 1-year higher mean age",
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey30"),
    axis.title.x = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(
      size = 9,
      colour = "grey30",
      hjust = 0.5,
      margin = margin(t = 12)
    ),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  )

print(p)

dev.off()

##### brain wMAE (3D) -----

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

brain_wMAE_mv_global_ageint <- readRDS(
  file.path(
    results_dir,
    "perf_brain_wMAE_noDunedin_global_ageint_final/model_brain_wMAE_mv_global_ageint.rds"
  )
)

brain_wMAE_mv_modelwise_ageint <- readRDS(
  file.path(
    results_dir,
    "perf_brain_wMAE_noDunedin_modelwise_ageint_final/model_brain_wMAE_mv_modelwise_ageint.rds"
  )
)
# overall 
age_model_est <- tibble(
  term = names(coef(brain_wMAE_mv_global_ageint)),
  estimate_log = as.numeric(coef(brain_wMAE_mv_global_ageint)),
  se_log = brain_wMAE_mv_global_ageint$se,
  ci.lb_log = brain_wMAE_mv_global_ageint$ci.lb,
  ci.ub_log = brain_wMAE_mv_global_ageint$ci.ub,
  pval = brain_wMAE_mv_global_ageint$pval
) %>%
  mutate(
    estimate_raw = estimate_log,
    ci.lb_raw = ci.lb_log,
    ci.ub_raw = ci.ub_log
  )

# modelwise 
modelwise_age_est <- tibble(
  term = rownames(brain_wMAE_mv_modelwise_ageint$b),
  estimate_log = as.numeric(brain_wMAE_mv_modelwise_ageint$b),
  se_log = brain_wMAE_mv_modelwise_ageint$se,
  ci.lb_log = brain_wMAE_mv_modelwise_ageint$ci.lb,
  ci.ub_log = brain_wMAE_mv_modelwise_ageint$ci.ub,
  pval = brain_wMAE_mv_modelwise_ageint$pval
) %>%
  mutate(
    type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
    model = term %>%
      str_remove("^model") %>%
      str_remove(":mean_age_c") %>%
      str_remove("\\.mean_age_c")
  ) %>%
  mutate(
    estimate_raw = estimate_log,
    ci.lb_raw = ci.lb_log,
    ci.ub_raw = ci.ub_log
  )

# --- plotting model-specific age slopes + global mean age effect ---

# --- model-specific age slopes + pooled/global age slope row ---

pooled_age_slope <- age_model_est %>%
  filter(term == "mean_age_c") %>%
  transmute(
    model = "Pooled",
    generation = "Pooled",
    type = "age_slope",
    estimate_raw,
    ci.lb_raw,
    ci.ub_raw
  )

age_slopes <- modelwise_age_est %>%
  filter(type == "age_slope") %>%
  transmute(
    model = as.character(model),
    type,
    estimate_raw,
    ci.lb_raw,
    ci.ub_raw
  ) %>%
  mutate(
    generation = case_when(
      model %in% c(
        "Pyment", "DevBrainAge", "Centile2", "DBN",
        "Kaufmann", "PyBrainAge", "ENIGMA"
      ) ~ "Gen1",
      model %in% c("DunedinPACNI") ~ "Gen2-4",
      TRUE ~ "Other"
    )
  )

gen1_order <- age_slopes %>%
  filter(generation == "Gen1") %>%
  arrange(estimate_raw) %>%
  pull(model)

gen2_order <- age_slopes %>%
  filter(generation == "Gen2-4") %>%
  arrange(estimate_raw) %>%
  pull(model)

other_order <- age_slopes %>%
  filter(generation == "Other") %>%
  arrange(estimate_raw) %>%
  pull(model)

model_order <- c(
  "Pooled",
  rev(other_order),
  rev(gen2_order),
  rev(gen1_order)
)

plot_age_slopes <- bind_rows(
  age_slopes,
  pooled_age_slope
) %>%
  mutate(
    model = factor(model, levels = model_order)
  )

overall_divider_position <- 1.5
divider_position_other_gen2 <- 1.5 + length(other_order)
divider_position_gen2_gen1 <- 1.5 + length(other_order) + length(gen2_order)

x_range <- range(
  c(plot_age_slopes$ci.lb_raw, plot_age_slopes$ci.ub_raw),
  na.rm = TRUE
)

x_pad <- diff(x_range) * 0.08


png(
  filename = file.path(
    plots_dir,
    "plot_brain_wMAE_modelwise_age_slopes_plus_pooled_by_generation.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

p <- ggplot(plot_age_slopes, aes(x = estimate_raw, y = model)) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.5
  ) +
  
  geom_segment(
    data = filter(plot_age_slopes, model != "Pooled"),
    aes(
      x = ci.lb_raw,
      xend = ci.ub_raw,
      y = model,
      yend = model,
      colour = model
    ),
    linewidth = 0.9,
    alpha = 0.6
  ) +
  
  geom_point(
    data = filter(plot_age_slopes, model != "Pooled"),
    aes(colour = model),
    shape = 18,
    size = 4,
    alpha = 1
  ) +
  
  geom_segment(
    data = filter(plot_age_slopes, model == "Pooled"),
    aes(
      x = ci.lb_raw,
      xend = ci.ub_raw,
      y = model,
      yend = model
    ),
    colour = "black",
    linewidth = 0.9,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = filter(plot_age_slopes, model == "Pooled"),
    aes(x = estimate_raw, y = model),
    colour = "black",
    shape = 18,
    size = 5,
    alpha = 1,
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
  
  scale_y_discrete(
    limits = model_order
  ) +
  
  scale_x_continuous(
    limits = c(x_range[1] - x_pad, x_range[2] + x_pad)
  ) +
  
  scale_colour_manual(
    values = c(model_palette, "Pooled" = "black"),
    guide = "none"
  ) +
  
  labs(
    title = "Model-specific age moderation",
    caption = "negative = improved performance with increasing age\n positive = worse performance with increasing age",
    x = "Change in wMAE per 1-year higher mean age",
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey30"),
    axis.title.x = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(
      size = 9,
      colour = "grey30",
      hjust = 0.5,
      margin = margin(t = 12)
    ),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  )

print(p)

dev.off()

##### brain Pearson age slope plot -----

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


brain_pearson_mv_global_ageint <- readRDS(
  file.path(
    results_dir,
    "perf_brain_pearson_global_ageint_final/model_brain_pearson_mv_global_ageint.rds"
  )
)

brain_pearson_mv_modelwise_ageint <- readRDS(
  file.path(
    results_dir,
    "perf_brain_pearson_modelwise_ageint_final/model_brain_pearson_mv_modelwise_ageint.rds"
  )
)

age_model_est <- tibble(
  term = names(coef(brain_pearson_mv_global_ageint)),
  estimate_z = as.numeric(coef(brain_pearson_mv_global_ageint)),
  se_z = brain_pearson_mv_global_ageint$se,
  ci.lb_z = brain_pearson_mv_global_ageint$ci.lb,
  ci.ub_z = brain_pearson_mv_global_ageint$ci.ub,
  pval = brain_pearson_mv_global_ageint$pval
) %>%
  mutate(
    estimate_r = if_else(term == "intrcpt", transf.ztor(estimate_z), NA_real_),
    ci.lb_r = if_else(term == "intrcpt", transf.ztor(ci.lb_z), NA_real_),
    ci.ub_r = if_else(term == "intrcpt", transf.ztor(ci.ub_z), NA_real_)
  )

modelwise_age_est <- tibble(
  term = rownames(brain_pearson_mv_modelwise_ageint$b),
  estimate_z = as.numeric(brain_pearson_mv_modelwise_ageint$b),
  se_z = brain_pearson_mv_modelwise_ageint$se,
  ci.lb_z = brain_pearson_mv_modelwise_ageint$ci.lb,
  ci.ub_z = brain_pearson_mv_modelwise_ageint$ci.ub,
  pval = brain_pearson_mv_modelwise_ageint$pval
) %>%
  mutate(
    type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
    model = term %>%
      str_remove("^model") %>%
      str_remove(":mean_age_c") %>%
      str_remove("\\.mean_age_c"),
    
    # intercepts can be back-transformed to Pearson r
    estimate_r = if_else(type == "intercept", transf.ztor(estimate_z), NA_real_),
    ci.lb_r = if_else(type == "intercept", transf.ztor(ci.lb_z), NA_real_),
    ci.ub_r = if_else(type == "intercept", transf.ztor(ci.ub_z), NA_real_)
  )



# pooled/global age slope
pooled_age_slope <- age_model_est %>%
  filter(term == "mean_age_c") %>%
  transmute(
    model = "Pooled",
    generation = "Pooled",
    type = "age_slope",
    slope_z = estimate_z,
    ci.lb_z = ci.lb_z,
    ci.ub_z = ci.ub_z
  )

# model-specific age slopes
age_slopes <- modelwise_age_est %>%
  filter(type == "age_slope") %>%
  transmute(
    model = as.character(model),
    type,
    slope_z = estimate_z,
    ci.lb_z,
    ci.ub_z
  ) %>%
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
      ) ~ "Gen1",
      
      model %in% c("DunedinPACNI") ~ "Gen2-4",
      
      TRUE ~ "Other"
    )
  )

gen1_order <- age_slopes %>%
  filter(generation == "Gen1") %>%
  arrange(slope_z) %>%
  pull(model)

gen2_order <- age_slopes %>%
  filter(generation == "Gen2-4") %>%
  arrange(slope_z) %>%
  pull(model)

other_order <- age_slopes %>%
  filter(generation == "Other") %>%
  arrange(slope_z) %>%
  pull(model)

model_order <- c(
  "Pooled",
  other_order,
  gen2_order,
  gen1_order
)

plot_age_slopes <- bind_rows(
  age_slopes,
  pooled_age_slope
) %>%
  mutate(
    model = factor(model, levels = model_order)
  )

overall_divider_position <- 1.5
divider_position_other_gen2 <- 1.5 + length(other_order)
divider_position_gen2_gen1 <- 1.5 + length(other_order) + length(gen2_order)

x_range <- range(
  c(plot_age_slopes$ci.lb_z, plot_age_slopes$ci.ub_z),
  na.rm = TRUE
)

x_pad <- diff(x_range) * 0.08

png(
  filename = file.path(
    plots_dir,
    "plot_brain_pearson_modelwise_age_slopes_plus_pooled_by_generation.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

p <- ggplot(plot_age_slopes, aes(x = slope_z, y = model)) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.5
  ) +
  
  geom_segment(
    data = filter(plot_age_slopes, model != "Pooled"),
    aes(
      x = ci.lb_z,
      xend = ci.ub_z,
      y = model,
      yend = model,
      colour = model
    ),
    linewidth = 0.9,
    alpha = 0.6
  ) +
  
  geom_point(
    data = filter(plot_age_slopes, model != "Pooled"),
    aes(colour = model),
    shape = 18,
    size = 4
  ) +
  
  geom_segment(
    data = filter(plot_age_slopes, model == "Pooled"),
    aes(
      x = ci.lb_z,
      xend = ci.ub_z,
      y = model,
      yend = model
    ),
    colour = "black",
    linewidth = 0.9,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = filter(plot_age_slopes, model == "Pooled"),
    aes(x = slope_z, y = model),
    colour = "black",
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
  
  scale_y_discrete(limits = model_order) +
  
  scale_x_continuous(
    limits = c(x_range[1] - x_pad, x_range[2] + x_pad)
  ) +
  
  scale_colour_manual(
    values = c(model_palette, "Pooled" = "black"),
    guide = "none"
  ) +
  
  labs(
    title = "Model-specific age moderation",
    caption = "negative = weaker age–brain-age correlation with increasing age\npositive = stronger age–brain-age correlation with increasing age",
    x = "Change in Fisher-z Pearson correlation per 1-year higher mean age",
    y = NULL
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(
      size = 9,
      colour = "grey30",
      hjust = 0.5,
      margin = margin(t = 12)
    ),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  ) +
  
  coord_cartesian(clip = "off")

print(p)

dev.off()

##### epi Pearson age slope plot -----

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

epi_pearson_mv_global_ageint <- readRDS(
  file.path(
    results_dir,
    "perf_epi_pearson_global_ageint_final/model_epi_pearson_mv_global_ageint.rds"
  )
)

epi_pearson_mv_modelwise_ageint <- readRDS(
  file.path(
    results_dir,
    "perf_epi_pearson_modelwise_ageint_final/model_epi_pearson_mv_modelwise_ageint.rds"
  )
)

age_model_est <- tibble(
  term = names(coef(epi_pearson_mv_global_ageint)),
  estimate_z = as.numeric(coef(epi_pearson_mv_global_ageint)),
  se_z = epi_pearson_mv_global_ageint$se,
  ci.lb_z = epi_pearson_mv_global_ageint$ci.lb,
  ci.ub_z = epi_pearson_mv_global_ageint$ci.ub,
  pval = epi_pearson_mv_global_ageint$pval
) %>%
  mutate(
    estimate_r = if_else(term == "intrcpt", transf.ztor(estimate_z), NA_real_),
    ci.lb_r = if_else(term == "intrcpt", transf.ztor(ci.lb_z), NA_real_),
    ci.ub_r = if_else(term == "intrcpt", transf.ztor(ci.ub_z), NA_real_)
  )

modelwise_age_est <- tibble(
  term = rownames(epi_pearson_mv_modelwise_ageint$b),
  estimate_z = as.numeric(epi_pearson_mv_modelwise_ageint$b),
  se_z = epi_pearson_mv_modelwise_ageint$se,
  ci.lb_z = epi_pearson_mv_modelwise_ageint$ci.lb,
  ci.ub_z = epi_pearson_mv_modelwise_ageint$ci.ub,
  pval = epi_pearson_mv_modelwise_ageint$pval
) %>%
  mutate(
    type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
    model = term %>%
      str_remove("^model") %>%
      str_remove(":mean_age_c") %>%
      str_remove("\\.mean_age_c"),
    
    estimate_r = if_else(type == "intercept", transf.ztor(estimate_z), NA_real_),
    ci.lb_r = if_else(type == "intercept", transf.ztor(ci.lb_z), NA_real_),
    ci.ub_r = if_else(type == "intercept", transf.ztor(ci.ub_z), NA_real_)
  )

pooled_age_slope <- age_model_est %>%
  filter(term == "mean_age_c") %>%
  transmute(
    model = "Pooled",
    generation = "Pooled",
    type = "age_slope",
    slope_z = estimate_z,
    ci.lb_z = ci.lb_z,
    ci.ub_z = ci.ub_z
  )

age_slopes <- modelwise_age_est %>%
  filter(type == "age_slope") %>%
  transmute(
    model = as.character(model),
    type,
    slope_z = estimate_z,
    ci.lb_z,
    ci.ub_z
  ) %>%
  mutate(
    generation = case_when(
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
      
      TRUE ~ "Other"
    )
  )

gen1_order <- age_slopes %>%
  filter(generation == "Gen1") %>%
  arrange(slope_z) %>%
  pull(model)

gen2_order <- age_slopes %>%
  filter(generation == "Gen2-4") %>%
  arrange(slope_z) %>%
  pull(model)

other_order <- age_slopes %>%
  filter(generation == "Other") %>%
  arrange(slope_z) %>%
  pull(model)

model_order <- c(
  "Pooled",
  other_order,
  gen2_order,
  gen1_order
)

plot_age_slopes <- bind_rows(
  age_slopes,
  pooled_age_slope
) %>%
  mutate(
    model = factor(model, levels = model_order)
  )

overall_divider_position <- 1.5
divider_position_other_gen2 <- 1.5 + length(other_order)
divider_position_gen2_gen1 <- 1.5 + length(other_order) + length(gen2_order)

x_range <- range(
  c(plot_age_slopes$ci.lb_z, plot_age_slopes$ci.ub_z),
  na.rm = TRUE
)

x_pad <- diff(x_range) * 0.08

png(
  filename = file.path(
    plots_dir,
    "plot_epi_pearson_modelwise_age_slopes_plus_pooled_by_generation.png"
  ),
  width = 2340,
  height = 1950,
  res = 300
)

p <- ggplot(plot_age_slopes, aes(x = slope_z, y = model)) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.5
  ) +
  
  geom_segment(
    data = filter(plot_age_slopes, model != "Pooled"),
    aes(
      x = ci.lb_z,
      xend = ci.ub_z,
      y = model,
      yend = model,
      colour = model
    ),
    linewidth = 0.9,
    alpha = 0.6
  ) +
  
  geom_point(
    data = filter(plot_age_slopes, model != "Pooled"),
    aes(colour = model),
    shape = 18,
    size = 4
  ) +
  
  geom_segment(
    data = filter(plot_age_slopes, model == "Pooled"),
    aes(
      x = ci.lb_z,
      xend = ci.ub_z,
      y = model,
      yend = model
    ),
    colour = "black",
    linewidth = 0.9,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = filter(plot_age_slopes, model == "Pooled"),
    aes(x = slope_z, y = model),
    colour = "black",
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
  
  scale_y_discrete(limits = model_order) +
  
  scale_x_continuous(
    limits = c(x_range[1] - x_pad, x_range[2] + x_pad)
  ) +
  
  scale_colour_manual(
    values = c(model_palette, "Pooled" = "black"),
    guide = "none"
  ) +
  
  labs(
    title = "Model-specific age moderation",
    caption = "negative = weaker age–epigenetic-age correlation with increasing age\npositive = stronger age–epigenetic-age correlation with increasing age",
    x = "Change in Fisher-z Pearson correlation per 1-year higher mean age",
    y = NULL
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(
      size = 9,
      colour = "grey30",
      hjust = 0.5,
      margin = margin(t = 12)
    ),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20
    )
  ) +
  
  coord_cartesian(clip = "off")

print(p)

dev.off()

