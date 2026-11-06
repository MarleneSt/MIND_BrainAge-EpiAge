### plots for heatmaps 

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
library(corrplot)

#### define important folders ----
home <- "data" # [REPRODUCIBILITY] relative path
data_dir <- file.path(home)
results_dir <- "results/within_modality/" # [REPRODUCIBILITY] relative path
plots_dir_withinmodality <- "results/within_modality/plots" # [REPRODUCIBILITY] relative path
dir.create(plots_dir_withinmodality, recursive = TRUE, showWarnings = FALSE) # [REPRODUCIBILITY] create folder on a fresh run


# to think about - the unclustered plots could be ordered by generation


##### PAR brain -----
# read in file if needed 
pair_est <- read.csv(
  file.path(results_dir,
    "within_modality_corr_brain_PAR_final/pairwise_pooled_correlations.csv"
  )
)


#### BRAIN PAR ----

#pair_est <- res_PAR_brain$pair_est
# or load csv from results_dir
analysis_name <- "brain-PAR"

# make symmetric matrix of pooled r values
corr_long <- pair_est %>%
  select(model_a, model_b, pooled_r) %>%
  bind_rows(
    pair_est %>%
      select(model_a, model_b, pooled_r) %>%
      rename(model_a = model_b, model_b = model_a)
  )

models <- sort(unique(c(corr_long$model_a, corr_long$model_b)))

diag_df <- tibble(
  model_a = models,
  model_b = models,
  pooled_r = 1
)

corr_long <- bind_rows(corr_long, diag_df)

corrmat <- corr_long %>%
  pivot_wider(
    names_from = model_b,
    values_from = pooled_r
  ) %>%
  column_to_rownames("model_a") %>%
  as.matrix()

corrmat <- corrmat[models, models]

###### unclustered color ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_unclustered_colour_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "color",
  order = "original",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = TRUE,
  addgrid.col = "grey85",
  rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nunclustered\n", analysis_name)
)
dev.off()

###### clustered color ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_clustered_color_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)


corrplot(
  corrmat,
  type = "lower",
  method = "color",
  order = "hclust",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = TRUE,
  addgrid.col = "grey85",
  rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nclustered\n", analysis_name)
)

dev.off()

###### unclustered circle ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_unclustered_circle_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "circle",
  order = "original",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = FALSE,
  # addgrid.col = "grey85",
  # rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nunclustered\n", analysis_name)
)

dev.off()

###### clustered circle -----
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_clustered_circle_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "circle",
  order = "hclust",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = FALSE,
  # addgrid.col = "grey85",
  # rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nclustered\n", analysis_name)
)

dev.off()

#### BRAIN predage ----

#pair_est <- res_predage_brain$pair_est
# or load csv from results_dir
# read in file if needed 
pair_est <- read.csv(
  file.path(results_dir,
            "within_modality_corr_brain_predage_final/pairwise_pooled_correlations.csv"
  )
)

analysis_name <- "brain-predicted age"

# make symmetric matrix of pooled r values
corr_long <- pair_est %>%
  select(model_a, model_b, pooled_r) %>%
  bind_rows(
    pair_est %>%
      select(model_a, model_b, pooled_r) %>%
      rename(model_a = model_b, model_b = model_a)
  )

models <- sort(unique(c(corr_long$model_a, corr_long$model_b)))

diag_df <- tibble(
  model_a = models,
  model_b = models,
  pooled_r = 1
)

corr_long <- bind_rows(corr_long, diag_df)

corrmat <- corr_long %>%
  pivot_wider(
    names_from = model_b,
    values_from = pooled_r
  ) %>%
  column_to_rownames("model_a") %>%
  as.matrix()

corrmat <- corrmat[models, models]

###### unclustered color ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_unclustered_colour_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "color",
  order = "original",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = TRUE,
  addgrid.col = "grey85",
  rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nunclustered\n", analysis_name)
)
dev.off()

###### clustered color ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_clustered_color_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)


corrplot(
  corrmat,
  type = "lower",
  method = "color",
  order = "hclust",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = TRUE,
  addgrid.col = "grey85",
  rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nclustered\n", analysis_name)
)

dev.off()

###### unclustered circle ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_unclustered_circle_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "circle",
  order = "original",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = FALSE,
  # addgrid.col = "grey85",
  # rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nunclustered\n", analysis_name)
)

dev.off()

###### clustered circle -----
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_clustered_circle_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "circle",
  order = "hclust",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = FALSE,
  # addgrid.col = "grey85",
  # rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nclustered\n", analysis_name)
)


dev.off()

#### Epi PAR ----

#pair_est <- res_PAR_epi$pair_est
# or load csv from results_dir
pair_est <- read.csv(
  file.path(results_dir,
            "within_modality_corr_epi_PAR_final/pairwise_pooled_correlations.csv"
  ))

analysis_name <- "epi-PAR"

# make symmetric matrix of pooled r values
corr_long <- pair_est %>%
  select(model_a, model_b, pooled_r) %>%
  bind_rows(
    pair_est %>%
      select(model_a, model_b, pooled_r) %>%
      rename(model_a = model_b, model_b = model_a)
  )

models <- sort(unique(c(corr_long$model_a, corr_long$model_b)))

diag_df <- tibble(
  model_a = models,
  model_b = models,
  pooled_r = 1
)

corr_long <- bind_rows(corr_long, diag_df)

corrmat <- corr_long %>%
  pivot_wider(
    names_from = model_b,
    values_from = pooled_r
  ) %>%
  column_to_rownames("model_a") %>%
  as.matrix()

corrmat <- corrmat[models, models]

###### unclustered color ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_unclustered_colour_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "color",
  order = "original",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = TRUE,
  addgrid.col = "grey85",
  rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nunclustered\n", analysis_name)
)
dev.off()

###### clustered color ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_clustered_color_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)


corrplot(
  corrmat,
  type = "lower",
  method = "color",
  order = "hclust",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = TRUE,
  addgrid.col = "grey85",
  rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nclustered\n", analysis_name)
)

dev.off()

###### unclustered circle ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_unclustered_circle_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "circle",
  order = "original",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = FALSE,
  # addgrid.col = "grey85",
  # rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nunclustered\n", analysis_name)
)

dev.off()

###### clustered circle -----
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_clustered_circle_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "circle",
  order = "hclust",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = FALSE,
  # addgrid.col = "grey85",
  # rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nclustered\n", analysis_name)
)

dev.off()

#### Epi predage ----

#pair_est <- res_predage_epi$pair_est
# or load csv from results_dir
pair_est <- read.csv(
  file.path(results_dir,
            "within_modality_corr_epi_predage_final/pairwise_pooled_correlations.csv"
  ))


analysis_name <- "epi-predicted age"

# make symmetric matrix of pooled r values
corr_long <- pair_est %>%
  select(model_a, model_b, pooled_r) %>%
  bind_rows(
    pair_est %>%
      select(model_a, model_b, pooled_r) %>%
      rename(model_a = model_b, model_b = model_a)
  )

models <- sort(unique(c(corr_long$model_a, corr_long$model_b)))

diag_df <- tibble(
  model_a = models,
  model_b = models,
  pooled_r = 1
)

corr_long <- bind_rows(corr_long, diag_df)

corrmat <- corr_long %>%
  pivot_wider(
    names_from = model_b,
    values_from = pooled_r
  ) %>%
  column_to_rownames("model_a") %>%
  as.matrix()

corrmat <- corrmat[models, models]

###### unclustered color ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_unclustered_colour_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "color",
  order = "original",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = TRUE,
  addgrid.col = "grey85",
  rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nunclustered\n", analysis_name)
)
dev.off()

###### clustered color ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_clustered_color_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)


corrplot(
  corrmat,
  type = "lower",
  method = "color",
  order = "hclust",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = TRUE,
  addgrid.col = "grey85",
  rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nclustered\n", analysis_name)
)

dev.off()

###### unclustered circle ------
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_unclustered_circle_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "circle",
  order = "original",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = FALSE,
  # addgrid.col = "grey85",
  # rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nunclustered\n", analysis_name)
)

dev.off()

###### clustered circle -----
png(
  file.path(
    plots_dir_withinmodality,
    paste0("corrplot_modelcombination_r_clustered_circle_", analysis_name, ".png")
  ),
  width = 2340,
  height = 2340,
  res = 300
)

par(
  mar = c(0, 0, 3, 0),
  family = "sans",
  cex.main = 1,
  font.main = 1
)

corrplot(
  corrmat,
  type = "lower",
  method = "circle",
  order = "hclust",
  diag = FALSE,
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 50,
  tl.cex = 0.75,
  number.cex = 0.55,
  cl.cex = 0.75,
  outline = FALSE,
  # addgrid.col = "grey85",
  # rect.col = "grey35",
  rect.lwd = 1,
  mar = c(0, 0, 3, 0),
  title = paste0("Within-modality meta-analytic correlations (Pearson's r)\nclustered\n", analysis_name)
)

dev.off()






#### heatmap currently not used ----

# # symmetric long dataframe
# corr_plot_df <- pair_est %>%
#   select(model_a, model_b, pooled_r) %>%
#   bind_rows(
#     pair_est %>%
#       select(model_a, model_b, pooled_r) %>%
#       rename(
#         model_a = model_b,
#         model_b = model_a
#       )
#   )
# 
# models <- sort(unique(c(corr_plot_df$model_a, corr_plot_df$model_b)))
# 
# diag_df <- tibble(
#   model_a = models,
#   model_b = models,
#   pooled_r = 1
# )
# 
# corr_plot_df <- bind_rows(corr_plot_df, diag_df)
# 
# 
# # --- unclustered plot ---
# 
# corr_plot_df_unclustered <- corr_plot_df %>%
#   mutate(
#     model_a = factor(model_a, levels = models),
#     model_b = factor(model_b, levels = rev(models))
#   )
# 
# png(
#   filename = file.path(
#     plots_dir_withinmodality,
#     paste0("heatmap_modelcombination_r_unclustered_", analysis_name, ".png")
#   ),
#   width = 3000,
#   height = 2700,
#   res = 300
# )
# 
# ggplot(
#   corr_plot_df_unclustered,
#   aes(x = model_a, y = model_b, fill = pooled_r)
# ) +
#   geom_tile(color = "white", linewidth = 0.25) +
#   coord_equal() +
#   scale_fill_gradient2(
#     limits = c(-1, 1),
#     midpoint = 0,
#     name = "Pooled r"
#   ) +
#   labs(
#     x = NULL,
#     y = NULL,
#     title = paste0("Pooled model-combination correlations: ", analysis_name),
#     subtitle = "Unclustered model order"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid = element_blank()
#   )
# 
# dev.off()
# 
# 
# # --- clustered plot ---
# 
# corr_wide <- corr_plot_df %>%
#   pivot_wider(
#     names_from = model_b,
#     values_from = pooled_r
#   ) %>%
#   column_to_rownames("model_a") %>%
#   as.matrix()
# 
# corr_wide <- corr_wide[models, models]
# 
# hc <- hclust(
#   as.dist(1 - corr_wide),
#   method = "complete"
# )
# 
# clustered_models <- hc$labels[hc$order]
# 
# corr_plot_df_clustered <- corr_plot_df %>%
#   mutate(
#     model_a = factor(model_a, levels = clustered_models),
#     model_b = factor(model_b, levels = rev(clustered_models))
#   )
# 
# png(
#   filename = file.path(
#     plots_dir_withinmodality,
#     paste0("heatmap_modelcombination_r_clustered_", analysis_name, ".png")
#   ),
#   width = 3000,
#   height = 2700,
#   res = 300
# )
# 
# ggplot(
#   corr_plot_df_clustered,
#   aes(x = model_a, y = model_b, fill = pooled_r)
# ) +
#   geom_tile(color = "white", linewidth = 0.25) +
#   coord_equal() +
#   scale_fill_gradient2(
#     limits = c(-1, 1),
#     midpoint = 0,
#     name = "Pooled r"
#   ) +
#   labs(
#     x = NULL,
#     y = NULL,
#     title = paste0("Pooled model-combination correlations: ", analysis_name),
#     subtitle = "Hierarchically clustered model order"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid = element_blank()
#   )
# 
# dev.off()
# 
# write.csv(
#   tibble(model = clustered_models),
#   file.path(
#     plots_dir_withinmodality,
#     paste0("heatmap_clustered_order_", analysis_name, ".csv")
#   ),
#   row.names = FALSE
# )

