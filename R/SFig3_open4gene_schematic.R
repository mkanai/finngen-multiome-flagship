library(dplyr)
library(ggplot2)
source(here::here("R/const.R"))

set.seed(42)
n <- 500

# Switch: sharp logistic on/off, completely flat expression among expressors
sim_switch <- function(n) {
  acc <- rbeta(n, 2, 2)
  # Sharp logistic: nearly all off at low acc, nearly all on at high acc
  p_expr <- plogis(-8 + 16 * acc)
  expressed <- rbinom(n, 1, p_expr)
  # Among expressors, expression is constant (no trend), tighter variance
  expr <- ifelse(expressed, rlnorm(n, 3.4, 0.2), 0)
  data.frame(
    accessibility = acc,
    expression = expr,
    expressed = factor(expressed)
  )
}

# Rheostat: nearly all cells express, but level scales strongly with accessibility
sim_rheostat <- function(n) {
  acc <- rbeta(n, 2, 2)
  # Almost all cells express (high baseline, no accessibility effect)
  p_expr <- plogis(2.5)
  expressed <- rbinom(n, 1, p_expr)
  # Among expressors: high baseline (~25 at x=0) with moderate slope
  expr <- ifelse(expressed, rlnorm(n, 3.2 + 1.0 * acc, 0.25), 0)
  data.frame(
    accessibility = acc,
    expression = expr,
    expressed = factor(expressed)
  )
}

# Dual: sharp logistic on/off AND strong dose-response
sim_dual <- function(n) {
  acc <- rbeta(n, 2, 2)
  # Sharp logistic (same as switch)
  p_expr <- plogis(-8 + 16 * acc)
  expressed <- rbinom(n, 1, p_expr)
  # Among expressors, strong dose-response spanning similar range
  expr <- ifelse(expressed, rlnorm(n, 1.5 + 2.5 * acc, 0.3), 0)
  data.frame(
    accessibility = acc,
    expression = expr,
    expressed = factor(expressed)
  )
}

df_switch <- cbind(sim_switch(n), mode = "Switch")
df_rheo <- cbind(sim_rheostat(n), mode = "Rheostat")
df_dual <- cbind(sim_dual(n), mode = "Dual")

df <- rbind(df_switch, df_rheo, df_dual)
df$mode <- factor(df$mode, levels = c("Switch", "Rheostat", "Dual"))
df$acc_bin <- cut(df$accessibility,
                  breaks = seq(0, 1, 0.1),
                  include.lowest = TRUE)
prop_df <- aggregate(
  as.numeric(as.character(df$expressed)),
  by = list(acc_bin = df$acc_bin, mode = df$mode),
  FUN = mean
)
names(prop_df)[3] <- "prop_expressed"
prop_df$acc_mid <- seq(0.05, 0.95, 0.1)[as.numeric(prop_df$acc_bin)]

all_bins <- levels(df$acc_bin)
all_modes <- levels(df$mode)
complete_grid <- expand.grid(
  acc_bin = factor(all_bins, levels = all_bins),
  mode = factor(all_modes, levels = all_modes)
)

expr_df <- subset(df, expressed == 1)
expr_df$acc_bin <- cut(expr_df$accessibility,
                       breaks = seq(0, 1, 0.1),
                       include.lowest = TRUE)

mean_df_raw <- aggregate(
  expr_df$expression,
  by = list(acc_bin = expr_df$acc_bin, mode = expr_df$mode),
  FUN = mean
)
names(mean_df_raw)[3] <- "mean_expr"

mean_df <- merge(complete_grid,
                 mean_df_raw,
                 by = c("acc_bin", "mode"),
                 all.x = TRUE)

mean_df$acc_mid <- seq(0.05, 0.95, 0.1)[as.numeric(mean_df$acc_bin)]


plot_scatter = function(df,
                        mode,
                        legend.position = "none",
                        hide.ytitle = TRUE) {
  df = dplyr::filter(df, mode == .env$mode)
  
  ggplot(df, aes(x = accessibility, y = expression)) +
    geom_point(aes(color = expressed)) +
    geom_smooth(
      data = subset(df, expressed == 1),
      method = "lm",
      se = FALSE,
      color = "grey20",
      linetype = "dashed"
    ) +
    scale_color_manual(
      values = c("0" = "grey50", "1" = unname(link_mode_colors[mode])),
      labels = c("0" = "Not expressed", "1" = "Expressed"),
      name = NULL
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = "Peak accessibility", y = "Gene expression", title = mode) +
    locusviz::get_default_theme(
      legend.position = legend.position,
      legend.justification = legend.position,
      hide.ytitle = hide.ytitle
    )
}

plot_zero = function(prop_df, mode, hide.ytitle = TRUE) {
  prop_df = dplyr::filter(prop_df, mode == .env$mode)
  
  ggplot(prop_df, aes(x = acc_mid, y = prop_expressed)) +
    geom_point(color = link_mode_colors[mode]) +
    geom_smooth(method = "loess",
                se = FALSE,
                color = link_mode_colors[mode]) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Peak accessibility", y = "P(expressed)") +
    locusviz::get_default_theme(hide.ytitle = hide.ytitle)
}

plot_count = function(mean_df, mode, hide.ytitle = TRUE) {
  mean_df = dplyr::filter(mean_df, mode == .env$mode)
  max_na_acc_mid = dplyr::filter(mean_df, is.na(mean_expr)) %>%
    dplyr::pull(acc_mid) %>%
    max()
  
  ggplot(mean_df, aes(x = acc_mid, y = mean_expr)) +
    locusviz::or_missing(is.finite(max_na_acc_mid), geom_rect(
      aes(
        xmin = 0,
        xmax = max_na_acc_mid + 0.05,
        ymin = 0,
        ymax = Inf
      ),
      fill = "grey90"
    )) +
    geom_point(color = link_mode_colors[mode]) +
    geom_smooth(method = "loess",
                se = FALSE,
                color = link_mode_colors[mode]) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 65)) +
    labs(x = "Peak accessibility", y = "E[expr | expr > 0]") +
    locusviz::get_default_theme(hide.ytitle = hide.ytitle)
}

p.open4gene.schematic =
  list(
    plot_scatter(
      df,
      "Switch",
      legend.position = c(1, 1),
      hide.ytitle = FALSE
    ),
    plot_scatter(df, "Rheostat"),
    plot_scatter(df, "Dual"),
    plot_zero(prop_df, "Switch", hide.ytitle = FALSE),
    plot_zero(prop_df, "Rheostat"),
    plot_zero(prop_df, "Dual"),
    plot_count(mean_df, "Switch", hide.ytitle = FALSE),
    plot_count(mean_df, "Rheostat"),
    plot_count(mean_df, "Dual")
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig3_open4gene_schematic.pdf",
  p.open4gene.schematic,
  base_height = 170,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig3_open4gene_schematic.png",
  p.open4gene.schematic,
  base_height = 170,
  base_width = 180,
  units = "mm",
  dpi = 300
)

