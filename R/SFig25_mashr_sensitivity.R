library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

read_lfsr <- function(org_path, filt_path) {
  df.lfsr.org <- rgsutil::read_gsfile(org_path) %>%
    tidyr::pivot_longer(cols = -id,
                        names_to = "cell_type",
                        values_to = "lfsr")
  df.lfsr.filt <- rgsutil::read_gsfile(filt_path) %>%
    tidyr::pivot_longer(cols = -id,
                        names_to = "cell_type",
                        values_to = "lfsr")
  
  dplyr::inner_join(df.lfsr.org, df.lfsr.filt, by = c("id", "cell_type"))
}

df.lfsr.gex = read_lfsr(
  "gs://expansion_areas/multiome/batch1_5/mashr/integrated_gex_batch1_5.mashr.lfsr.tsv.gz",
  "gs://expansion_areas/multiome/batch1_5/mashr/sensitivity_0.5/integrated_gex_batch1_5.mashr.lfsr.tsv.gz"
) %>%
  dplyr::mutate(QTL = "eQTL")

df.lfsr.atac = read_lfsr(
  "gs://expansion_areas/multiome/batch1_5/mashr/integrated_atac_batch1_5.mashr.lfsr.tsv.gz",
  "gs://expansion_areas/multiome/batch1_5/mashr/sensitivity_0.5/integrated_atac_batch1_5.mashr.lfsr.tsv.gz"
) %>%
  dplyr::mutate(QTL = "caQTL")

df.lfsr.rate =
  dplyr::bind_rows(df.lfsr.gex, df.lfsr.atac) %>%
  dplyr::group_by(QTL, cell_type) %>%
  dplyr::summarize(locusviz::binom_ci(sum((lfsr.x < 0.05) == (lfsr.y < 0.05)), n()),
                   locusviz::spearman_ci(lfsr.x, lfsr.y),
                   .groups = "drop") %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))

with(df.lfsr.gex, locusviz::spearman_ci(lfsr.x, lfsr.y))
with(df.lfsr.atac, locusviz::spearman_ci(lfsr.x, lfsr.y))


plot_lfsr_density <- function(df,
                              title,
                              bw = c(0.01, 0.01),
                              hide.ylab = FALSE) {
  dens <- KernSmooth::bkde2D(cbind(df$lfsr.x, df$lfsr.y),
                             bandwidth = bw,
                             gridsize = c(500, 500))
  df$density <- fields::interp.surface(list(
    x = dens$x1,
    y = dens$x2,
    z = dens$fhat
  ),
  cbind(df$lfsr.x, df$lfsr.y))
  
  ggplot(df, aes(lfsr.x, lfsr.y, color = density)) +
    ggrastr::rasterize(geom_point(size = 1, stroke = 0), dpi = 300) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    scale_x_continuous(expand = expansion()) +
    scale_y_continuous(expand = expansion()) +
    scale_color_viridis_c(trans = "log10", option = "inferno", labels = scales::label_log()) +
    locusviz::get_default_theme(
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      hide.ylab = hide.ylab
    ) +
    theme(
      plot.title = element_text(hjust = 0.02),
      legend.direction = "horizontal",
      legend.key.width = unit(4, "mm"),
      legend.title = element_text(margin = margin(b = 1, unit = "mm")),
      legend.title.position = "top"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      x = "LFSR (original)",
      y = "LFSR (filtered)",
      color = "Density",
      title = title
    )
}

p.gex <- dplyr::filter(df.lfsr.gex, cell_type == "predicted.celltype.l1.PBMC") %>%
  plot_lfsr_density("eQTL (PBMC)")
p.gex

p.atac <- dplyr::filter(df.lfsr.atac, cell_type == "predicted.celltype.l1.PBMC") %>%
  plot_lfsr_density("caQTL (PBMC)", hide.ylab = TRUE)
p.atac

cell_type.levels =
  dplyr::group_by(df.lfsr.rate, cell_type) %>%
  dplyr::summarize(frac = mean(frac), .groups = "drop") %>%
  dplyr::arrange(frac) %>%
  dplyr::pull(cell_type)

p.sig.rate =
  dplyr::mutate(df.lfsr.rate, cell_type = factor(cell_type, cell_type.levels)) %>%
  ggplot(aes(frac, cell_type, color = QTL)) +
  geom_errorbar(aes(xmin = frac_lower, xmax = frac_upper),
                orientation = "y",
                width = 0) +
  geom_point() +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion(c(0, 0.05))) +
  scale_y_discrete(labels = l1.labels) +
  scale_color_manual(values = qtl.colors) +
  locusviz::get_default_theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    hide.ytitle = TRUE
  ) +
  labs(x = "% significance concordance")

p.rho =
  dplyr::mutate(df.lfsr.rate, cell_type = factor(cell_type, cell_type.levels)) %>%
  ggplot(aes(rho, cell_type, color = QTL)) +
  geom_errorbar(aes(xmin = rho_lower, xmax = rho_upper),
                orientation = "y",
                width = 0) +
  geom_point() +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(expand = expansion(c(0, 0.05))) +
  scale_color_manual(values = qtl.colors) +
  locusviz::get_default_theme(legend.position = "none", hide.ylab = TRUE) +
  labs(x = expression(paste("LFSR correlation ", italic(rho))))

plt = p.gex + p.atac + p.sig.rate + p.rho + patchwork::plot_layout(nrow = 2) + patchwork::plot_annotation(tag_levels = "a")
plt

cowplot::save_plot(
  "figures/SFig25_mashr_sensitivity.pdf",
  plt,
  base_height = 120,
  base_width = 120,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig25_mashr_sensitivity.png",
  plt,
  base_height = 120,
  base_width = 120,
  units = "mm",
  dpi = 300
)
