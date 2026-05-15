library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.gex.batch = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/sensitivity/batch/integrated_gex_batch1_5.fgid.qc.predicted.celltype.l1.PBMC.summary.tsv.gz"
) %>%
  dplyr::mutate(QTL = "eQTL")
df.atac.batch = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/sensitivity/batch/integrated_atac_batch1_5.fgid.predicted.celltype.l1.PBMC.summary.tsv.gz"
)  %>%
  dplyr::mutate(QTL = "caQTL")

df.batch = dplyr::bind_rows(df.gex.batch, df.atac.batch) %>%
  dplyr::mutate(
    batch_col = factor(
      batch_col,
      levels = rev(c("pool_name", "flow_cell_id", "sequencer")),
      labels = c("Multiome channel", "Flow cell", "Sequencer")
    ),
    phase = factor(
      phase,
      levels = c("pre", "post"),
      labels = c("Pre-adjustment", "Post-adjustment")
    ),
    y = interaction(batch_col, QTL)
  ) %>%
  tidyr::drop_na()

p.batch =
  ggplot(df.batch, aes(adj_median, y, color = phase)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbar(aes(xmin = adj_q25, xmax = adj_q75), width = 0) +
  geom_point() +
  locusviz::get_default_theme(hide.ytitle = TRUE) +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = BuenColors::jdb_palette("corona")[c(3, 4)]) +
  scale_y_discrete(guide = legendry::guide_axis_nested(
    key = legendry::key_range_auto(sep = "\\."),
    levels_text = list(NULL, element_text(angle = 90, hjust = 0.5))
  )) +
  labs(x = expression(paste("Adjusted ", italic(R)^2)))

cowplot::save_plot(
  "figures/SFig23_batch_variance.pdf",
  p.batch,
  base_height = 90,
  base_width = 90,
  units = "mm"
)

cowplot::save_plot(
  "figures/SFig23_batch_variance.png",
  p.batch,
  base_height = 90,
  base_width = 90,
  units = "mm",
  dpi = 300
)

################################################################################

df.gex.lisi = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/sensitivity/batch/integrated_gex_batch1_5.fgid.qc.umap.lisi.seed42.per_cell.tsv.gz"
) %>%
  dplyr::mutate(Atlas = "snRNA")
df.atac.lisi = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/sensitivity/batch/integrated_atac_batch1_5.fgid.umap.lisi.seed42.per_cell.tsv.gz"
)  %>%
  dplyr::mutate(Atlas = "snATAC")

df.lisi = dplyr::bind_rows(df.gex.lisi, df.atac.lisi) %>%
  dplyr::select(Atlas,
                pool_name,
                flow_cell_id,
                sequencer,
                predicted.celltype.l1) %>%
  tidyr::pivot_longer(-Atlas, names_to = "batch_col", values_to = "lisi")


plot_lisi = function(df) {
  ggplot(df, aes(lisi, Atlas, fill = Atlas)) +
    ggridges::geom_density_ridges(alpha = 0.8, linewidth = 0) +
    locusviz::get_default_theme(hide.ytitle = TRUE, legend.position = "none") +
    scale_fill_manual(values = BuenColors::jdb_palette("corona")) +
    scale_x_continuous(expand = expansion()) +
    scale_y_discrete(expand = expansion(mult = 0, add = c(0, 2))) +
    coord_cartesian(xlim = c(0, ceiling(max(df$lisi))))
}

p.lisi = list(
  dplyr::filter(df.lisi, batch_col == "pool_name") %>%
    plot_lisi() +
    labs(x = "iLISI", title = "Multiome channel"),
  dplyr::filter(df.lisi, batch_col == "flow_cell_id") %>%
    plot_lisi() +
    labs(x = "iLISI", title = "Flow cell"),
  dplyr::filter(df.lisi, batch_col == "sequencer") %>%
    plot_lisi() +
    labs(x = "iLISI", title = "Sequencer"),
  dplyr::filter(df.lisi, batch_col == "predicted.celltype.l1") %>%
    plot_lisi() +
    labs(x = "cLISI", title = "L1 cell type")
) %>%
  purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig24_batch_lisi.pdf",
  p.lisi,
  base_height = 120,
  base_width = 120,
  units = "mm"
)

cowplot::save_plot(
  "figures/SFig24_batch_lisi.png",
  p.lisi,
  base_height = 120,
  base_width = 120,
  units = "mm",
  dpi = 300
)