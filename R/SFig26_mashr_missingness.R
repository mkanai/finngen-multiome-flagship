library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

gex.data.strong <- readRDS("data/integrated_gex_batch1_5.mashr.data.strong.rds")
gex.is_missing  <- (gex.data.strong$Bhat == 0) & (gex.data.strong$Shat == 1000)
gex.per_ct_missing <- tibble::tibble(
  QTL = "eQTL",
  cell_type = colnames(gex.is_missing),
  per_ct_missing = colMeans(gex.is_missing)
)

atac.data.strong <- readRDS("data/integrated_atac_batch1_5.mashr.data.strong.rds")
atac.is_missing  <- (atac.data.strong$Bhat == 0) & (atac.data.strong$Shat == 1000)
atac.per_ct_missing <- tibble::tibble(
  QTL = "caQTL",
  cell_type = colnames(atac.is_missing),
  per_ct_missing = colMeans(atac.is_missing)
)

df.freq = read.table("tables/ST2_cell_type_freqs.tsv", TRUE, sep = "\t") %>%
  dplyr::mutate(QTL = ifelse(Atlas == "snRNA-Seq", "eQTL", "caQTL")) %>%
  dplyr::left_join(dplyr::bind_rows(gex.per_ct_missing, atac.per_ct_missing), by = c("QTL", "cell_type")) %>%
  tidyr::drop_na(per_ct_missing) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type),
                cell_type = factor(cell_type, levels = l1.order))

p.freq =
  ggplot(df.freq, aes(frac, per_ct_missing, color = cell_type, shape = QTL)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 1)) +
  scale_color_manual(values = l1.colors, labels = l1.labels) +
  scale_shape_manual(values = qtl.shapes) +
  locusviz::get_default_theme() +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Cell type frequency", y = "Per-cell-type missingness", color = "Cell type")

p.freq

cowplot::save_plot(
  "figures/SFig26_mashr_missingness.pdf",
  p.freq,
  base_height = 60,
  base_width = 60,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig26_mashr_missingness.png",
  p.freq,
  base_height = 60,
  base_width = 60,
  units = "mm",
  dpi = 300
)
