library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.eqtl.max_pip <- rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/saige_qtl/susie/full.max_pip.annot.most_severe.txt.bgz"
) %>%
  dplyr::transmute(variant_id = variant,
                   max_pip = max_pip,
                   QTL = "eQTL")

df.caqtl.max_pip <- rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/atac_results/susie/full.max_pip.annot.most_severe.txt.bgz"
) %>%
  dplyr::transmute(variant_id = variant,
                   max_pip = max_pip,
                   QTL = "caQTL")

df.microarray.meta = data.table::fread(
  "https://zenodo.org/record/7808390/files/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz",
  data.table = FALSE
) %>%
  dplyr::select(phenotype_id, gene_id)

df.stim <- rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/eQTL_Catalogue/finemapped_resources_finngen_version_20250513.stimulated.tsv.gz"
) %>%
  dplyr::filter(data_type == "eQTL" &
                  quant_method %in% c("ge", "microarray")) %>%
  dplyr::left_join(df.microarray.meta, by = c("trait" = "phenotype_id")) %>%
  dplyr::mutate(
    is_naive = stringr::str_starts(condition_label, "naive"),
    gene_id = ifelse(
      is.na(gene_id),
      stringr::str_extract(trait, "^ENSG[0-9]+"),
      gene_id
    )
  )

df.max_pip.stim <-
  dplyr::filter(df.stim, !is_naive) %>%
  dplyr::group_by(variant_id) %>%
  dplyr::summarize(max_pip_stim = max(pip), .groups = "drop")

df.max_pip.naive <-
  dplyr::filter(df.stim, is_naive) %>%
  dplyr::group_by(variant_id) %>%
  dplyr::summarize(max_pip_naive = max(pip), .groups = "drop")

df.max_pip.stim_naive = dplyr::full_join(df.max_pip.naive, df.max_pip.stim, by = "variant_id") %>%
  dplyr::mutate(
    category = dplyr::case_when(
      (max_pip_naive > 0.5) & (max_pip_stim > 0.5) ~ "Shared",
      (max_pip_stim > 0.5) &
        (is.na(max_pip_naive) |
           max_pip_naive < 0.01) ~ "Stimulus-responsive",
      (max_pip_naive > 0.5) &
        (is.na(max_pip_stim) | max_pip_stim < 0.01) ~ "Naive-only",
      TRUE ~ NA_character_
    ),
    category = factor(category, levels = rev(c("Shared", "Naive-only", "Stimulus-responsive")))
  )

table(df.max_pip.stim_naive$category)
# Stimulus-responsive          Naive-only              Shared 
#                4626                4151                2498 

pip_bin_breaks <- c(-Inf, 0.01, 0.1, 0.5, 0.9, 1.0)

df.max_pip = dplyr::bind_rows(df.eqtl.max_pip, df.caqtl.max_pip) %>%
  dplyr::full_join(df.max_pip.stim_naive, by = "variant_id") %>%
  dplyr::mutate(
    max_pip_bin = cut(max_pip, pip_bin_breaks),
    max_pip_bin = forcats::fct_recode(max_pip_bin, "[0,0.01]" = "(-Inf,0.01]"),
    max_pip_bin = forcats::fct_na_value_to_level(max_pip_bin, "N.S."),
    max_pip_bin = forcats::fct_relevel(max_pip_bin, "N.S.")
  )

dplyr::filter(df.max_pip, category == "Stimulus-responsive") %>%
  dplyr::mutate(n_total = n()) %>%
  dplyr::group_by(max_pip_bin) %>%
  dplyr::summarize(frac = n() / n_total[1], .groups = "drop")
  

df.rates =
  dplyr::group_by(df.max_pip, QTL, max_pip_bin) %>%
  dplyr::mutate(n_total = n()) %>%
  dplyr::group_by(QTL, max_pip_bin, category) %>%
  dplyr::summarize(
    n_total = n_total[1],
    n_annot = sum(!is.na(category)),
    locusviz::binom_ci(n_annot, n_total),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(category))

pd = position_dodge(width = 0.9)
eqtl.status.colors = c(
  "Shared" = BuenColors::jdb_palette("corona")[5],
  "Stimulus-responsive" = BuenColors::jdb_palette("corona")[3],
  "Naive-only" = "grey50"
)

p.eqtl.rate =
  dplyr::filter(df.rates, QTL == "eQTL" | is.na(QTL)) %>%
  ggplot(aes(frac, max_pip_bin, color = category)) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = frac_lower, xmax = frac_upper), height = 0, position = pd) +
  geom_point(position = pd) +
  scale_x_continuous(expand = expansion(), labels = scales::label_percent()) +
  scale_y_discrete() +
  scale_color_manual(values = eqtl.status.colors) +
  locusviz::get_default_theme(legend.position = c(1, 0.22),
                              legend.justification = c(1, 0.22)) +
  guides(color = guide_legend(reverse = TRUE)) +
  labs(x = "% max PIP > 0.5 (eQTL Catalogue)", y = "FinnGen max PIP bin", color = "eQTL status", title = "eQTL") +
  coord_cartesian(xlim = c(0, 0.07))

p.caqtl.rate =
  dplyr::filter(df.rates, QTL == "caQTL" | is.na(QTL)) %>%
  ggplot(aes(frac, max_pip_bin, color = category)) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = frac_lower, xmax = frac_upper), height = 0, position = pd) +
  geom_point(position = pd) +
  scale_x_continuous(expand = expansion(), labels = scales::label_percent()) +
  scale_y_discrete() +
  scale_color_manual(values = eqtl.status.colors) +
  locusviz::get_default_theme(legend.position = "none", hide.ylab = TRUE) +
  guides(color = guide_legend(reverse = TRUE)) +
  labs(x = "% max PIP > 0.5 (eQTL Catalogue)", y = "Max PIP bin", title = "caQTL") +
  coord_cartesian(xlim = c(0, 0.012))


p.n_stim =
  dplyr::filter(df.max_pip, category == "Stimulus-responsive") %>%
  dplyr::group_by(max_pip_bin, QTL) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>%
  ggplot(aes(n, max_pip_bin, fill = QTL)) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey50") +
  geom_col(position = pd) +
  geom_text(aes(label = scales::comma(n)), size = 2, hjust = -0.5, position = pd) +
  scale_x_continuous(labels = scales::label_comma(),  expand = expansion(c(0, 0.2))) +
  scale_fill_manual(values = qtl.colors, na.value = "grey50", breaks = names(qtl.colors)) +
  locusviz::get_default_theme(hide.ylab = TRUE) +
  labs(x = "# stimulus-responsive eQTLs", y = "Max PIP bin")

plt =  p.eqtl.rate + p.caqtl.rate + p.n_stim +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(tag_levels = "a")
plt

cowplot::save_plot(
  "figures/SFig10_stimulus_responsive_QTL.pdf",
  plt,
  base_height = 60,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig10_stimulus_responsive_QTL.png",
  plt,
  base_height = 60,
  base_width = 180,
  units = "mm",
  dpi = 300
)
