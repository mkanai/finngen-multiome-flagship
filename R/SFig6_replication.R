library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

pip_bin_levels <- c("[0,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]", "(0.9,1]")

read_rep_metrics <- function(path, study) {
  rgsutil::read_gsfiles(path, combine = "row") %>%
    dplyr::mutate(
      level = ifelse(stringr::str_detect(path, "_l1\\."), "L1", "L2"),
      study = sprintf("%s (%s)", .env$study, level),
      pip_bin = dplyr::recode(pip_bin, "(-Inf,0.01]" = "[0,0.01]", "(0,0.01]" = "[0,0.01]"),
      pip_bin = factor(pip_bin, levels = pip_bin_levels)
    ) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::filter(any(n_lead > 0)) %>%
    dplyr::ungroup()
}

read_rho <- function(path, study) {
  rgsutil::read_gsfile(path) %>%
    dplyr::mutate(
      level = ifelse(stringr::str_detect(path, "_l1\\."), "L1", "L2"),
      study = sprintf("%s (%s)", .env$study, level),
      pip_bin = dplyr::recode(pip_bin, "(-Inf,0.01]" = "[0,0.01]", "(0,0.01]" = "[0,0.01]"),
      pip_bin = factor(pip_bin, levels = pip_bin_levels)
    )
}


df.rep.tk <- read_rep_metrics(
  "gs://expansion_areas/multiome/misc/replication/TenK10K/TenK10K_replication_l*.per_pip_bin.tsv",
  study = "TenK10K"
)

df.rep.ok <- read_rep_metrics(
  "gs://expansion_areas/multiome/misc/replication/OneK1K_eQTL_Catalogue/OneK1K_replication_l*.per_pip_bin.tsv",
  study = "OneK1K"
)

df.rep.cima <- read_rep_metrics(
  "gs://expansion_areas/multiome/misc/replication/CIMA/CIMA_eqtl_replication_l2.per_pip_bin.tsv",
  study = "CIMA"
)
df.rep.aida <- read_rep_metrics(
  "gs://expansion_areas/multiome/misc/replication/AIDA/AIDA_eqtl_replication_l2.per_pip_bin.tsv",
  study = "AIDA"
)
df.rep.bbj <- read_rep_metrics(
  "gs://expansion_areas/multiome/misc/replication/bbj/BBJ_eqtl_replication_l1.per_pip_bin.tsv",
  study = "BBJ"
)

df.rho.tk = read_rho(
  "gs://expansion_areas/multiome/misc/replication/TenK10K/TenK10K_replication_l2.rho_per_pip.tsv",
  study = "TenK10K"
)
df.rho.ok = read_rho(
  "gs://expansion_areas/multiome/misc/replication/OneK1K_eQTL_Catalogue/OneK1K_replication_l2.rho_per_pip.tsv",
  study = "OneK1K"
)

df.rho.lead.tk = read_rho(
  "gs://expansion_areas/multiome/misc/replication/TenK10K/TenK10K_replication_l2.rho_per_pip_lead.tsv",
  study = "TenK10K"
)
df.rho.lead.ok = read_rho(
  "gs://expansion_areas/multiome/misc/replication/OneK1K_eQTL_Catalogue/OneK1K_replication_l2.rho_per_pip_lead.tsv",
  study = "OneK1K"
)
df.rho.lead.aida <- read_rho(
  "gs://expansion_areas/multiome/misc/replication/AIDA/AIDA_eqtl_replication_l2.rho_per_pip_lead.tsv",
  study = "AIDA"
)
df.rho.lead.cima <- read_rho(
  "gs://expansion_areas/multiome/misc/replication/CIMA/CIMA_eqtl_replication_l2.rho_per_pip_lead.tsv",
  study = "CIMA"
)
df.rho.lead.bbj <- read_rho(
  "gs://expansion_areas/multiome/misc/replication/bbj/BBJ_eqtl_replication_l1.rho_per_pip_lead.tsv",
  study = "BBJ"
)



base_sig_cols = c(
  "Sign" = "n_concordant_sign",
  "P < 0.05" = "n_sig_0.05",
  "P < 1e-4" = "n_sig_0.0001",
  "P < 1e-6" = "n_sig_1e-06",
  "P < 5e-8" = "n_sig_5e-08",
  "In 95% CS" = "n_cs",
  "PIP > 0.1" = "n_pip"
)
sig.colors = c(
  "Sign" = BuenColors::jdb_palette("brewer_fire")[8],
  "P < 0.05" = BuenColors::jdb_palette("brewer_heat")[7],
  "P < 1e-4" = BuenColors::jdb_palette("brewer_heat")[6],
  "P < 1e-6" = BuenColors::jdb_palette("brewer_heat")[5],
  "P < 5e-8" = BuenColors::jdb_palette("brewer_heat")[4],
  "In 95% CS" = BuenColors::jdb_palette("brewer_purple")[7],
  "PIP > 0.1" = BuenColors::jdb_palette("brewer_purple")[6]
)

studies = c("TenK10K", "OneK1K", "CIMA", "AIDA", "BBJ")
levels = c("L1", "L2")
study.levels = as.vector(outer(levels, studies, function(l, s)
  paste0(s, " (", l, ")")))
base.cols = setNames(BuenColors::jdb_palette("corona")[1:5], studies)
study.cols = setNames(rep(base.cols, each = length(levels)), study.levels)

df.rep = dplyr::bind_rows(df.rep.tk, df.rep.ok, df.rep.aida, df.rep.cima, df.rep.bbj) %>%
  dplyr::mutate(study = factor(study, levels = rev(study.levels)))

df.rho = dplyr::bind_rows(
  df.rho.lead.tk,
  df.rho.lead.ok,
  df.rho.lead.aida,
  df.rho.lead.cima,
  df.rho.lead.bbj
) %>%
  dplyr::mutate(study = factor(study, levels = rev(study.levels)))

munge_rep_rate = function(df,
                          sig_cols = base_sig_cols,
                          denom = "n") {
  if (is.null(names(sig_cols))) {
    names(sig_cols) = sig_cols
  }
  
  tidyr::pivot_longer(
    df,
    cols = dplyr::all_of(unname(sig_cols)),
    names_to = "threshold",
    values_to = "n_sig"
  ) %>%
    dplyr::mutate(threshold =
                    forcats::fct_recode(threshold, !!!sig_cols) %>%
                    forcats::fct_relevel(rev(names(sig_cols)))) %>%
    dplyr::group_by(study, pip_bin, threshold) %>%
    dplyr::summarise(n_sig = sum(n_sig),
                     n = sum(.data[[denom]]),
                     .groups = "drop") %>%
    dplyr::mutate(locusviz::binom_ci(n_sig, n))
}

plot_rep_rate = function(df, hide.ylab = FALSE) {
  pd = position_dodge(width = 0.9)
  
  ggplot(df, aes(frac, pip_bin, color = threshold)) +
    geom_errorbar(aes(xmin = frac_lower, xmax = frac_upper),
                  width = 0,
                  position = pd) +
    geom_point(position = pd) +
    locusviz::get_default_theme(hide.ylab = hide.ylab) +
    theme(plot.title = element_text(hjust = 0.02)) +
    scale_x_continuous(labels = scales::label_percent()) +
    scale_color_manual(values = sig.colors) +
    guides(color = guide_legend(reverse = TRUE)) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "% variants", y = "FinnGen PIP bin", color = "Replication")
}

plot_lead_overlap = function(df) {
  ggplot(df, aes(frac, pip_bin, color = study)) +
    geom_errorbar(aes(xmin = frac_lower, xmax = frac_upper),
                  width = 0,
                  position = pd) +
    geom_point(position = pd) +
    locusviz::get_default_theme(legend.position = c(1, 0),
                                legend.justification = c(1, 0)) +
    theme(plot.title = element_text(hjust = 0.02)) +
    scale_x_continuous(labels = scales::label_percent()) +
    scale_color_manual(values = study.cols) +
    guides(color = guide_legend(reverse = TRUE)) +
    labs(x = "% lead variants", y = "FinnGen PIP bin", color = "Study")
}

plot_sign_concordance = function(df) {
  ggplot(df, aes(frac, pip_bin, color = study)) +
    geom_errorbar(aes(xmin = frac_lower, xmax = frac_upper),
                  width = 0,
                  position = pd) +
    geom_point(position = pd) +
    locusviz::get_default_theme(legend.position = "none", hide.ylab = TRUE) +
    scale_x_continuous(labels = scales::label_percent()) +
    scale_color_manual(values = study.cols) +
    guides(color = guide_legend(reverse = TRUE)) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "% sign concordance", y = "FinnGen PIP bin", color = "Study")
}

plot_rho = function(df) {
  ggplot(df, aes(rho, pip_bin, color = study)) +
    geom_errorbar(aes(xmin = rho_lower, xmax = rho_upper),
                  width = 0,
                  position = pd) +
    geom_point(position = pd) +
    locusviz::get_default_theme(legend.position = "none", hide.ylab = TRUE) +
    scale_x_continuous() +
    scale_color_manual(values = study.cols) +
    guides(color = guide_legend(reverse = TRUE)) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = expression(paste("Effect correlation ", italic(rho))),
         y = "FinnGen PIP bin",
         color = "Study")
}

pd = position_dodge(width = 0.9)

p.tk_l1 =
  dplyr::filter(df.rep.tk, filter_l1_cell_types(cell_type)) %>%
  munge_rep_rate() %>%
  plot_rep_rate() +
  labs(title = "TenK10K eQTL (L1)")

p.tk_l2 =
  dplyr::filter(df.rep.tk, !filter_l1_cell_types(cell_type)) %>%
  munge_rep_rate() %>%
  plot_rep_rate(hide.ylab = TRUE) +
  labs(title = "TenK10K eQTL (L2)")

p.ok_l1 =
  dplyr::filter(df.rep.ok, filter_l1_cell_types(cell_type)) %>%
  munge_rep_rate() %>%
  plot_rep_rate() +
  labs(title = "OneK1K eQTL (L1)")

p.ok_l2 =
  dplyr::filter(df.rep.ok, !filter_l1_cell_types(cell_type)) %>%
  munge_rep_rate() %>%
  plot_rep_rate(hide.ylab = TRUE) +
  labs(title = "OneK1K eQTL (L2)")

p.lead =
  dplyr::filter(df.rep, level == "L2" |
                  study == "BBJ (L1)") %>%
  munge_rep_rate(sig_cols = "n_lead") %>%
  plot_lead_overlap() +
  coord_cartesian(xlim = c(0, 0.35)) +
  labs(title = "Lead eQTL variants")

p.sign =
  dplyr::filter(df.rep, level == "L2" |
                  study == "BBJ (L1)") %>%
  munge_rep_rate("n_lead_concordant_sign", denom = "n_lead") %>%
  plot_sign_concordance()

p.rho = plot_rho(df.rho)

#################################

df.rep.cima.caQTL <- read_rep_metrics(
  "gs://expansion_areas/multiome/misc/replication/CIMA/CIMA_caqtl_replication_l2.per_pip_bin.tsv",
  study = "CIMA"
)
df.rep.bbj.caQTL <- read_rep_metrics(
  "gs://expansion_areas/multiome/misc/replication/bbj/BBJ_caqtl_replication_l1.per_pip_bin.tsv",
  study = "BBJ"
)

df.rho.lead.cima.caQTL <- read_rho(
  "gs://expansion_areas/multiome/misc/replication/CIMA/CIMA_caqtl_replication_l2.rho_per_pip_lead.tsv",
  study = "CIMA"
)
df.rho.lead.bbj.caQTL <- read_rho(
  "gs://expansion_areas/multiome/misc/replication/bbj/BBJ_caqtl_replication_l1.rho_per_pip_lead.tsv",
  study = "BBJ"
)

df.rep.caQTL = dplyr::bind_rows(df.rep.cima.caQTL, df.rep.bbj.caQTL) %>%
  dplyr::mutate(study = factor(study, levels = rev(study.levels)))

df.rho.caQTL = dplyr::bind_rows(df.rho.lead.cima.caQTL, df.rho.lead.bbj.caQTL) %>%
  dplyr::mutate(study = factor(study, levels = rev(study.levels)))


p.lead.caQTL =
  munge_rep_rate(df.rep.caQTL, sig_cols = "n_lead") %>%
  plot_lead_overlap() +
  coord_cartesian(xlim = c(0, 0.07)) +
  labs(title = "Lead caQTL variants")

p.sign.caQTL =
  munge_rep_rate(df.rep.caQTL, "n_lead_concordant_sign", denom = "n_lead") %>%
  plot_sign_concordance()

p.rho.caQTL = plot_rho(df.rho.caQTL)

list(p.tk_l1, p.tk_l2, p.ok_l1, p.ok_l2) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 2, guides = "collect") +
  patchwork::plot_annotation(tag_levels = "a")

list(p.lead, p.sign, p.rho, p.lead.caQTL, p.sign.caQTL, p.rho.caQTL) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout() +
  patchwork::plot_annotation(tag_levels = "a")


################################################################################
df.cross <- rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/TenK10K/TenK10K_cross_celltype.tsv.bgz"
)

df.cross_matrix <-
  dplyr::group_by(df.cross, fg_cell_type, tk_cell_type) %>%
  dplyr::filter(dplyr::n() > 20) %>%
  dplyr::summarize(
    n = dplyr::n(),
    rho = cor(fg_beta, tk_beta, method = "spearman"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(fg_level = ifelse(filter_l1_cell_types(fg_cell_type), "l1", "l2"))

tk_cell_type_labels <- c(
  "ASDC" = "ASDC",
  "B_intermediate" = "B Intermediate",
  "B_memory" = "B Memory",
  "B_naive" = "B Naive",
  "Plasmablast" = "Plasmablast",
  "CD4_CTL" = "CD4+ CTL",
  "CD4_Naive" = "CD4+ Naive",
  "CD4_Proliferating" = "CD4+ Proliferating",
  "CD4_TCM" = "CD4+ TCM",
  "CD4_TEM" = "CD4+ TEM",
  "Treg" = "Treg",
  "CD8_Naive" = "CD8+ Naive",
  "CD8_Proliferating" = "CD8+ Proliferating",
  "CD8_TCM" = "CD8+ TCM",
  "CD8_TEM" = "CD8+ TEM",
  "MAIT" = "MAIT",
  "CD14_Mono" = "CD14+ Mono",
  "CD16_Mono" = "CD16+ Mono",
  "cDC1" = "cDC1",
  "cDC2" = "cDC2",
  "pDC" = "pDC",
  "NK" = "NK",
  "NK_CD56bright" = "NK CD56bright",
  "NK_Proliferating" = "NK Proliferating",
  "HSPC" = "HSPC",
  "ILC" = "ILC",
  "dnT" = "dnT",
  "gdT" = "gdT"
)

plot_cross_rho = function(df.cross) {
  rho_wide_fg <-
    tidyr::pivot_wider(
      df.cross,
      id_cols = fg_cell_type,
      names_from = tk_cell_type,
      values_from = rho,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("fg_cell_type")
  fg_order <- rownames(rho_wide_fg)[hclust(dist(rho_wide_fg))$order]
  
  rho_wide_tk <-
    tidyr::pivot_wider(
      df.cross,
      id_cols = tk_cell_type,
      names_from = fg_cell_type,
      values_from = rho,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("tk_cell_type")
  tk_order <- rownames(rho_wide_tk)[hclust(dist(rho_wide_tk))$order]
  
  dplyr::mutate(
    df.cross,
    fg_cell_type = factor(fg_cell_type, levels = rev(fg_order)),
    tk_cell_type = factor(tk_cell_type, levels = tk_order)
  ) %>%
    ggplot(aes(tk_cell_type, fg_cell_type, fill = rho)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = BuenColors::jdb_palette("solar_extra"),
      limits = c(0, 1),
      oob = scales::squish,
      labels = scales::label_number(drop0trailing = TRUE)
    ) +
    scale_x_discrete(labels = tk_cell_type_labels) +
    scale_y_discrete(labels = l2.labels) +
    locusviz::get_default_theme(legend.position = "bottom",
                                legend.justification = "bottom") +
    theme(
      legend.key.width = unit(5, "mm"),
      legend.title = element_text(margin = margin(
        t = 22, b = 1, unit = "mm"
      )),
      legend.title.position = "top"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "TenK10K",
         y = "FinnGen (L2)",
         fill = expression(paste("Effect correlation ", italic(rho))))
}

p.cross.l2 =
  dplyr::filter(df.cross_matrix, !filter_l1_cell_types(fg_cell_type)) %>%
  dplyr::mutate(fg_cell_type = remove_cell_type_prefix(fg_cell_type)) %>%
  plot_cross_rho()

layout = "
AB#
CDE
FFF
"

plt =
  list(p.tk_l1,
       p.tk_l2,
       p.ok_l1,
       p.ok_l2,
       patchwork::guide_area(),
       (patchwork::free(p.cross.l2, side = "l"))) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(
    design = layout,
    guides = "collect",
    widths = c(2, 2, 1),
    heights = c(2, 2, 3)
  ) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig8_TenK10K_OneK1K_replication.pdf",
  plt,
  base_height = 170,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig8_TenK10K_OneK1K_replication.png",
  plt,
  base_height = 170,
  base_width = 180,
  units = "mm",
  dpi = 300
)

plt.lead =
  list(p.lead, p.sign, p.rho, p.lead.caQTL, p.sign.caQTL, p.rho.caQTL) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout() +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig9_lead_replication.pdf",
  plt.lead,
  base_height = 120,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig9_lead_replication.png",
  plt.lead,
  base_height = 120,
  base_width = 180,
  units = "mm",
  dpi = 300
)

################################################################################
bbj_to_fg_l1 <- c(
  B = "predicted.celltype.l1.B",
  CD4T = "predicted.celltype.l1.CD4_T",
  CD8T = "predicted.celltype.l1.CD8_T",
  DC = "predicted.celltype.l1.DC",
  Mono = "predicted.celltype.l1.Mono",
  NK = "predicted.celltype.l1.NK"
)
df.bbj.peak_map <- rgsutil::read_gsfile("gs://expansion_areas/multiome/misc/replication/bbj/bbj_to_fg_peak_map.tsv.gz")

df.bbj.open4gene = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/misc/replication/bbj/*.open4gene.tsv.gz", function(df, path) {
  cell_type = stringr::str_replace(path, ".*/(.*)\\.open4gene\\.tsv\\.gz", "\\1")
  dplyr::mutate(
    df,
    cell_type = bbj_to_fg_l1[cell_type],
    bbj_peak_id = stringr::str_replace_all(peak_id, "-", "_")
  )
})  %>%
  tidyr::drop_na(hurdle_zero_nlog10p, hurdle_count_nlog10p) %>%
  dplyr::left_join(df.features, by = c("gene_id" = "symbol")) %>%
  dplyr::left_join(df.bbj.peak_map, by = "bbj_peak_id") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(
    fasthurdle::joint_score_test(10 ** -hurdle_zero_nlog10p, 10 ** -hurdle_count_nlog10p, alpha = 0.05)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(peak_id = fg_peak_id, gene_id = phenotype_id)

df.open4gene = readRDS("data/open4gene.rds")
df.open4gene.sig = dplyr::filter(df.open4gene, sig)

open4gene.pval.levels = c("P < 0.05", "P < 1e-4", "P < 5e-8", "P < 1e-10", "P < 1e-20")
open4gene.pval.colors = c(
  "P < 0.05" = BuenColors::jdb_palette("brewer_heat")[8],
  "P < 1e-4" = BuenColors::jdb_palette("brewer_heat")[7],
  "P < 5e-8" = BuenColors::jdb_palette("brewer_heat")[6],
  "P < 1e-10" = BuenColors::jdb_palette("brewer_heat")[5],
  "P < 1e-20" = BuenColors::jdb_palette("brewer_heat")[4]
)


df.open4gene.merged =
  dplyr::inner_join(df.open4gene.sig,
                    df.bbj.open4gene,
                    by = c("gene_id", "peak_id", "cell_type")) %>%
  dplyr::mutate(
    fg_pval_bin = dplyr::case_when(
      p_joint.x < 1e-20 ~ "P < 1e-20",
      p_joint.x < 1e-10 ~ "P < 1e-10",
      p_joint.x < 5e-8 ~ "P < 5e-8",
      p_joint.x < 1e-4 ~ "P < 1e-4",
      p_joint.x < 0.05 ~ "P < 0.05",
      TRUE ~ NA_character_
    ),
    fg_pval_bin = factor(fg_pval_bin, levels = open4gene.pval.levels)
  )

dplyr::count(df.open4gene.merged, fg_pval_bin, mode.x, mode.y) %>% clipr::write_clip()

df.open4gene.merged.rate =
  dplyr::group_by(df.open4gene.merged, fg_pval_bin) %>%
  dplyr::reframe(
    n_total = n(),
    n_count = c(
      sum(p_joint.y < 0.05),
      sum(p_joint.y < 1e-4),
      sum(p_joint.y < 5e-8),
      sum(p_joint.y < 1e-10),
      sum(p_joint.y < 1e-20)
    ),
    threshold = c("P < 0.05", "P < 1e-4", "P < 5e-8", "P < 1e-10", "P < 1e-20"),
    threshold = factor(threshold, levels = rev(open4gene.pval.levels))
  ) %>%
  dplyr::group_by(threshold) %>%
  dplyr::mutate(locusviz::binom_ci(n_count, n_total))

df.open4gene.merged.rho =
  dplyr::group_by(df.open4gene.merged, fg_pval_bin) %>%
  dplyr::summarize(
    locusviz::binom_ci(sum(
      hurdle_zero_beta.x * hurdle_zero_beta.y > 0
    ), n(), colname = "frac_zero"),
    locusviz::binom_ci(
      sum(hurdle_count_beta.x * hurdle_count_beta.y > 0),
      n(),
      colname = "frac_count"
    ),
    locusviz::spearman_ci(hurdle_zero_beta.x, hurdle_zero_beta.y, colname = "rho_zero"),
    locusviz::spearman_ci(hurdle_count_beta.x, hurdle_count_beta.y, colname = "rho_count"),
    .groups = "drop"
  ) %>%
  dplyr::rename_with(~ sub("^(rho_(zero|count))$", "\\1_value", .x)) %>%
  dplyr::rename_with(~ sub("^(frac_(zero|count))$", "\\1_value", .x)) %>%
  tidyr::pivot_longer(
    cols = tidyselect::matches("^(rho|frac)_(zero|count)_"),
    names_to = c("metric", "model", ".value"),
    names_pattern = "(rho|frac)_(zero|count)_(value|lower|upper|p)"
  ) %>%
  tidyr::pivot_wider(
    names_from = metric,
    values_from = c(value, lower, upper, p),
    names_glue = "{metric}_{.value}"
  ) %>%
  dplyr::rename(rho = rho_value, frac = frac_value) %>%
  dplyr::mutate(model = factor(model, levels = rev(c("zero", "count"))))


df.open4gene.merged.mode <-
  dplyr::filter(
    df.open4gene.merged,
    mode.x %in% c("rheostat", "dual", "switch"),
    mode.y %in% c("rheostat", "dual", "switch"),
    !is.na(fg_pval_bin)
  ) %>%
  dplyr::mutate(mode.x = factor(mode.x, levels = c("rheostat", "dual", "switch")),
                mode.y = factor(mode.y, levels = c("rheostat", "dual", "switch"))) %>%
  dplyr::count(fg_pval_bin, mode.x, mode.y, .drop = FALSE) %>%
  dplyr::group_by(fg_pval_bin, mode.x) %>%
  dplyr::mutate(n_total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_total > 0) %>%
  dplyr::mutate(locusviz::binom_ci(n, n_total))

pd = position_dodge(width = 0.9)

p.open4gene.sig =
  ggplot(df.open4gene.merged.rate,
         aes(frac, fg_pval_bin, color = threshold)) +
  geom_errorbar(aes(xmin = frac_lower, xmax = frac_upper),
                width = 0,
                position = pd) +
  geom_point(position = pd) +
  locusviz::get_default_theme(legend.position = c(1, 0),
                              legend.justification = c(1, 0)) +
  theme(plot.title = element_text(hjust = 0.02)) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_color_manual(values = open4gene.pval.colors) +
  guides(color = guide_legend(reverse = TRUE)) +
  coord_cartesian(xlim = c(0, 0.61)) +
  labs(
    x = "% peak-gene pairs",
    y = expression(paste("FinnGen ", italic(P)[joint], " bin")),
    color = expression(paste("BBJ ", italic(P)[joint])),
    title = "Peak-gene links"
  )



p.open4gene.sign =
  ggplot(df.open4gene.merged.rho, aes(frac, fg_pval_bin, color = model)) +
  geom_errorbar(aes(xmin = frac_lower, xmax = frac_upper),
                width = 0,
                position = pd) +
  geom_point(position = pd) +
  locusviz::get_default_theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    hide.ylab = TRUE
  ) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_color_manual(values = model.cols) +
  guides(color = guide_legend(reverse = TRUE)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "% sign concordance",
       y = expression(paste("FinnGen ", italic(P)[joint], " bin")),
       color = "Model")

p.open4gene.rho =
  ggplot(df.open4gene.merged.rho, aes(rho, fg_pval_bin, color = model)) +
  geom_errorbar(aes(xmin = rho_lower, xmax = rho_upper),
                width = 0,
                position = pd) +
  geom_point(position = pd) +
  locusviz::get_default_theme(legend.position = "none", hide.ylab = TRUE) +
  scale_x_continuous() +
  scale_color_manual(values = model.cols) +
  guides(color = guide_legend(reverse = TRUE)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = expression(paste("Effect correlation ", italic(rho))),
       y = expression(paste("FinnGen ", italic(P)[joint], " bin")),
       color = "Model")

plot_mode_concordance <- function(df,
                                  fg_mode,
                                  hide.ylab = TRUE,
                                  legend.position = "none") {
  dplyr::filter(df, mode.x == fg_mode) %>%
    ggplot(aes(frac, fg_pval_bin, color = mode.y)) +
    geom_errorbar(aes(xmin = frac_lower, xmax = frac_upper),
                  width = 0,
                  position = pd) +
    geom_point(position = pd) +
    locusviz::get_default_theme(
      legend.position = legend.position,
      legend.justification = legend.position,
      hide.ylab = hide.ylab
    ) +
    scale_x_continuous(labels = scales::label_percent()) +
    scale_color_manual(values = link_mode_colors, labels = stringr::str_to_title) +
    guides(color = guide_legend(reverse = TRUE)) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      x = "% BBJ-replicating links",
      y = expression(paste("FinnGen ", italic(P)[joint], " bin")),
      color = "BBJ mode",
      title = sprintf("FinnGen %s", fg_mode)
    )
}

p.mode.switch   <- plot_mode_concordance(
  df.open4gene.merged.mode,
  "switch",
  hide.ylab = FALSE,
  legend.position =  c(1, 0)
)
p.mode.rheostat <- plot_mode_concordance(df.open4gene.merged.mode, "rheostat")
p.mode.dual     <- plot_mode_concordance(df.open4gene.merged.mode, "dual")

plt.open4gene =
  list(
    p.open4gene.sig,
    p.open4gene.sign,
    p.open4gene.rho,
    p.mode.switch,
    p.mode.rheostat,
    p.mode.dual
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout() +
  patchwork::plot_annotation(tag_levels = "a")

plt.open4gene

cowplot::save_plot(
  "figures/SFig6_peak_gene_replication.pdf",
  plt.open4gene,
  base_height = 120,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig6_peak_gene_replication.png",
  plt.open4gene,
  base_height = 120,
  base_width = 180,
  units = "mm",
  dpi = 300
)
