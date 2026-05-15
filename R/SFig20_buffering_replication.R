library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))


df.trait = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/metadata/phenotypes/FG_R12_endpoints.KANTA.summary.txt"
) %>%
  dplyr::filter(H2_Z > 2 & num_gw_significant > 10) %>%
  dplyr::select(-sumstats, -ldsc_munged) %>%
  dplyr::rename(
    ldsc_h2 = H2,
    ldsc_h2_se = SE,
    ldsc_h2_z = H2_Z,
    ldsc_int = INT,
    ldsc_int_se = INT_SE,
    ldsc_ratio = RATIO,
    ldsc_ratio_se = RATIO_SE
  )

df.loeuf.v4.rep = rgsutil::read_gsfile(
  "gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
) %>%
  dplyr::inner_join(df.features, by = c("gene_id" = "phenotype_id")) %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile)


df.coloc.fg = readRDS("data/df.coloc.all.rds") %>%
  dplyr::filter(trait1 %in% df.trait$phenocode & QTL == "eQTL")
df.acat.fg = data.table::fread("data/df.gex.with_lead_beta.tsv.gz", data.table = FALSE)

df.fg.n_egenes = dplyr::filter(df.acat.fg, ACAT_q < 0.05) %>%
  dplyr::distinct(phenotype_id) %>%
  dplyr::summarize(study = "FinnGen", n_egenes = n())

df.fg.l2.n_egenes = dplyr::filter(df.acat.fg, ACAT_q < 0.05 &
                                    !filter_l1_cell_types(cell_type)) %>%
  dplyr::distinct(phenotype_id, cell_type) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarize(study = "FinnGen", n_egenes = n())

df.acat.tk = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/TenK10K/TenK10K_ACAT_lead_with_beta.tsv.bgz"
)
df.tk.n_egenes = dplyr::filter(df.acat.tk, ACAT_q < 0.05) %>%
  dplyr::distinct(gene) %>%
  dplyr::summarize(study = "TenK10K",
                   n_egenes = n(),
                   .groups = "drop")
df.tk.l2.n_egenes = dplyr::filter(df.acat.tk, ACAT_q < 0.05) %>%
  dplyr::distinct(cell_type, gene) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarize(study = "TenK10K",
                   n_egenes = n(),
                   .groups = "drop")

df.coloc.tk.l2 = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/TenK10K/coloc_100kb/TenK10K_coloc.l2.tsv.bgz"
)

df.coloc.eqtl_catalogue = readRDS("data/df.coloc.eqtl_catalogue.rds")
valid_dataset_ids = readRDS("data/valid_dataset_ids.rds")

df.eq.permuted = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/eQTL_Catalogue/eQTL_Catalogue.r8_beta.ge.permuted.tsv.gz"
) %>%
  tidyr::drop_na(molecular_trait_id) %>%
  dplyr::filter(dataset_id %in% valid_dataset_ids)

df.eq.n_egenes =
  dplyr::group_by(df.eq.permuted, dataset_id) %>%
  dplyr::summarize(n_egenes = sum(qval < 0.05), .groups = "drop")

if (!file.exists("data/df.sb.meta.rds")) {
  df.sb.meta = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/misc/replication/SingleBrain/meta_analysis/*_eqtl_top_assoc.tsv.gz", function(df, path) {
    dplyr::mutate(
      df,
      cell_type = stringr::str_replace(basename(path), "_eqtl_top_assoc.tsv.gz", ""),
      feature = stringr::str_remove(feature, "\\.[0-9]+$")
    ) %>%
      dplyr::select(cell_type, feature, variant_id, fixed_beta, Fixed_P, qval)
  }) %>%
    dplyr::rename(gene_id = feature, beta = fixed_beta)
  saveRDS(df.sb.meta, "data/df.sb.meta.rds")
}
if (!file.exists("data/df.sb.saige.rds")) {
  df.sb.saige = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/misc/replication/SingleBrain/SAIGE_QTL/*_top_assoc.tsv.gz", function(df, path) {
    dplyr::mutate(
      df,
      cell_type = stringr::str_replace(basename(path), "_top_assoc.tsv.gz", ""),
      feature = stringr::str_remove(feature, "\\.[0-9]+$")
    ) %>%
      dplyr::mutate(qval = qvalue::qvalue(ACAT_p)$qvalues)
  }) %>%
    dplyr::rename(gene_id = feature, beta = BETA)
  saveRDS(df.sb.saige, "data/df.sb.saige.rds")
}


df.sb.meta = readRDS("data/df.sb.meta.rds")
df.sb.meta.n_egenes =
  dplyr::filter(df.sb.meta, qval < 0.05) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarize(study = "SingleBrain (meta)",
                   n_egenes = n(),
                   .groups = "drop")

df.sb.saige = readRDS("data/df.sb.saige.rds")
df.sb.saige.n_egenes =
  dplyr::filter(df.sb.saige, qval < 0.05) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarize(study = "SingleBrain (SAIGE-QTL)",
                   n_egenes = n(),
                   .groups = "drop")

df.fg.rho =
  dplyr::filter(df.acat.fg, ACAT_q < 0.05) %>%
  dplyr::inner_join(df.loeuf.v4.rep, by = c("phenotype_id" = "gene_id")) %>%
  dplyr::left_join(df.fg.l2.n_egenes, by = "cell_type") %>%
  dplyr::group_by(cell_type, n_egenes) %>%
  dplyr::summarize(
    study = "FinnGen",
    method = "SAIGE-QTL",
    locusviz::spearman_ci(abs(lead_beta), lof.oe_ci.upper),
    .groups = "drop"
  )

df.tk.rho =
  dplyr::filter(df.acat.tk, ACAT_q < 0.05) %>%
  dplyr::inner_join(df.loeuf.v4.rep, by = c("gene" = "gene_id")) %>%
  dplyr::left_join(df.tk.l2.n_egenes, by = "cell_type") %>%
  dplyr::group_by(cell_type, n_egenes) %>%
  dplyr::summarize(
    study = "TenK10K",
    method = "SAIGE-QTL",
    locusviz::spearman_ci(abs(lead_beta), lof.oe_ci.upper),
    .groups = "drop"
  )

df.sb.saige.rho =
  dplyr::filter(df.sb.saige, qval < 0.05) %>%
  dplyr::inner_join(df.loeuf.v4.rep, by = "gene_id") %>%
  dplyr::left_join(df.sb.saige.n_egenes, by = "cell_type") %>%
  dplyr::group_by(cell_type, n_egenes) %>%
  dplyr::summarize(
    study = "SingleBrain",
    method = "SAIGE-QTL",
    locusviz::spearman_ci(abs(beta), lof.oe_ci.upper),
    .groups = "drop"
  )


df.eq.beta.rho =
  dplyr::filter(df.eq.permuted, qval < 0.05) %>%
  dplyr::rename(gene_id = molecular_trait_id) %>%
  dplyr::inner_join(df.loeuf.v4.rep, by = "gene_id") %>%
  dplyr::left_join(df.eq.n_egenes, by = "dataset_id") %>%
  dplyr::filter(n_egenes > 100) %>%
  dplyr::group_by(dataset_id, n_egenes) %>%
  dplyr::summarize(
    study = "eQTL Catalogue",
    method = "QTLtools (RINT)",
    locusviz::spearman_ci(abs(beta), lof.oe_ci.upper),
    .groups = "drop"
  )

study.colors = c(
  "FinnGen" = BuenColors::jdb_palette("corona")[14],
  "FinnGen (L2)" = BuenColors::jdb_palette("corona")[14],
  "TenK10K" = BuenColors::jdb_palette("corona")[9],
  "eQTL Catalogue" = BuenColors::jdb_palette("corona")[4],
  "SingleBrain" = BuenColors::jdb_palette("corona")[10]
)

df.beta.rho = dplyr::bind_rows(df.eq.beta.rho, df.fg.rho, df.tk.rho, df.sb.saige.rho) %>%
  dplyr::filter(n_egenes > 100)

dplyr::summarize(
  df.beta.rho,
  n_total = n(),
  n_positive = sum(rho > 0),
  frac_positive = mean(rho > 0),
  median_rho = median(rho),
  min_rho = min(rho),
  max_rho = max(rho)
)

dplyr::filter(df.beta.rho, study != "eQTL Catalogue") %>% clipr::write_clip()

df.beta.rho.label =
  dplyr::filter(
    df.beta.rho,
    study %in% c("FinnGen", "TenK10K", "SingleBrain") |
      dataset_id == "QTD001000"
  ) %>%
  dplyr::group_by(study) %>%
  dplyr::filter(n_egenes == max(n_egenes)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label = ifelse(study == "eQTL Catalogue", "INTERVAL", study),
    hjust = ifelse(study == "eQTL Catalogue", 0.75, -0.25),
    vjust = dplyr::case_when(study == "eQTL Catalogue" ~ 3, study == "SingleBrain" ~ -3, TRUE ~ 0.5),
  )

p.beta.rho =
  dplyr::mutate(
    df.beta.rho,
    w = 1 / (rho_upper - rho_lower)^2,
    study = factor(study, level = names(study.colors)),
    method = factor(method, levels = c("SAIGE-QTL", "QTLtools (RINT)"))
  ) %>%
  dplyr::arrange(dplyr::desc(study)) %>%
  ggplot(aes(n_egenes, rho, color = study, shape = method)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_smooth(
    aes(weight = w, fill = study),
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,
    alpha = 0.1
  ) +
  geom_point() +
  geom_text(
    aes(
      label = label,
      hjust = hjust,
      vjust = vjust
    ),
    size = 2,
    data = df.beta.rho.label,
    show.legend = FALSE
  ) +
  coord_cartesian(ylim = c(-0.05, 0.4)) +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous(expand = expansion()) +
  scale_color_manual(values = study.colors) +
  scale_fill_manual(values = study.colors) +
  scale_shape_manual(values = c("SAIGE-QTL" = 17, "QTLtools (RINT)" = 16)) +
  locusviz::get_default_theme(legend.position = c(1, 0.1),
                              legend.justification = c(1, 0.1)) +
  labs(
    x = "# eGenes",
    y = expression(paste(
      "|", italic(beta)[eQTL], "|-LOEUF Spearman's ", italic(rho)
    )),
    color = "Study",
    shape = "Method",
    fill = "Study"
  ) +
  guides(
    color = guide_legend(order = 1),
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )
p.beta.rho


df.coloc.tk.rate =
  dplyr::distinct(df.acat.tk, gene) %>%
  dplyr::left_join(
    dplyr::filter(df.coloc.tk.l2, PP.H4.abf > 0.8) %>%
      dplyr::distinct(gene) %>%
      dplyr::mutate(coloc = TRUE),
    by = "gene"
  ) %>%
  dplyr::mutate(coloc = !is.na(coloc)) %>%
  dplyr::rename(gene_id = gene) %>%
  dplyr::left_join(df.loeuf.v4.rep, by = "gene_id") %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile) %>%
  dplyr::group_by(lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(
    study = "TenK10K",
    dataset_id = "TenK10K",
    method = "SAIGE-QTL",
    n_coloc = sum(coloc),
    n_total = n(),
    locusviz::binom_ci(n_coloc, n_total),
    .groups = "drop"
  ) %>%
  dplyr::left_join(df.tk.n_egenes, by = "study")

df.coloc.fg.rate =
  dplyr::distinct(df.acat.fg, phenotype_id) %>%
  dplyr::left_join(
    dplyr::distinct(df.coloc.fg, trait2) %>%
      dplyr::mutate(coloc = TRUE),
    by = c("phenotype_id" = "trait2")
  ) %>%
  dplyr::mutate(coloc = !is.na(coloc)) %>%
  dplyr::rename(gene_id = phenotype_id) %>%
  dplyr::left_join(df.loeuf.v4.rep, by = "gene_id") %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile) %>%
  dplyr::group_by(lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(
    study = "FinnGen",
    dataset_id = "FinnGen",
    method = "SAIGE-QTL",
    n_coloc = sum(coloc),
    n_total = n(),
    locusviz::binom_ci(n_coloc, n_total),
    .groups = "drop"
  ) %>%
  dplyr::left_join(df.fg.n_egenes, by = "study")

df.coloc.eq.rate =
  dplyr::distinct(df.eq.permuted, dataset_id, molecular_trait_id) %>%
  dplyr::left_join(
    dplyr::filter(df.coloc.eqtl_catalogue, PP.H4.abf > 0.8) %>%
      dplyr::distinct(dataset_id, trait2) %>%
      dplyr::mutate(coloc = TRUE),
    by = c("dataset_id", "molecular_trait_id" = "trait2")
  ) %>%
  dplyr::mutate(coloc = !is.na(coloc)) %>%
  dplyr::rename(gene_id = molecular_trait_id) %>%
  dplyr::inner_join(df.loeuf.v4.rep, by = "gene_id") %>%
  dplyr::group_by(dataset_id, lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(
    study = "eQTL Catalogue",
    method = "QTLtools (RINT)",
    n_coloc = sum(coloc),
    n_total = n(),
    locusviz::binom_ci(n_coloc, n_total),
    .groups = "drop"
  ) %>%
  dplyr::left_join(df.eq.n_egenes, by = "dataset_id")

dplyr::filter(df.coloc.eq.rate, n_egenes > 100 &
                round(frac, 2) == 0.01) %>%
  summary()



df.coloc.rate =
  dplyr::bind_rows(df.coloc.eq.rate, df.coloc.tk.rate, df.coloc.fg.rate)

df.coloc.rate.rho =
  dplyr::group_by(df.coloc.rate, study, dataset_id, n_egenes, method) %>%
  dplyr::summarize(
    n_coloc = sum(n_coloc),
    locusviz::spearman_ci(frac, lof.oe_ci.upper_bin_decile),
    .groups = "drop"
  )

df.coloc.rate.rho.label =
  dplyr::filter(
    df.coloc.rate.rho,
    study %in% c("FinnGen", "TenK10K", "SingleBrain") |
      dataset_id == "QTD001000"
  ) %>%
  dplyr::group_by(study) %>%
  dplyr::filter(n_egenes == max(n_egenes)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label = ifelse(study == "eQTL Catalogue", "INTERVAL", study),
    hjust = ifelse(study == "eQTL Catalogue", -0.15, 0.65),
    vjust = ifelse(study == "eQTL Catalogue", -0.5, 2),
  )


df.coloc.rate.hlines =
  dplyr::filter(df.coloc.rate.rho, study != "eQTL Catalogue")

p.coloc.rate.rho =
  dplyr::filter(df.coloc.rate.rho, n_egenes > 100) %>%
  dplyr::mutate(w = 1 / (rho_upper - rho_lower)^2,
                study = factor(study, level = names(study.colors))) %>%
  ggplot(aes(n_egenes, rho, color = study, shape = method)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_smooth(
    aes(weight = w, fill = study),
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,
    alpha = 0.1
  ) +
  geom_point() +
  geom_text(aes(
    label = label,
    hjust = hjust,
    vjust = vjust
  ),
  size = 2,
  data = df.coloc.rate.rho.label) +
  # geom_hline(aes(yintercept = rho, color = study),
  #            linetype = "dashed",
  #            data = df.coloc.rate.hlines) +
  # geom_errorbar(aes(ymin = rho_lower, ymax = rho_upper, color = study), width = 0) +
  coord_cartesian(ylim = c(-1.25, 1)) +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous(expand = expansion()) +
  scale_color_manual(values = study.colors) +
  scale_fill_manual(values = study.colors, guide = "none") +
  scale_shape_manual(values = c("SAIGE-QTL" = 17, "QTLtools (RINT)" = 16)) +
  locusviz::get_default_theme(legend.position = "none") +
  labs(
    x = "# eGenes",
    y = expression(paste("% coloc-LOEUF Spearman's ", italic(rho))),
    color = "Study",
    fill = "Study"
  ) +
  guides(color = guide_legend(order = 1))

plt =
  p.beta.rho + p.coloc.rate.rho +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig20_buffering_replication.pdf",
  plt,
  base_height = 60,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig20_buffering_replication.png",
  plt,
  base_height = 60,
  base_width = 180,
  units = "mm",
  dpi = 300
)
