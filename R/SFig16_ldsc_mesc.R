library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.mesc = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/mesc/results_indiv/*.mesc.all.h2med.gz", function(df, path) {
  dplyr::mutate(
    df,
    Z = Estimate / `SE(Estimate)`,
    QTL = dplyr::case_when(
      basename(path) == "integrated_atac_batch1_5.fgid.mesc.all.h2med.gz" ~ "caQTL",
      basename(path) == "integrated_gex_batch1_5.fgid.qc.mesc.all.h2med.gz" ~ "eQTL",
      basename(path) == "integrated_gex_batch1_5.fgid.qc.on_atac.mesc.all.h2med.gz" ~ "eQTL on caQTL",
      basename(path) == "integrated_atac_batch1_5.fgid.on_gex.mesc.all.h2med.gz" ~ "caQTL on eQTL",
      basename(path) == "Olink_pQTL_2023_10_11.mesc.all.h2med.gz" ~ "pQTL",
      basename(path) == "caQTL_eQTL.mesc.all.h2med.gz" ~ "caQTL & eQTL",
      basename(path) == "caQTL_pQTL.mesc.all.h2med.gz" ~ "caQTL & pQTL",
      basename(path) == "eQTL_pQTL.mesc.all.h2med.gz" ~ 'eQTL & pQTL',
      basename(path) == "caQTL_eQTL_pQTL.mesc.all.h2med.gz" ~ "caQTL & eQTL & pQTL",
      TRUE ~ basename(path)
    ),
    modality = factor(
      QTL,
      levels = c(
        "caQTL",
        "eQTL",
        "pQTL",
        "caQTL & eQTL",
        "caQTL & pQTL",
        "eQTL & pQTL",
        "caQTL & eQTL & pQTL"
      )
    )
  )
}) %>%
  dplyr::filter(!QTL %in% c("eQTL on caQTL", "caQTL on eQTL"))

meta_reml <- function(df, group_vars) {
  if (nrow(df) > 1) {
    fit <- metafor::rma(yi = df$beta,
                        sei = df$se,
                        method = "REML")
    result <- tibble::tibble(beta_meta = coef(fit), se_meta = fit$se)
  } else {
    result <- tibble::tibble(beta_meta = df$beta[1], se_meta = df$se[1])
  }
  
  group_cols <- df %>%
    dplyr::select(all_of(group_vars)) %>%
    dplyr::slice(1)
  
  dplyr::bind_cols(group_cols, result)
}

group_vars.mesc <- c("Category", "Cell_Type", "modality")
df.mesc.meta <-
  dplyr::filter(df.mesc,
                Quantity == "h2med" &
                  Cell_Type %in% c("predicted.celltype.l1.PBMC", "bulk")) %>%
  dplyr::mutate(
    Cell_Type = remove_cell_type_prefix(Cell_Type),
    Category = dplyr::case_when(
      Phenotype %in% immune_traits ~ "Immune",
      Phenotype %in% cardiometabolic_traits ~ "Cardiometabolic",
      Phenotype %in% cancer_traits ~ "Cancer",
      TRUE ~ NA_character_
    ),
    Category = factor(Category, levels = c("Immune", "Cancer", "Cardiometabolic")),
    beta = Estimate_over_h2,
    se = `SE(Estimate_over_h2)`
  ) %>%
  dplyr::filter(!is.na(Category)) %>%
  dplyr::group_split(!!!syms(group_vars.mesc), .keep = TRUE) %>%
  purrr::map_dfr(~ meta_reml(.x, group_vars.mesc))

if (!file.exists("tables/ST30_mesc_meta.tsv")) {
  dplyr::rename(df.mesc.meta,
                Estimate_over_h2 = beta_meta,
                SE_Estimate_over_h2 = se_meta) %>%
    export_table("tables/ST30_mesc_meta.tsv", "ST30")
}

plot_mesc = function(df) {
  ggplot(df, aes(modality, beta_meta, fill = modality)) +
    geom_hline(
      aes(yintercept = beta_meta),
      data = dplyr::filter(df, modality == "caQTL"),
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_col() +
    geom_errorbar(
      aes(ymin = beta_meta - se_meta, ymax = beta_meta + se_meta),
      width = 0,
      color = "grey20"
    ) +
    scale_fill_manual(values = qtl.colors) +
    locusviz::get_default_theme(legend.position = "none", hide.xlab = TRUE) +
    theme(strip.background = element_blank()) +
    scale_y_continuous(expand = expansion(),
                       labels = scales::label_number(drop0trailing = TRUE))
}

p.mesc <-
  plot_mesc(df.mesc.meta) +
  facet_grid(~ Category) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = expression("Meta-analyzed " * italic(h)[med]^2 / italic(h)[g]^2))
p.mesc

p.mesc.trait <-
  dplyr::filter(
    df.mesc,
    Quantity == "h2med" &
      Cell_Type %in% c("predicted.celltype.l1.PBMC", "bulk") &
      Phenotype %in% c("AUTOIMMUNE", "C3_SKIN_EXALLC", "T2D")
  ) %>%
  dplyr::mutate(
    beta_meta = Estimate_over_h2,
    se_meta = `SE(Estimate_over_h2)`,
    Phenotype = factor(
      Phenotype,
      levels = c("AUTOIMMUNE", "C3_SKIN_EXALLC", "T2D"),
      labels = c("Autoimmune", "Skin Cancer", "T2D")
    )
  ) %>%
  plot_mesc() +
  facet_grid(~ Phenotype) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = expression("Meta-analyzed " * italic(h)[med]^2 / italic(h)[g]^2))
p.mesc.trait

modalities <- c("caQTL", "eQTL", "pQTL")

combos <- purrr::map(seq_along(modalities), function(i) {
  combn(modalities, i, simplify = FALSE)
}) %>%
  purrr::flatten()

mesc.upset <- purrr::map_dfr(seq_len(length(combos) * 3), function(i) {
  data.frame(
    i = (i - 1) %% length(combos) + 1,
    facet = (i - 1) %/% length(combos),
    modality = modalities,
    value = as.integer(modalities %in% combos[[(i - 1) %% length(combos) + 1]])
  )
}) %>%
  dplyr::mutate(value = factor(value),
                modality = factor(modality, levels = rev(modalities)))

mesc.upset.segments <-
  dplyr::filter(mesc.upset, value == 1) %>%
  dplyr::group_by(i, facet) %>%
  dplyr::summarize(
    ymin = min(as.numeric(modality)),
    ymax = max(as.numeric(modality)),
    n = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n > 1)

p.mesc.upset =
  ggplot() +
  geom_point(aes(factor(i), modality, color = value), data = mesc.upset) +
  geom_segment(aes(
    x = i,
    xend = i,
    y = ymin,
    yend = ymax
  ), data = mesc.upset.segments) +
  scale_color_manual(values = c("0" = "grey80", "1" = "black")) +
  # scale_y_discrete(expand = expansion()) +
  locusviz::get_default_theme(
    legend.position = "none",
    hide.ytitle = TRUE,
    hide.xlab = TRUE
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  facet_grid(~ facet)

list((p.mesc + theme(axis.text.x = element_blank())), p.mesc.upset) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1, heights = c(1, 0.4))


################################################################################
if (!file.exists("data/df.ldsc.rds")) {
  df.ldsc = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/ldsc/*.all.results.gz", function(df, path) {
    dplyr::filter(df, Category == "L2_0") %>%
      dplyr::mutate(annot = stringr::str_split_fixed(basename(path), "\\.", 2)[, 1])
  }) %>%
    dplyr::mutate(
      annot = forcats::fct_recode(
        annot,
        "All peaks" = "FinnGen_CA",
        "CA-Link+" = "FinnGen_CA_Link",
        "CA-Link-" = "FinnGen_CA_NoLink",
        "Switch" = "FinnGen_CA_Switch_Link",
        "Rheostat" = "FinnGen_CA_Rheostat_Link",
        "Dual" = "FinnGen_CA_Dual_Link",
        "caPeak" = "FinnGen_CA_sig"
      ) %>%
        forcats::fct_relevel(
          "All peaks",
          "CA-Link+",
          "CA-Link-",
          "Switch",
          "Rheostat",
          "Dual",
          "caPeak"
        )
    )
  saveRDS(df.ldsc, "data/df.ldsc.rds")
}

df.ldsc = readRDS("data/df.ldsc.rds") %>%
  dplyr::filter(annot != "caPeak")

group_vars.ldsc <- c("Category", "Cell_Type", "annot")
df.ldsc.category <-
  dplyr::filter(df.ldsc, Cell_Type == "predicted.celltype.l1.PBMC") %>%
  dplyr::mutate(
    Cell_Type = remove_cell_type_prefix(Cell_Type),
    Category = dplyr::case_when(
      Phenotype %in% immune_traits ~ "Immune",
      Phenotype %in% cardiometabolic_traits ~ "Cardiometabolic",
      Phenotype %in% cancer_traits ~ "Cancer",
      TRUE ~ NA_character_
    ),
    Category = factor(Category, levels = c("Immune", "Cancer", "Cardiometabolic"))
  ) %>%
  dplyr::filter(!is.na(Category))

df.ldsc.h2.meta =
  dplyr::mutate(df.ldsc.category, beta = Prop._h2, se = Prop._h2_std_error) %>%
  dplyr::group_split(!!!syms(group_vars.ldsc), .keep = TRUE) %>%
  purrr::map_dfr(~ meta_reml(.x, group_vars.ldsc))


df.ldsc.enrichment.meta =
  dplyr::mutate(df.ldsc.category, beta = Enrichment, se = Enrichment_std_error) %>%
  dplyr::group_split(!!!syms(group_vars.ldsc), .keep = TRUE) %>%
  purrr::map_dfr(~ meta_reml(.x, group_vars.ldsc))

if (!file.exists("tables/ST31_ldsc_meta.tsv")) {
  dplyr::left_join(
    dplyr::rename(
      df.ldsc.h2.meta,
      Prop_h2 = beta_meta,
      Prop_h2_std_error = se_meta
    ),
    dplyr::rename(
      df.ldsc.enrichment.meta,
      Enrichment = beta_meta,
      Enrichment_std_error = se_meta
    ),
    by = c("Category", "Cell_Type", "annot")
  ) %>%
    export_table("tables/ST31_ldsc_meta.tsv", "ST31")
}


plot_ldsc = function(df) {
  pd = position_dodge(width = 0.9)
  ggplot(df, aes(Category, beta_meta, fill = annot)) +
    geom_col(position = pd) +
    geom_errorbar(
      aes(ymin = beta_meta - se_meta, ymax = beta_meta + se_meta),
      width = 0,
      position = pd
    ) +
    locusviz::get_default_theme(hide.xtitle = TRUE) +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values = c(
      link_colors,
      link_mode_colors,
      "All peaks" = BuenColors::jdb_palette("brewer_blue")[9]
    ))
}

p.ldsc.h2.meta =
  dplyr::filter(df.ldsc.h2.meta, !annot %in% c("Switch", "Rheostat", "Dual")) %>%
  plot_ldsc() +
  scale_y_continuous(expand = expansion(c(0, 0.1)), labels = scales::label_percent()) +
  labs(y = expression(paste("Proportion of ", italic(h)[g]^2)))

p.ldsc.h2.meta.mode =
  dplyr::filter(df.ldsc.h2.meta, annot %in% c("Switch", "Rheostat", "Dual")) %>%
  plot_ldsc() +
  scale_y_continuous(expand = expansion(c(0, 0.1)), labels = scales::label_percent()) +
  labs(y = expression(paste("Proportion of ", italic(h)[g]^2)))

p.ldsc.enrichment.meta =
  dplyr::filter(df.ldsc.enrichment.meta,
                !annot %in% c("Switch", "Rheostat", "Dual")) %>%
  plot_ldsc() +
  scale_y_continuous(expand = expansion(c(0, 0.1))) +
  labs(y = expression(paste("Enrichment of ", italic(h)[g]^2))) +
  theme(legend.position = "none")

p.ldsc.enrichment.meta.mode =
  dplyr::filter(df.ldsc.enrichment.meta,
                annot %in% c("Switch", "Rheostat", "Dual")) %>%
  plot_ldsc() +
  scale_y_continuous(expand = expansion(c(0, 0.1))) +
  labs(y = expression(paste("Enrichment of ", italic(h)[g]^2))) +
  theme(legend.position = "none")

p.ldsc.h2.trait =
  dplyr::filter(
    df.ldsc,
    Phenotype %in% c("AUTOIMMUNE", "C3_SKIN_EXALLC", "T2D", "HEIGHT_IRN") &
      Cell_Type == "predicted.celltype.l1.PBMC"
  ) %>%
  dplyr::mutate(Category = Phenotype,
                beta_meta = Prop._h2,
                se_meta = Prop._h2_std_error) %>%
  plot_ldsc() +
  scale_y_continuous(expand = expansion(c(0, 0.1)), labels = scales::label_percent()) +
  labs(y = expression(paste("Proportion of ", italic(h)[g]^2)))

p.ldsc.enrichment.trait =
  dplyr::filter(
    df.ldsc,
    Phenotype %in% c("AUTOIMMUNE", "C3_SKIN_EXALLC", "T2D", "HEIGHT_IRN") &
      Cell_Type == "predicted.celltype.l1.PBMC"
  ) %>%
  dplyr::mutate(
    beta_meta = Enrichment,
    se_meta = Enrichment_std_error,
    Category = factor(
      Phenotype,
      levels = c("AUTOIMMUNE", "C3_SKIN_EXALLC", "T2D", "HEIGHT_IRN"),
      labels = c("Autoimmune", "Skin Cancer", "T2D", "Height")
    )
  ) %>%
  plot_ldsc() +
  scale_y_continuous(expand = expansion(c(0, 0.1))) +
  labs(y = expression(paste("Enrichment of ", italic(h)[g]^2)))
p.ldsc.enrichment.trait

list((p.mesc + theme(axis.text.x = element_blank())), p.mesc.upset) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1, heights = c(1, 0.4))

layout = "
AC
BC
"

p.mesc_ldsc =
  list(
    (p.mesc + labs(tag = "a") + theme(axis.text.x = element_blank())),
    patchwork::free(p.mesc.upset, side = "b"),
    patchwork::free((p.ldsc.h2.meta + labs(tag = "b")) + (p.ldsc.enrichment.meta + labs(tag = "c")),
                    side = "t"
    )
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(design = layout,
                         widths = c(0.6, 1),
                         height = c(1, 0.2))

p.mesc.combined =
  list((p.mesc + theme(axis.text.x = element_blank())), p.mesc.upset) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1, height = c(1, 0.2))

p.ldsc = list(
  p.ldsc.h2.meta,
  p.ldsc.enrichment.meta,
  p.ldsc.h2.meta.mode,
  p.ldsc.enrichment.meta.mode
) %>%
  purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig16_mesc_upset.pdf",
  p.mesc.combined,
  base_height = 90,
  base_width = 90,
  units = "mm"
)

cowplot::save_plot(
  "figures/SFig16_mesc_upset.png",
  p.mesc.combined,
  base_height = 90,
  base_width = 90,
  units = "mm",
  dpi = 300
)

cowplot::save_plot(
  "figures/SFig17_ldsc.pdf",
  p.ldsc,
  base_height = 120,
  base_width = 120,
  units = "mm"
)

cowplot::save_plot(
  "figures/SFig17_ldsc.png",
  p.ldsc,
  base_height = 120,
  base_width = 120,
  units = "mm",
  dpi = 300
)
