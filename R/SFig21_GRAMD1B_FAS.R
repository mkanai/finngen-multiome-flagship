library(dplyr)
library(ggplot2)
source(here::here("R/const.R"))

################################################################################
# rs35923643, chr11_123484683_A_G
# rs2147420, chr10_88999856_A_G
# tabix -h gs://finngen-public-data-r12/meta_analysis/mvp_ukbb/summary_stats/C3_CLL_EXALLC_meta_out.tsv.gz 10:88899856-89099856 11:122484683-124484683

start.GRAMD1B <- 123484683 - 5e4
end.GRAMD1B <- 123484683 + 5e4
highlight_pos.GRAMD1B <- c(123484683)

start.FAS <- 88999856 - 5e4
end.FAS <- 88999856 + 5e4
highlight_pos.FAS <- c(88999856)

df.cll.GRAMD1B <-
  rgsutil:::fread_wrapper("data/C3_CLL_EXALLC_meta_out.GRAMD1B.FAS.txt.bgz") %>%
  dplyr::filter(`#CHR` == 11) %>%
  dplyr::mutate(
    chromosome = stringr::str_c("chr", `#CHR`),
    position = POS,
    variant = stringr::str_c("chr", SNP),
    r2 = 0,
    pip = NA,
    cs_id = NA
  ) %>%
  locusviz::preprocess(beta_col = "all_inv_var_meta_beta",
                       se_col = "all_inv_var_meta_sebeta",
                       pvalue_col = "all_inv_var_meta_p") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.GRAMD1B) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.cll.FAS <-
  rgsutil:::fread_wrapper("data/C3_CLL_EXALLC_meta_out.GRAMD1B.FAS.txt.bgz") %>%
  dplyr::filter(`#CHR` == 10) %>%
  dplyr::mutate(
    chromosome = stringr::str_c("chr", `#CHR`),
    position = POS,
    variant = stringr::str_c("chr", SNP),
    r2 = 0,
    pip = NA,
    cs_id = NA
  ) %>%
  locusviz::preprocess(beta_col = "all_inv_var_meta_beta",
                       se_col = "all_inv_var_meta_sebeta",
                       pvalue_col = "all_inv_var_meta_p") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.FAS) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.gex.GRAMD1B = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  l1.cell_types.primary,
  "ENSG00000023171"
)

df.gex.GRAMD1B.B =
  dplyr::filter(df.gex.GRAMD1B, cell_type == "predicted.celltype.l1.B") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.GRAMD1B) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.gex.FAS =
  munge_zfile(
    "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
    c(
      l1.cell_types.primary,
      "predicted.celltype.l2.CD14_Mono",
      "predicted.celltype.l2.CD16_Mono"
    ),
    "ENSG00000026103"
  )

df.gex.FAS.B =
  dplyr::filter(df.gex.FAS, cell_type == "predicted.celltype.l1.B") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.FAS) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

p.gex.GRAMD1B.forest =
  dplyr::filter(df.gex.GRAMD1B, rsid == "chr11_123484683_A_G" &
                  p < 5e-8) %>%
  plot_forest(hide.xtitle = TRUE) +
  theme(plot.title = element_text(hjust = 0.15)) +
  labs(title = "rs35923643")

p.gex.FAS.forest =
  dplyr::filter(df.gex.FAS, rsid == "chr10_88999856_A_G" &
                  p < 5e-8) %>%
  plot_forest(
    colors = c(l1.colors, l2.colors[c("CD14_Mono", "CD16_Mono")]),
    labels = c(l1.labels, l2.labels[c("CD14_Mono", "CD16_Mono")]),
    hide.xtitle = FALSE
  ) +
  theme(plot.title = element_text(hjust = 0.15)) +
  labs(title = "rs2147420", x = "Effect size (ALT)")


p.cll.GRAMD1B =
  locusviz::plot_manhattan_panel(
    df.cll.GRAMD1B,
    highlight_pos = highlight_pos.GRAMD1B,
    xlim = c(start.GRAMD1B, end.GRAMD1B),
    title = "CLL (FinnGen + MVP + UKBB)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = c(1, 1.2)) +
  geom_text(
    aes(
      x = df.cll.GRAMD1B$position[df.cll.GRAMD1B$lead_variant],
      y = df.cll.GRAMD1B$nlog10p[df.cll.GRAMD1B$lead_variant],
      label = "rs35923643"
    ),
    vjust = -1,
    size = 2.5,
    data = NULL
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.38), 0)) +
  theme(plot.margin = margin(t = 4))

p.gex.GRAMD1B =
  locusviz::plot_manhattan_panel(
    df.gex.GRAMD1B.B,
    highlight_pos = highlight_pos.GRAMD1B,
    xlim = c(start.GRAMD1B, end.GRAMD1B),
    title = "GRAMD1B eQTL (B cells)"
  ) +
  locusviz::get_default_theme(hide.xtitle = TRUE, legend.position = "none")

p.gene.GRAMD1B =
  locusviz::plot_gene_panel(
    "chr11",
    start.GRAMD1B,
    end.GRAMD1B,
    genome_build = "hg38",
    highlight_pos = highlight_pos.GRAMD1B
  )

p.cll.FAS =
  locusviz::plot_manhattan_panel(
    df.cll.FAS,
    highlight_pos = highlight_pos.FAS,
    xlim = c(start.FAS, end.FAS),
    title = "CLL (FinnGen + MVP + UKBB)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none") +
  geom_text(
    aes(
      x = df.cll.FAS$position[df.cll.FAS$lead_variant],
      y = df.cll.FAS$nlog10p[df.cll.FAS$lead_variant],
      label = "rs2147420"
    ),
    vjust = -1,
    size = 2.5,
    data = NULL
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.35), 0)) +
  theme(plot.margin = margin(t = 4))

p.gex.FAS =
  locusviz::plot_manhattan_panel(
    df.gex.FAS.B,
    highlight_pos = highlight_pos.FAS,
    xlim = c(start.FAS, end.FAS),
    title = "FAS eQTL (B cells)"
  ) +
  locusviz::get_default_theme(hide.xtitle = TRUE, legend.position = "none")

p.gene.FAS =
  locusviz::plot_gene_panel("chr10",
                            start.FAS,
                            end.FAS,
                            genome_build = "hg38",
                            highlight_pos = highlight_pos.FAS)

layout = "
AG
BG
CG
DG
EG
FG
"

p.GRAMD1B_FAS =
  list(
    p.cll.GRAMD1B + labs(tag = "a"),
    p.gex.GRAMD1B,
    p.gene.GRAMD1B,
    p.cll.FAS  + labs(tag = "c"),
    p.gex.FAS,
    p.gene.FAS,
    patchwork::free((p.gex.GRAMD1B.forest + labs(tag = "b")) / (p.gex.FAS.forest + labs(tag = "d")),
                    side = "l"
    )
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(design = layout,
                         widths = c(3, 1),
                         heights = rep(c(0.35, 0.35, 0.02), 2))

cowplot::save_plot(
  "figures/SFig21_GRAMD1B_FAS.pdf",
  p.GRAMD1B_FAS,
  base_height = 90,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig21_GRAMD1B_FAS.png",
  p.GRAMD1B_FAS,
  base_height = 90,
  base_width = 180,
  units = "mm",
  dpi = 300
)
