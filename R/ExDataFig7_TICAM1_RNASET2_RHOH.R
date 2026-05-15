library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.open4gene = readRDS("data/open4gene.rds")
df.open4gene.sig = readRDS("data/open4gene.sig.rds")

# rs10425559 - chr19_4837475_A_G
# rs10424978 _ chr19_4837545_C_A
start.TICAM1 <- 4837475 - 8e3
end.TICAM1 <- 4837475 + 5e3
highlight_pos.TICAM1 <- c(4837475, 4837545)
peak.TICAM1.start =  4837478
peak.TICAM1.end = 4837977

# rs3798307 - chr6_166950460_G_A
start.RNASET2 <- 166950460 - 3e4
end.RNASET2 <- 166950460 + 1e4
highlight_pos.RNASET2 <- c(166950460)
peak.RNASET2.start <- 166950262
peak.RNASET2.end <- 166950954

# rs13105678 - chr4_40258759_C_A
start.RHOH <- 40258759 - 75000
end.RHOH <- 40258759 + 15000
highlight_pos.RHOH <- c(40258759)
peak.RHOH.start <- 40258141
peak.RHOH.end <- 40259058


# tabix -h gs://finngen-production-library-green/finngen_R12/finngen_R12_analysis_data/finemap/full/susie/E4_HYTHY_AI_STRICT.SUSIE.snp.bgz chr19:3837475-5837475 | cat SUSIE.snp.header.txt - | bgzip -c > E4_HYTHY_AI_STRICT.TICAM1.SUSIE.snp.gz
# tabix -h gs://finngen-production-library-green/finngen_R12/finngen_R12_analysis_data/finemap/full/susie/E4_HYTHY_AI_STRICT.SUSIE.snp.bgz chr6:166850460-167050460 | cat SUSIE.snp.header.txt - | bgzip -c > E4_HYTHY_AI_STRICT.RNASET2.SUSIE.snp.gz
# tabix -h gs://finngen-production-library-green/finngen_R12/finngen_R12_analysis_data/finemap/full/susie/E4_HYTHY_AI_STRICT.SUSIE.snp.bgz chr4:39258759-41258759 | cat SUSIE.snp.header.txt - | bgzip -c > E4_HYTHY_AI_STRICT.RHOH.SUSIE.snp.gz

df.aiht.TICAM1 <-
  rgsutil:::fread_wrapper("data/E4_HYTHY_AI_STRICT.TICAM1.SUSIE.snp.gz") %>%
  dplyr::mutate(r2 = 0) %>%
  locusviz::preprocess(variant_col = "rsid",
                       pip_col = "prob",
                       cs_id_col = "cs") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.TICAM1[1]) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.aiht.RNASET2 <-
  rgsutil:::fread_wrapper("data/E4_HYTHY_AI_STRICT.RNASET2.SUSIE.snp.gz") %>%
  dplyr::mutate(r2 = 0) %>%
  locusviz::preprocess(variant_col = "rsid",
                       pip_col = "prob",
                       cs_id_col = "cs") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.RNASET2) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.aiht.RHOH <-
  rgsutil:::fread_wrapper("data/E4_HYTHY_AI_STRICT.RHOH.SUSIE.snp.gz") %>%
  dplyr::mutate(r2 = 0) %>%
  locusviz::preprocess(variant_col = "rsid",
                       pip_col = "prob",
                       cs_id_col = "cs") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.RHOH) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.gex.TICAM1 = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  c(l1.cell_types.primary, "predicted.celltype.l2.CD16_Mono"),
  "ENSG00000127666"
)

df.gex.TICAM1.CD16_Mono =
  dplyr::filter(df.gex.TICAM1, cell_type == "predicted.celltype.l2.CD16_Mono") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.TICAM1[1]) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.atac.TICAM1 = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/atac_results/cis_nominal/zfiles/%s/integrated_atac_batch1_5.fgid.%s.sum.inv.%s.%s.z",
  setdiff(
    c(l1.cell_types.primary, "predicted.celltype.l2.CD16_Mono"),
    c("predicted.celltype.l1.B")
  ),
  "chr19-4837478-4837977"
)

df.link.TICAM1 =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% c(l1.cell_types.primary, "predicted.celltype.l2.CD16_Mono") &
      gene_id == "ENSG00000127666" &
      .env$start.TICAM1 < peak_end & peak_start < .env$end.TICAM1
  ) %>%
  dplyr::mutate(
    # MANE select TSS
    tss = 4831712,
    peak_mid = (peak_start + peak_end) / 2,
    start = ifelse(tss < peak_mid, tss, peak_mid),
    end = ifelse(tss < peak_mid, peak_mid, tss),
    group = peak_id,
    score = hurdle_count_beta,
    beta =  hurdle_count_beta
  )

df.gex.RNASET2 = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  l1.cell_types.primary,
  "ENSG00000026297"
)

df.gex.RNASET2.Mono =
  dplyr::filter(df.gex.RNASET2, cell_type == "predicted.celltype.l1.Mono") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.RNASET2) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.atac.RNASET2 = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/atac_results/cis_nominal/zfiles/%s/integrated_atac_batch1_5.fgid.%s.sum.inv.%s.%s.z",
  l1.cell_types.primary,
  "chr6-166950262-166950954"
)

df.link.RNASET2 =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% l1.cell_types.primary &
      gene_id == "ENSG00000026297" &
      .env$start.RNASET2 < peak_end & peak_start < .env$end.RNASET2
  ) %>%
  dplyr::mutate(
    peak_mid = (peak_start + peak_end) / 2,
    start = ifelse(tss < peak_mid, tss, peak_mid),
    end = ifelse(tss < peak_mid, peak_mid, tss),
    group = peak_id,
    score = hurdle_count_beta,
    beta =  hurdle_count_beta
  )


df.gex.RHOH = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  l1.cell_types.primary,
  "ENSG00000168421"
)

df.gex.RHOH.CD4_T =
  dplyr::filter(df.gex.RHOH, cell_type == "predicted.celltype.l1.CD4_T") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.RHOH) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.atac.RHOH = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/atac_results/cis_nominal/zfiles/%s/integrated_atac_batch1_5.fgid.%s.sum.inv.%s.%s.z",
  l1.cell_types.primary,
  "chr4-40258141-40259058"
)

df.link.RHOH =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% l1.cell_types.primary &
      gene_id == "ENSG00000168421" &
      .env$start.RHOH < peak_end & peak_start < .env$end.RHOH
  ) %>%
  dplyr::mutate(
    # MANE select
    tss = 40193663,
    peak_mid = (peak_start + peak_end) / 2,
    start = ifelse(tss < peak_mid, tss, peak_mid),
    end = ifelse(tss < peak_mid, peak_mid, tss),
    group = peak_id,
    score = hurdle_count_beta,
    beta =  hurdle_count_beta
  )

################################################################################
p.aiht.TICAM1 =
  locusviz::plot_manhattan_panel(
    df.aiht.TICAM1,
    highlight_pos = highlight_pos.TICAM1,
    background.layers = peak_background(peak.TICAM1.start, peak.TICAM1.end),
    xlim = c(start.TICAM1, end.TICAM1),
    title = "Autoimmune hypothyroidism (FinnGen R12)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = c(1, 1.2)) +
  geom_text(
    aes(
      x = df.aiht.TICAM1$position[df.aiht.TICAM1$position %in% highlight_pos.TICAM1],
      y = df.aiht.TICAM1$nlog10p[df.aiht.TICAM1$position %in% highlight_pos.TICAM1],
      label = c("rs10425559", "rs10424978")
    ),
    hjust = c(1.2, -0.2),
    size = 2.5,
    data = NULL
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.38), 0))

p.gex.TICAM1 =
  locusviz::plot_manhattan_panel(
    df.gex.TICAM1.CD16_Mono,
    highlight_pos = highlight_pos.TICAM1,
    background.layers = peak_background(peak.TICAM1.start, peak.TICAM1.end),
    xlim = c(start.TICAM1, end.TICAM1),
    title = "TICAM1 eQTL (CD16+ Monocytes)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none") +
  scale_y_continuous(expand = expansion(c(0, 0.38), 0))


p.link.TICAM1 =
  dplyr::mutate(df.link.TICAM1,
                highlight = (
                  peak_id == "chr19-4837478-4837977" &
                    gene_id == "ENSG00000127666"
                )) %>%
  plot_links(
    start.TICAM1,
    end.TICAM1,
    highlight_pos.TICAM1,
    # xbreaks = setdiff(seq(52e5, 56e5, by = 1e5), 54e5),
    hide.xtitle = FALSE,
    background.layers = peak_background(peak.TICAM1.start, peak.TICAM1.end, flip_y = TRUE),
  ) +
  labs(x = "Chromosome 19", title = "TICAM1 Link")

p.gene.TICAM1 =
  locusviz::plot_gene_panel(
    "chr19",
    start.TICAM1,
    end.TICAM1,
    genome_build = "hg38",
    highlight_pos = highlight_pos.TICAM1,
    background.layers = peak_background(peak.TICAM1.start, peak.TICAM1.end),
  ) +
  theme(axis.title.x = element_blank())

p.aiht.RNASET2 =
  locusviz::plot_manhattan_panel(
    df.aiht.RNASET2,
    highlight_pos = highlight_pos.RNASET2,
    background.layers = peak_background(peak.RNASET2.start, peak.RNASET2.end),
    xlim = c(start.RNASET2, end.RNASET2),
    title = "Autoimmune hypothyroidism (FinnGen R12)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = c(1, 1.2)) +
  geom_text(
    aes(
      x = df.aiht.RNASET2$position[df.aiht.RNASET2$lead_variant],
      y = df.aiht.RNASET2$nlog10p[df.aiht.RNASET2$lead_variant],
      label = "rs3798307"
    ),
    vjust = -1,
    size = 2.5,
    data = NULL
  ) +
  theme(plot.margin = margin(t = 4))

p.gex.RNASET2 =
  locusviz::plot_manhattan_panel(
    df.gex.RNASET2.Mono,
    highlight_pos = highlight_pos.RNASET2,
    background.layers = peak_background(peak.RNASET2.start, peak.RNASET2.end),
    xlim = c(start.RNASET2, end.RNASET2),
    title = "RNASET2 eQTL (Monocytes)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none")

p.link.RNASET2 =
  dplyr::mutate(
    df.link.RNASET2,
    highlight = (
      peak_id == "chr6-166950262-166950954" &
        gene_id == "ENSG00000026297"
    )
  ) %>%
  plot_links(
    start.RNASET2,
    end.RNASET2,
    highlight_pos.RNASET2,
    # xbreaks = setdiff(seq(52e5, 56e5, by = 1e5), 54e5),
    hide.xtitle = FALSE,
    background.layers = peak_background(peak.RNASET2.start, peak.RNASET2.end, flip_y = TRUE),
  ) +
  labs(x = "Chromosome 6", title = "RNASET2 Link")

p.gene.RNASET2 =
  locusviz::plot_gene_panel(
    "chr6",
    start.RNASET2,
    end.RNASET2,
    genome_build = "hg38",
    highlight_pos = highlight_pos.RNASET2,
    background.layers = peak_background(peak.RNASET2.start, peak.RNASET2.end),
  ) +
  theme(axis.title.x = element_blank())

p.aiht.RHOH =
  locusviz::plot_manhattan_panel(
    df.aiht.RHOH,
    highlight_pos = highlight_pos.RHOH,
    background.layers = peak_background(peak.RHOH.start, peak.RHOH.end),
    xlim = c(start.RHOH, end.RHOH),
    title = "Autoimmune hypothyroidism (FinnGen R12)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = c(1, 1.2)) +
  geom_text(
    aes(
      x = df.aiht.RHOH$position[df.aiht.RHOH$lead_variant],
      y = df.aiht.RHOH$nlog10p[df.aiht.RHOH$lead_variant],
      label = "rs13105678"
    ),
    vjust = -1,
    size = 2.5,
    data = NULL
  ) +
  theme(plot.margin = margin(t = 4))

p.gex.RHOH =
  locusviz::plot_manhattan_panel(
    df.gex.RHOH.CD4_T,
    highlight_pos = highlight_pos.RHOH,
    background.layers = peak_background(peak.RHOH.start, peak.RHOH.end),
    xlim = c(start.RHOH, end.RHOH),
    title = "RHOH eQTL (CD4+ T)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none")

p.link.RHOH =
  dplyr::mutate(
    df.link.RHOH,
    highlight = (
      peak_id == "chr4-40258141-40259058" &
        gene_id == "ENSG00000168421"
    )
  ) %>%
  plot_links(
    start.RHOH,
    end.RHOH,
    highlight_pos.RHOH,
    # xbreaks = setdiff(seq(52e5, 56e5, by = 1e5), 54e5),
    hide.xtitle = FALSE,
    background.layers = peak_background(peak.RHOH.start, peak.RHOH.end, flip_y = TRUE),
  ) +
  labs(x = "Chromosome 4", title = "RHOH Link")

p.gene.RHOH =
  locusviz::plot_gene_panel(
    "chr4",
    start.RHOH,
    end.RHOH,
    genome_build = "hg38",
    highlight_pos = highlight_pos.RHOH,
    background.layers = peak_background(peak.RHOH.start, peak.RHOH.end),
  ) +
  theme(axis.title.x = element_blank())

p.TICAM1_RNASET2_RHOH =
  list(
    p.aiht.TICAM1,
    p.gex.TICAM1,
    p.gene.TICAM1,
    p.link.TICAM1,
    p.aiht.RNASET2,
    p.gex.RNASET2,
    p.gene.RNASET2,
    p.link.RNASET2,
    p.aiht.RHOH,
    p.gex.RHOH,
    p.gene.RHOH,
    p.link.RHOH
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1, heights = rep(c(rep(1, 2), 0.2, 0.7), 3))

cowplot::save_plot(
  "figures/ExtendedDataFig7_TICAM1_locuszoom.pdf",
  p.TICAM1_RNASET2_RHOH,
  base_height = 170,
  base_width = 180,
  units = "mm"
)

################################################################################

# rs10425559 - chr19_4837475_A_G
# rs10424978 _ chr19_4837545_C_A
p.gex.TICAM1.forest =
  dplyr::filter(df.gex.TICAM1, rsid == "chr19_4837475_A_G") %>%
  plot_forest(
    colors = c(l1.colors, l2.colors["CD16_Mono"]),
    labels = c(l1.labels, l2.labels["CD16_Mono"]),
    legend.position = c(1, 0),
    legend.justification = c(1, 0)
  ) +
  labs(title = "TICAM1 eQTL (rs10425559)")

p.atac.TICAM1.forest =
  dplyr::filter(df.atac.TICAM1, rsid == "chr19_4837475_A_G") %>%
  plot_forest(
    colors = c(l1.colors, l2.colors["CD16_Mono"]),
    labels = c(l1.labels, l2.labels["CD16_Mono"]),
    legend.position = c(1, 0),
    legend.justification = c(1, 0)
  ) +
  labs(title = "caQTL (rs10425559)")

p.link.TICAM1.forest =
  dplyr::filter(df.open4gene,
                peak_id == "chr19-4837478-4837977" &
                  gene_id == "ENSG00000127666") %>%
  dplyr::rename(beta = hurdle_count_beta, se = hurdle_count_se) %>%
  dplyr::filter(cell_type %in% c(l1.cell_types.primary, "predicted.celltype.l2.CD16_Mono")) %>%
  dplyr::mutate(sig = ifelse(sig, "FDR < 0.05", "FDR >= 0.05")) %>%
  plot_forest(
    colors = c(l1.colors, l2.colors["CD16_Mono"]),
    labels = c(l1.labels, l2.labels["CD16_Mono"]),
    legend.position = c(1, 0),
    legend.justification = c(1, 0)
  ) +
  labs(title = "Link")

# rs3798307 - chr6_166950460_G_A
p.gex.RNASET2.forest =
  dplyr::filter(df.gex.RNASET2, rsid == "chr6_166950460_G_A") %>%
  plot_forest() +
  labs(title = "RNASET2 eQTL (rs3798307)")

p.atac.RNASET2.forest =
  dplyr::filter(df.atac.RNASET2, rsid == "chr6_166950460_G_A") %>%
  plot_forest() +
  labs(title = "caQTL (rs3798307)")

p.link.RNASET2.forest =
  dplyr::filter(df.open4gene,
                peak_id == "chr6-166950262-166950954" &
                  gene_id == "ENSG00000026297") %>%
  dplyr::rename(beta = hurdle_count_beta, se = hurdle_count_se) %>%
  dplyr::filter(cell_type %in% l1.cell_types.primary) %>%
  dplyr::mutate(sig = ifelse(sig, "FDR < 0.05", "FDR >= 0.05")) %>%
  plot_forest() +
  labs(title = "Link")

# rs13105678 - chr4_40258759_C_A
p.gex.RHOH.forest =
  dplyr::filter(df.gex.RHOH, rsid == "chr4_40258759_C_A") %>%
  plot_forest() +
  labs(title = "RHOH eQTL (rs13105678)")

p.atac.RHOH.forest =
  dplyr::filter(df.atac.RHOH, rsid == "chr4_40258759_C_A") %>%
  plot_forest() +
  labs(title = "caQTL (rs13105678)")

p.link.RHOH.forest =
  dplyr::filter(df.open4gene,
                peak_id == "chr4-40258141-40259058" &
                  gene_id == "ENSG00000168421") %>%
  dplyr::rename(beta = hurdle_count_beta, se = hurdle_count_se) %>%
  dplyr::filter(cell_type %in% l1.cell_types.primary) %>%
  dplyr::mutate(sig = ifelse(sig, "FDR < 0.05", "FDR >= 0.05")) %>%
  plot_forest() +
  labs(title = "Link")


p.forest = list(
  p.gex.TICAM1.forest,
  p.atac.TICAM1.forest,
  p.link.TICAM1.forest,
  p.gex.RNASET2.forest,
  p.atac.RNASET2.forest,
  p.link.RNASET2.forest,
  p.gex.RHOH.forest,
  p.atac.RHOH.forest,
  p.link.RHOH.forest
) %>%
  purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a")
p.forest

cowplot::save_plot(
  "figures/SFig12_TICAM1_forest.pdf",
  p.forest,
  base_height = 170,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig12_TICAM1_forest.png",
  p.forest,
  base_height = 170,
  base_width = 180,
  units = "mm",
  dpi = 300
)
