library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.open4gene = readRDS("data/open4gene.rds")
df.open4gene.sig = readRDS("data/open4gene.sig.rds")

df.ibd.TNRC18 <-
  rgsutil:::fread_wrapper(
    "data/K11_IBD_STRICT.SUSIE.TNRC18.snp.bgz"
  ) %>%
  dplyr::mutate(r2 = 0) %>%
  locusviz::preprocess(variant_col = "rsid",
                       pip_col = "prob",
                       cs_id_col = "cs") %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

immune_traits.TRNC18 = c(
  "K11_IBD_STRICT" = "IBD",
  "M13_ANKYLOSPON" = "AS",
  "L12_PSORIASIS" = "PS",
  "T1D" = "T1D",
  "E4_HYTHY_AI_STRICT" = "AIHT",
  "E4_GRAVES_STRICT" = "GD"
)

start.TNRC18 <- 5397122 - 2e5
end.TNRC18 <- 5397122 + 2e5
highlight_pos.TNRC18 <- c(5397122)
peak.TNRC18.start = 5397119
peak.TNRC18.end = 5397618


df.ibd.forest = data.table::fread(
  "data/7_5397122_C_T_phenotype_associations.tsv",
  data.table = F
) %>%
  dplyr::mutate(
    sebeta = beta / qnorm(pval, lower.tail = F),
    trait = immune_traits.TRNC18[phenocode],
    trait = factor(trait, levels = trait[order(beta)]),
    category = dplyr::case_when(
      trait == "T1D" ~ "IV Endocrine, nutritional and metabolic diseases (E4_)",
      TRUE ~ category
    )
  ) %>%
  dplyr::filter(phenocode %in% names(immune_traits.TRNC18))

p.ibd.forest =
  ggplot(df.ibd.forest, aes(beta, trait)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.3
  ) +
  geom_errorbar(
    aes(xmin = beta - sebeta, xmax = beta + sebeta),
    width = 0,
    linewidth = 0.3,
    color = "grey50"
  ) +
  geom_point() +
  locusviz::get_default_theme(
    legend.position = "none",
    hide.xtitle = TRUE,
    hide.ytitle = TRUE
  ) +
  labs(title = "GWAS")

df.gex.TNRC18 = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  l1.cell_types.primary,
  "ENSG00000182095"
)

df.gex.TNRC18.CD4_T =
  dplyr::filter(df.gex.TNRC18, cell_type == "predicted.celltype.l1.CD4_T") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.atac.TNRC18 = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/atac_results/cis_nominal/zfiles/%s/integrated_atac_batch1_5.fgid.%s.sum.inv.%s.%s.z",
  l1.cell_types.primary,
  "chr7-5397119-5397618"
)

p.gex.TNRC18.forest =
  dplyr::filter(df.gex.TNRC18, rsid == "chr7_5397122_C_T") %>%
  plot_forest(hide.xtitle = TRUE, legend.position = c(0.9, 0), legend.justification = c(0.9, 0)) +
  labs(title = "eQTL")


p.atac.TNRC18.forest =
  dplyr::filter(df.atac.TNRC18, rsid == "chr7_5397122_C_T") %>%
  plot_forest() +
  labs(title = "caQTL", x = "Effect size (ALT)")

df.atac.TNRC18.tracks = import_peak_track("data/fragment.%s.chr7_5397122_C_T.%d.bw", cell_types = l1.cell_types.primary)

df.link.TNRC18 =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% l1.cell_types.primary &
      gene_id == "ENSG00000182095" &
      .env$start.TNRC18 < peak_end & peak_start < .env$end.TNRC18
  ) %>%
  dplyr::mutate(
    peak_mid = (peak_start + peak_end) / 2,
    start = ifelse(tss < peak_mid, tss, peak_mid),
    end = ifelse(tss < peak_mid, peak_mid, tss),
    group = peak_id,
    score = hurdle_count_beta,
    beta =  hurdle_count_beta
  )

p.link.TNRC18.forest =
  dplyr::filter(df.open4gene,
                peak_id == "chr7-5397119-5397618" &
                  gene_id == "ENSG00000182095") %>%
  dplyr::rename(beta = hurdle_count_beta, se = hurdle_count_se) %>%
  dplyr::filter(cell_type %in% l1.cell_types.primary) %>%
  dplyr::mutate(sig = ifelse(sig, "FDR < 0.05", "FDR >= 0.05")) %>%
  plot_forest(hide.xtitle = FALSE, legend.position = c(0, 0.8), legend.justification = c(0, 0.8)) +
  labs(title = "Link")


p.ibd.TNRC18 =
  locusviz::plot_manhattan_panel(
    df.ibd.TNRC18,
    highlight_pos = highlight_pos.TNRC18,
    background.layers = peak_background(peak.TNRC18.start, peak.TNRC18.end),
    xlim = c(start.TNRC18, end.TNRC18),
    title = "IBD (FinnGen R12)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = c(1, 1.2)) +
  geom_text(
    aes(
      x = df.ibd.TNRC18$position[df.ibd.TNRC18$lead_variant],
      y = df.ibd.TNRC18$nlog10p[df.ibd.TNRC18$lead_variant],
      label = "rs748670681"
    ),
    hjust = 1.2,
    size = 2.5,
    data = NULL
  ) + geom_text(
    aes(
      x = df.ibd.TNRC18$position[df.ibd.TNRC18$lead_variant],
      y = df.ibd.TNRC18$nlog10p[df.ibd.TNRC18$lead_variant],
      label = "(PIP = 1.0)"
    ),
    hjust = -0.2,
    size = 2.5,
    data = NULL
  )

p.gex.TNRC18 =
  locusviz::plot_manhattan_panel(
    df.gex.TNRC18.CD4_T,
    highlight_pos = highlight_pos.TNRC18,
    background.layers = peak_background(peak.TNRC18.start, peak.TNRC18.end),
    xlim = c(start.TNRC18, end.TNRC18),
    title = "TNRC18 eQTL (CD4+ T)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none")

p.link.TNRC18 =
  dplyr::mutate(
    df.link.TNRC18,
    highlight = (
      peak_id == "chr7-5397119-5397618" &
        gene_id == "ENSG00000182095"
    )
  ) %>%
  plot_links(
    start.TNRC18,
    end.TNRC18,
    highlight_pos.TNRC18,
    xbreaks = setdiff(seq(52e5, 56e5, by = 1e5), 54e5),
    hide.xtitle = FALSE,
    background.layers = peak_background(peak.TNRC18.start, peak.TNRC18.end, flip_y = TRUE)
  ) +
  labs(x = "Chromosome 7", title = "TNRC18 Link")
p.link.TNRC18

p.gene.TNRC18 =
  locusviz::plot_gene_panel(
    "chr7",
    start.TNRC18,
    end.TNRC18,
    genome_build = "hg38",
    highlight_pos = highlight_pos.TNRC18,
    highlight_pos_y = 2,
    background.layers = peak_background(peak.TNRC18.start, peak.TNRC18.end)
  ) +
  theme(axis.title.x = element_blank(), plot.tag = element_text(face = "bold"))

plts.TNRC18 = list(
  A = p.ibd.TNRC18,
  B = p.gex.TNRC18,
  C = p.gene.TNRC18,
  D = p.link.TNRC18,
  E = p.ibd.forest,
  F = p.gex.TNRC18.forest,
  G = p.link.TNRC18.forest
)

################################################################################
# IL21R

df.asthma.IL21R <-
  rgsutil:::fread_wrapper("data/J10_ASTHMA_EXMORE.IL21R.SUSIE.snp.bgz") %>%
  dplyr::mutate(r2 = 0) %>%
  locusviz::preprocess(variant_col = "rsid",
                       pip_col = "prob",
                       cs_id_col = "cs") %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.aiht.IL21R <-
  rgsutil:::fread_wrapper("data/E4_HYTHY_AI_STRICTE.IL21R.SUSIE.snp.bgz") %>%
  dplyr::mutate(r2 = 0) %>%
  locusviz::preprocess(variant_col = "rsid",
                       pip_col = "prob",
                       cs_id_col = "cs") %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)


if (!file.exists("data/df.bbj.IL21R.rds")) {
  df.bbj.IL21R <- rgsutil:::fread_wrapper("data/ATC_H03A.IL21R.txt") %>%
    dplyr::mutate(
      variant = locusviz::variant_str(chrom, pos, ref, alt),
      locusviz::liftover_variant(stringr::str_replace_all(variant, "_", ":"), genome_build = "hg19")
    ) %>%
    dplyr::select(-variant) %>%
    dplyr::mutate(pip = NA, cs_id = NA, r2 = 0) %>%
    locusviz::preprocess(
      chromosome_col = "new_chromosome",
      position_col = "new_position",
      variant_col = "new_variant",
      se_col = "sebeta",
      pvalue_col = "pval"
    ) %>%
    locusviz::annotate_r2(reference_panel = "1000G", window = 2e6, population = "EAS")
  saveRDS(df.bbj.IL21R, "data/df.bbj.IL21R.rds")
}
df.bbj.IL21R = readRDS("data/df.bbj.IL21R.rds")

df.IL21R.r2 = 
  purrr::map_dfr(seq(3), function(i) {
  data =tibble::tibble(
    variant = c("chr16_27344903_G_A", "chr16_27384341_C_CT", "chr16_27386721_T_C"),
    lead_variant = FALSE
  )
  data$lead_variant[i] = TRUE
  data = locusviz::annotate_r2(data, reference_panel = "sisu42", window = 100000) %>%
    dplyr::mutate(variant2 = variant[lead_variant])
  return(data)
}) %>%
  dplyr::mutate(r2 = ifelse(is.na(r2), 1, r2))

rsid.map = c("chr16_27344903_G_A" = "rs144651842", "chr16_27384341_C_CT" = "rs201121732", "chr16_27386721_T_C" = "rs76715626")

p.IL21R.r2 = 
  dplyr::mutate(
  df.IL21R.r2,
  rsid = rsid.map[variant],
  rsid2 = rsid.map[variant2]
) %>%
ggplot(aes(rsid, rsid2, fill = r2)) +
  geom_tile() +
  geom_text(aes(label = scales::number(r2, accuracy = 0.01)), size = 2) +
  # scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
  scale_fill_gradientn(colors = BuenColors::jdb_palette("brewer_red"), limits = c(0, 1)) +
  scale_x_discrete(expand = expansion()) +
  scale_y_discrete(expand = expansion()) +
  locusviz::get_default_theme(hide.xtitle = TRUE, hide.ytitle = TRUE, legend.position = "none") +
  theme(axis.line = element_blank(),
        plot.title = element_text(margin = margin())) +
  labs(title = "r2")

df.gex.IL21R = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  l1.cell_types.primary,
  "ENSG00000103522"
)

df.gex.IL21R.NK =
  dplyr::filter(df.gex.IL21R, cell_type == "predicted.celltype.l1.NK") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.variant = readRDS("data/df.variant.rds")
cascade.peaks =
  dplyr::filter(df.variant, variant_id == "chr16_27384341_C_CT") %>%
  dplyr::pull(cascade_peak_genes) %>%
  stringr::str_split(",") %>%
  unlist() %>%
  stringr::str_remove("->ENSG00000103522")

p.gex.IL21R.forest =
  dplyr::filter(df.gex.IL21R, rsid == "chr16_27384341_C_CT") %>%
  plot_forest(hide.xtitle = FALSE, legend.position = c(1, 0), legend.justification = c(1, 0)) +
  labs(title = "eQTL")

start.IL21R <- 27270341
end.IL21R <- 27464341
highlight_pos.IL21R <- c(27384341, 27344903, 27386721)
peak.IL21R.start = c(27374916, 27377587, 27383475, 27407901, 27426358, 27463601)
peak.IL21R.end = c(27375415, 27378086, 27384329, 27408400, 27426857, 27464100)
peaks.IL21R = c("chr16-27374916-27375415", "chr16-27377587-27378086", "chr16-27383475-27384329", "chr16-27407901-27408400", "chr16-27426358-27426857", "chr16-27463601-27464100")

df.link.IL21R =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% l1.cell_types.primary &
      gene_id == "ENSG00000103522" &
      .env$start.IL21R < peak_end & peak_start < .env$end.IL21R
  ) %>%
  dplyr::mutate(
    peak_mid = (peak_start + peak_end) / 2,
    start = ifelse(tss < peak_mid, tss, peak_mid),
    end = ifelse(tss < peak_mid, peak_mid, tss),
    group = peak_id,
    score = hurdle_count_beta,
    beta =  hurdle_count_beta
  )

p.asthma.IL21R =
  locusviz::plot_manhattan_panel(
    df.asthma.IL21R,
    highlight_pos = highlight_pos.IL21R,
    background.layers = peak_background(peak.IL21R.start, peak.IL21R.end),
    xlim = c(start.IL21R, end.IL21R),
    title = "Asthma (FinnGen R12)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = c(1, 1.2)) +
  geom_text(
    aes(
      x = 27344903,
      y = df.asthma.IL21R$nlog10p[df.asthma.IL21R$position == 27344903],
      label = c("rs144651842", "(IL4R:p.Ala82Thr)")
    ),
    hjust = 1.2,
    vjust = c(0.5, 2),
    size = 2.5,
    data = NULL
  ) + geom_text(
    aes(
      x = 27344903,
      y = df.asthma.IL21R$nlog10p[df.asthma.IL21R$position == 27344903],
      label = "(PIP = 0.44)"
    ),
    hjust = -0.2,
    size = 2.5,
    data = NULL
  ) + scale_y_continuous(expand = expansion(c(0, 0.3))) +
  theme(plot.margin = margin(t = 4))

p.aiht.IL21R =
  locusviz::plot_manhattan_panel(
    df.aiht.IL21R,
    highlight_pos = highlight_pos.IL21R,
    background.layers = peak_background(peak.IL21R.start, peak.IL21R.end),
    xlim = c(start.IL21R, end.IL21R),
    title = "AIHT (FinnGen R12)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none") +
  geom_text(
    aes(
      x = c(27384341, 27386721),
      y = df.aiht.IL21R$nlog10p[df.aiht.IL21R$position %in% c(27384341, 27386721)],
      label = c("rs201121732", "rs76715626")
    ),
    hjust = c(1.2, -0.2),
    size = 2.5,
    data = NULL
  ) + geom_text(
    aes(
      x = c(27384341, 27386721),
      y = df.aiht.IL21R$nlog10p[df.aiht.IL21R$position %in% c(27384341, 27386721)],
      label = c("(PIP = 0.80)", "(PIP = 0.20)")
    ),
    hjust = c(1.2, -0.2),
    vjust = 2,
    size = 2.5,
    data = NULL
  )

p.bbj.IL21R =
  locusviz::plot_manhattan_panel(
    df.bbj.IL21R,
    highlight_pos = highlight_pos.IL21R,
    background.layers = peak_background(peak.IL21R.start, peak.IL21R.end),
    xlim = c(start.IL21R, end.IL21R),
    title = "Thyroid preparations (BBJ)"
  ) + 
  locusviz::get_default_theme(hide.xlab = FALSE, legend.position = "none") +
  geom_text(
    aes(
      x = c(27384341),
      y = df.bbj.IL21R$nlog10p[df.bbj.IL21R$position %in% c(27384341)],
      label = c("rs201121732")
    ),
    hjust = 1.2,
    vjust = 1,
    size = 2.5,
    data = NULL
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.2)))

p.gex.IL21R =
  locusviz::plot_manhattan_panel(
    df.gex.IL21R.NK,
    highlight_pos = highlight_pos.IL21R,
    background.layers = peak_background(peak.IL21R.start, peak.IL21R.end),
    xlim = c(start.IL21R, end.IL21R),
    title = "IL21R eQTL (NK)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none")

peak.ranges.IL21R = legendry::key_range_manual(start = peak.IL21R.start,
                                               end = peak.IL21R.end,
                                               name = c("i", "ii", "iii", "iv", "v", "vi"))

df.IL21R.pip =
  dplyr::bind_rows(
  df.asthma.IL21R,
  df.aiht.IL21R
) %>%
  dplyr::filter(variant %in% c("chr16_27344903_G_A", "chr16_27384341_C_CT", "chr16_27386721_T_C")) %>%
  dplyr::mutate(
    rsid = factor(dplyr::case_when(
      variant == "chr16_27344903_G_A" ~ "rs144651842",
      variant == "chr16_27384341_C_CT" ~ "rs201121732",
      variant == "chr16_27386721_T_C" ~ "rs76715626",
      TRUE ~ NA_character_
    ), levels = c("rs76715626", "rs201121732", "rs144651842"))
  )

p.IL21R.pip.asthma = 
  dplyr::filter(df.IL21R.pip, trait == "J10_ASTHMA_EXMORE") %>%
  ggplot(aes(pip, rsid)) +
  geom_point() +
  locusviz::get_default_theme(
    legend.position = "none",
    hide.xtitle = TRUE,
    hide.ytitle = TRUE
  ) +
  labs(title = "Asthma")

p.IL21R.pip.aiht =
  dplyr::filter(df.IL21R.pip, trait == "E4_HYTHY_AI_STRICT") %>%
  ggplot(aes(pip, rsid)) +
  geom_point() +
  locusviz::get_default_theme(
    legend.position = "none",
    hide.xtitle = FALSE,
    hide.ytitle = TRUE
  ) +
  labs(title = "AIHT", x = "PIP")

df.IL21R.forest =
  dplyr::bind_rows(
    df.asthma.IL21R,
    df.aiht.IL21R
  ) %>%
  dplyr::filter(variant %in% c("chr16_27344903_G_A", "chr16_27384341_C_CT")) %>%
  dplyr::mutate(
    rsid = factor(dplyr::case_when(
      variant == "chr16_27344903_G_A" ~ "rs144651842",
      variant == "chr16_27384341_C_CT" ~ "rs201121732",
      TRUE ~ NA_character_
    ), levels = c("rs201121732", "rs144651842"))
  ) %>%
  dplyr::group_by(trait, rsid, variant) %>%
  dplyr::reframe(
    beta = c(beta, mean),
    sebeta = c(se, sd),
    beta_type = c(rep("Marginal", n()), rep("Posterior", n()))
  )

pd = position_dodge(width = 0.5)
dplyr::filter(df.IL21R.forest, trait == "J10_ASTHMA_EXMORE") %>%
  ggplot(aes(beta, rsid, shape = beta_type)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.3
  ) +
  geom_errorbarh(
    aes(xmin = beta - sebeta, xmax = beta + sebeta),
    height = 0,
    linewidth = 0.3,
    color = "grey50",
    position = pd
  ) +
  geom_point(aes(shape = beta_type), position = pd) +
  locusviz::get_default_theme(
    legend.position = "none",
    hide.xtitle = TRUE,
    hide.ytitle = TRUE
  ) +
  labs(title = "Asthma")

p.link.IL21R =
  dplyr::mutate(df.link.IL21R,
                highlight = (
                  peak_id %in% peaks.IL21R &
                    gene_id == "ENSG00000103522"
                )) %>%
  plot_links(
    start.IL21R,
    end.IL21R,
    highlight_pos.IL21R,
    peak.ranges.IL21R,
    xbreaks = c(27300000, 27350000, 27450000),
    background.layers = peak_background(peak.IL21R.start, peak.IL21R.end, flip_y = TRUE),
  ) +
  labs(x = "Chromosome 16", title = "IL21R Link")

p.gene.IL21R =
  locusviz::plot_gene_panel("chr16",
                            start.IL21R,
                            end.IL21R,
                            genome_build = "hg38",
                            highlight_pos = highlight_pos.IL21R,
                            background.layers = peak_background(peak.IL21R.start, peak.IL21R.end)) +
  theme(axis.title.x = element_blank())



MEF2.motif = readRDS("data/MEF2.rds")
reverseComplementMotif <- function(pwm) {
  rows <- rownames(pwm)
  pwm <- pwm[4:1, ncol(pwm):1]
  rownames(pwm) <- rows
  return(pwm)
}

MEF2D = MEF2.motif[["Hsapiens-HOCOMOCO-MEF2D_f1"]]
MEF2C_rc = reverseComplementMotif(MEF2.motif[["Hsapiens-HOCOMOCO-MEF2C_f1"]])


plot_seqlogo = function(pwm) {
  ggplot() +
    ggseqlogo::geom_logo(pwm) +
    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      data = data.frame(
        xmin = 1.5,
        xmax = 3.5,
        ymin = 0,
        ymax = 2
      ),
      fill = NA,
      linetype = "dashed",
      color = "grey50"
    ) +
    scale_y_continuous(breaks = seq(0, 2),
                       labels = scales::label_number(drop0trailing = TRUE)) +
    ggseqlogo::theme_logo(base_size = 8) +
    theme(
      plot.background = element_blank(),
      plot.margin = margin(0, 0.1, 0, 0.1, unit = "cm"),
      plot.tag = element_text(face = "bold"),
      plot.title = element_text(
        hjust = 4e-3,
        margin = margin(),
        size = 8
      )
    ) +
    labs(y = "Bits")
}

p.MEF2D = plot_seqlogo(MEF2D) + labs(title = "MEF2D")
p.MEF2C_rc = plot_seqlogo(MEF2C_rc) + labs(title = "MEF2C (reverse complement)")

plts.IL21R = list(
  A = p.asthma.IL21R,
  B = p.aiht.IL21R,
  C = p.gex.IL21R,
  D = patchwork::free(p.gene.IL21R, side = "t"),
  E = p.link.IL21R,
  F = p.bbj.IL21R,
  G = p.IL21R.pip.asthma,
  H = p.IL21R.pip.aiht,
  K = p.gex.IL21R.forest,
  L = patchwork::free(p.MEF2D / p.MEF2C_rc, side = "lr")
)

layout.combined =
  c(
    # TNRC18
    # IBD
    patchwork::area(1, 1, 1, 4),
    # TNRC18 eQTL
    patchwork::area(2, 1, 2, 4),
    # Gene
    patchwork::area(3, 1, 3, 4),
    # Link
    patchwork::area(4, 1, 4, 4),
    # Forest plots
    patchwork::area(1, 5, 1, 5),
    patchwork::area(2, 5, 2, 5),
    patchwork::area(3, 5, 4, 5),
    # IL21R
    # Asthma
    patchwork::area(5, 1, 5, 4),
    # AIHT
    patchwork::area(6, 1, 6, 4),
    # IL21R eQTL
    patchwork::area(7, 1, 7, 4),
    # Gene
    patchwork::area(8, 1, 8, 4),
    # Link
    patchwork::area(9, 1, 9, 4),
    # BBJ
    patchwork::area(10, 1, 10, 4),
    # PIP
    patchwork::area(5, 5, 5, 5),
    patchwork::area(6, 5, 6, 5),
    # Forest eQTL
    patchwork::area(7, 5, 7, 5),
    # Seq logos
    patchwork::area(8, 5, 10, 5)
  )
plot(layout.combined)

plt.combined =
  purrr::reduce(c(plts.TNRC18, plts.IL21R), `+`) +
  patchwork::plot_layout(design = layout.combined,
                         heights =
                           c(rep(1, 2), 0.2, 0.7, rep(1, 3), 0.02, 0.7, 1))

cowplot::save_plot(
  "figures/Fig6_locuszoom.pdf",
  plt.combined,
  base_height = 170,
  base_width = 180,
  units = "mm"
)
