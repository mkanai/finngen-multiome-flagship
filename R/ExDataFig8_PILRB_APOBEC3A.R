library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.open4gene = readRDS("data/open4gene.rds")
df.open4gene.sig = readRDS("data/open4gene.sig.rds")

df.variant = readRDS("data/df.variant.rds")

# rs6979218 - chr7_100295525_C_G
start.PILRB <- 100295525 - 5e4
end.PILRB <- 100295525 + 7e4
highlight_pos.PILRB <- c(100295525)
peak.PILRB.start <- 100295229
peak.PILRB.end <- 100295728

# rs1244787406 - chr19_35901079_T_G
start.TYROBP <- 35901079 - 1e4
end.TYROBP <- 35901079 + 12000
highlight_pos.TYROBP <- c(35901079)

# rs12628403 - chr22_38962032_A_C
start.APOBEC3A <- 38962032 - 1e5
end.APOBEC3A <- 38962032 + 1e5
highlight_pos.APOBEC3A <- c(38962032)

df.alz.PILRB <-
  rgsutil:::fread_wrapper("data/GCST90027158.Alzheimer.PILRB.txt.gz") %>%
  dplyr::mutate(
    chromosome = stringr::str_c("chr", chromosome),
    position = base_pair_location,
    variant = locusviz::variant_str(chromosome, position, other_allele, effect_allele),
    variant2 = locusviz::variant_str(chromosome, position, effect_allele, other_allele),
    r2 = 0,
    pip = NA,
    cs_id = NA
  ) %>%
  locusviz::preprocess(beta_col = "beta",
                       se_col = "standard_error",
                       pvalue_col = "p_value") %>%
  dplyr::mutate(
    lead_variant = position == highlight_pos.PILRB,
    variant = ifelse(lead_variant, "chr7_100295525_C_G", variant)
  ) %>%
  locusviz::annotate_r2(reference_panel = "1000G", window = 2e6)

df.alz.TYROBP <-
  rgsutil:::fread_wrapper("data/G6_ALZHEIMER.TYROBP.SUSIE.snp.gz") %>%
  dplyr::mutate(r2 = 0) %>%
  locusviz::preprocess(variant_col = "rsid",
                       pip_col = "prob",
                       cs_id_col = "cs") %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.gex.PILRB = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  l1.cell_types.primary,
  "ENSG00000121716"
)

df.gex.PILRB.Mono =
  dplyr::filter(df.gex.PILRB, cell_type == "predicted.celltype.l1.Mono") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.PILRB) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.atac.PILRB = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/atac_results/cis_nominal/zfiles/%s/integrated_atac_batch1_5.fgid.%s.sum.inv.%s.%s.z",
  c("predicted.celltype.l1.CD4_T", "predicted.celltype.l1.Mono"),
  "chr7-100295229-100295728"
)

df.link.PILRB =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% l1.cell_types.primary &
      gene_id == "ENSG00000121716" &
      .env$start.PILRB < peak_end & peak_start < .env$end.PILRB
  ) %>%
  dplyr::mutate(
    peak_mid = (peak_start + peak_end) / 2,
    start = ifelse(tss < peak_mid, tss, peak_mid),
    end = ifelse(tss < peak_mid, peak_mid, tss),
    group = peak_id,
    score = hurdle_count_beta,
    beta =  hurdle_count_beta
  )


df.gex.TYROBP.Mono =
  rgsutil:::fread_wrapper("data/rs1244787406.predicted.celltype.l1.Mono.SAIGE.step2.gz") %>%
  dplyr::mutate(
    chromosome = stringr::str_c("chr", CHR),
    r2 = 0,
    pip = NA,
    cs_id = NA
  ) %>%
  locusviz::preprocess(
    position_col = "POS",
    variant_col = "MarkerID",
    beta_col = "BETA",
    se_col = "SE",
    pvalue_col = "p.value"
  ) %>%
  dplyr::mutate(lead_variant = position == highlight_pos.TYROBP) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

# gsutil cat gs://finngen-production-library-green/finngen_R12/finngen_R12_analysis_data/finemap/full/susie/C3_BREAST_EXALLC.SUSIE.snp.bgz | zcat | head -n1 > SUSIE.snp.header.txt
# tabix -h gs://finngen-production-library-green/finngen_R12/finngen_R12_analysis_data/finemap/full/susie/C3_BREAST_EXALLC.SUSIE.snp.bgz chr22:38862032-39062032 | cat SUSIE.snp.header.txt - | bgzip -c > C3_BREAST_EXALLC.APOBEC3A.SUSIE.snp.gz

# curl http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000359/QTD000359.cc.tsv.gz | zcat | head -n1 > APOBEC3A.header.txt
# tabix -h http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000359/QTD000359.cc.tsv.gz 22:38862032-39062032 | cat APOBEC3A.header.txt - | bgzip -c > eQTL_Catalogue_GTEx_Whole_Blood.APOBEC3A.txrev.gz

APOBEC3A_CBX6.cascade.peaks =
  dplyr::filter(df.variant, variant_id == "chr22_38962032_A_C") %>%
  dplyr::pull(cascade_peak_genes) %>%
  stringr::str_split(",") %>%
  unlist() %>%
  stringr::str_split_fixed("->", 2) %>%
  .[, 1] %>%
  unique()

df.brc.APOBEC3A <-
  rgsutil:::fread_wrapper("data/C3_BREAST_EXALLC.APOBEC3A.SUSIE.snp.gz") %>%
  dplyr::mutate(r2 = 0) %>%
  locusviz::preprocess(variant_col = "rsid",
                       pip_col = "prob",
                       cs_id_col = "cs") %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.gex.APOBEC3A = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  setdiff(l1.cell_types.primary, "predicted.celltype.l1.B"),
  "ENSG00000128383"
)

df.gex.APOBEC3A.Mono =
  dplyr::filter(df.gex.APOBEC3A, cell_type == "predicted.celltype.l1.Mono") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.APOBEC3A) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.link.APOBEC3A =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% l1.cell_types.primary &
      gene_id == "ENSG00000128383" &
      .env$start.APOBEC3A < peak_end &
      peak_start < .env$end.APOBEC3A
  ) %>%
  dplyr::mutate(
    peak_mid = (peak_start + peak_end) / 2,
    start = ifelse(tss < peak_mid, tss, peak_mid),
    end = ifelse(tss < peak_mid, peak_mid, tss),
    group = peak_id,
    score = hurdle_count_beta,
    beta =  hurdle_count_beta
  )

df.gex.CBX6 = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  l1.cell_types.primary,
  "ENSG00000183741"
)

df.gex.CBX6.Mono =
  dplyr::filter(df.gex.CBX6, cell_type == "predicted.celltype.l1.Mono") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.APOBEC3A) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.link.CBX6 =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% l1.cell_types.primary &
      gene_id == "ENSG00000183741" &
      .env$start.APOBEC3A < peak_end &
      peak_start < .env$end.APOBEC3A
  ) %>%
  dplyr::mutate(
    peak_mid = (peak_start + peak_end) / 2,
    start = ifelse(tss < peak_mid, tss, peak_mid),
    end = ifelse(tss < peak_mid, peak_mid, tss),
    group = peak_id,
    score = hurdle_count_beta,
    beta =  hurdle_count_beta
  )

df.txrev.APOBEC3A = rgsutil:::fread_wrapper("data/eQTL_Catalogue_GTEx_Whole_Blood.APOBEC3A.txrev.gz")

df.txrev.APOBEC3A.downstream = 
  dplyr::filter(df.txrev.APOBEC3A, molecular_trait_id == "ENSG00000128383.grp_1.downstream.ENST00000495988") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "variant") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.APOBEC3A) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.txrev.APOBEC3B.downstream = 
  dplyr::filter(df.txrev.APOBEC3A, molecular_trait_id == "ENSG00000179750.grp_1.downstream.ENST00000333467") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "variant") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.APOBEC3A) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

################################################################################
p.alz.PILRB =
  locusviz::plot_manhattan_panel(
    df.alz.PILRB,
    highlight_pos = highlight_pos.PILRB,
    background.layers = peak_background(peak.PILRB.start, peak.PILRB.end),
    xlim = c(start.PILRB, end.PILRB),
    title = "Alzheimer's disease (Bellenguez, C et al, 2022)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = c(1, 1.2)) +
  geom_text(
    aes(
      x = df.alz.PILRB$position[df.alz.PILRB$lead_variant],
      y = df.alz.PILRB$nlog10p[df.alz.PILRB$lead_variant],
      label = "rs6979218"
    ),
    vjust = -1,
    size = 2.5,
    data = NULL
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.38), 0))

p.gex.PILRB =
  locusviz::plot_manhattan_panel(
    df.gex.PILRB.Mono,
    highlight_pos = highlight_pos.PILRB,
    background.layers = peak_background(peak.PILRB.start, peak.PILRB.end),
    xlim = c(start.PILRB, end.PILRB),
    title = "PILRB eQTL (Monocytes)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none")

p.link.PILRB =
  dplyr::mutate(
    df.link.PILRB,
    highlight = (
      peak_id == "chr7-100295229-100295728" &
        gene_id == "ENSG00000121716"
    )
  ) %>%
  plot_links(
    start.PILRB,
    end.PILRB,
    highlight_pos.PILRB,
    # xbreaks = setdiff(seq(52e5, 56e5, by = 1e5), 54e5),
    hide.xtitle = FALSE,
    background.layers = peak_background(peak.PILRB.start, peak.PILRB.end, flip_y = TRUE),
  ) +
  labs(x = "Chromosome 7", title = "PILRB Link")

p.gene.PILRB =
  locusviz::plot_gene_panel(
    "chr7",
    start.PILRB,
    end.PILRB,
    genome_build = "hg38",
    highlight_pos = highlight_pos.PILRB,
    background.layers = peak_background(peak.PILRB.start, peak.PILRB.end)
  ) +
  theme(axis.title.x = element_blank())

p.alz.TYROBP =
  locusviz::plot_manhattan_panel(
    df.alz.TYROBP,
    highlight_pos = highlight_pos.TYROBP,
    xlim = c(start.TYROBP, end.TYROBP),
    title = "Alzheimer's disease (FinnGen R12)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = c(1, 1.2)) +
  geom_text(
    aes(
      x = df.alz.TYROBP$position[df.alz.TYROBP$lead_variant],
      y = df.alz.TYROBP$nlog10p[df.alz.TYROBP$lead_variant],
      label = "rs1244787406"
    ),
    vjust = -1,
    size = 2.5,
    data = NULL
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.38), 0)) +
  guides(x.sec = legendry::primitive_bracket(
    legendry::key_range_manual(
      start = 35905673,
      end = 35905673 + 5200,
      name = c("5.2 kb deletion")
    )
  )) + theme(plot.title = element_text(margin = margin(t = 8, b = -12)))
p.alz.TYROBP

p.gex.TYROBP =
  locusviz::plot_manhattan_panel(
    df.gex.TYROBP.Mono,
    highlight_pos = highlight_pos.TYROBP,
    xlim = c(start.TYROBP, end.TYROBP),
    title = "TYROBP eQTL (Monocytes)"
  ) +
  locusviz::get_default_theme(hide.xtitle = TRUE, legend.position = "none")

p.gene.TYROBP =
  locusviz::plot_gene_panel(
    "chr19",
    start.TYROBP,
    end.TYROBP,
    genome_build = "hg38",
    highlight_pos = highlight_pos.TYROBP,
    highlight_pos_y = 2
  )

p.brc.APOBEC3A =
  locusviz::plot_manhattan_panel(
    df.brc.APOBEC3A,
    highlight_pos = highlight_pos.APOBEC3A,
    # background.layers = peak_background(peak.APOBEC3A.start, peak.APOBEC3A.end),
    xlim = c(start.APOBEC3A, end.APOBEC3A),
    title = "Breast cancer (FinnGen R12)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = c(1, 1.2)) +
  geom_text(
    aes(
      x = df.brc.APOBEC3A$position[df.brc.APOBEC3A$lead_variant],
      y = df.brc.APOBEC3A$nlog10p[df.brc.APOBEC3A$lead_variant],
      label = "rs12628403"
    ),
    vjust = -1,
    size = 2.5,
    data = NULL
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.38), 0)) + theme(plot.margin = margin(t = 4))

p.gex.APOBEC3A =
  locusviz::plot_manhattan_panel(
    df.gex.APOBEC3A.Mono,
    highlight_pos = highlight_pos.APOBEC3A,
    xlim = c(start.APOBEC3A, end.APOBEC3A),
    title = "APOBEC3A eQTL (Monocytes)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none")

p.gex.CBX6 =
  locusviz::plot_manhattan_panel(
    df.gex.CBX6.Mono,
    highlight_pos = highlight_pos.APOBEC3A,
    xlim = c(start.APOBEC3A, end.APOBEC3A),
    title = "CBX6 eQTL (Monocytes)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none")

p.link.APOBEC3A =
  dplyr::mutate(
    df.link.APOBEC3A,
    highlight = (
      peak_id %in% APOBEC3A_CBX6.cascade.peaks &
        gene_id == "ENSG00000128383"
    )
  ) %>%
  plot_links(start.APOBEC3A,
             end.APOBEC3A,
             highlight_pos.APOBEC3A,
             # xbreaks = setdiff(seq(52e5, 56e5, by = 1e5), 54e5),
             hide.xtitle = FALSE) +
  labs(x = "Chromosome 22", title = "APOBEC3A Link")

p.link.CBX6 =
  dplyr::mutate(
    df.link.CBX6,
    highlight = (
      peak_id %in% APOBEC3A_CBX6.cascade.peaks &
        gene_id == "ENSG00000183741"
    )
  ) %>%
  plot_links(start.APOBEC3A,
             end.APOBEC3A,
             highlight_pos.APOBEC3A,
             # xbreaks = setdiff(seq(52e5, 56e5, by = 1e5), 54e5),
             hide.xtitle = FALSE) +
  labs(x = "Chromosome 22", title = "CBX6 Link")

p.gene.APOBEC3A =
  locusviz::plot_gene_panel(
    "chr22",
    start.APOBEC3A,
    end.APOBEC3A,
    genome_build = "hg38",
    highlight_pos = highlight_pos.APOBEC3A
  )


p.PILRB_APOBEC3A =
  list(
    p.alz.PILRB,
    p.gex.PILRB,
    p.gene.PILRB,
    p.link.PILRB,
    p.alz.TYROBP,
    p.gex.TYROBP,
    p.gene.TYROBP,
    p.brc.APOBEC3A,
    p.gex.APOBEC3A,
    p.gex.CBX6,
    p.gene.APOBEC3A + theme(axis.title.x = element_blank()),
    p.link.APOBEC3A,
    p.link.CBX6
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1,
                         heights = c(rep(1, 2), 0.2, 0.7, rep(1, 2), 0.2, rep(1, 3), 0.2, 0.7, 0.7))

cowplot::save_plot(
  "figures/ExtendedDataFig8_PILRB_locuszoom.pdf",
  p.PILRB_APOBEC3A,
  base_height = 170,
  base_width = 180,
  units = "mm"
)

################################################################################
p.txrev.APOBEC3A.downstream =
  locusviz::plot_manhattan_panel(
    df.txrev.APOBEC3A.downstream,
    highlight_pos = highlight_pos.APOBEC3A,
    xlim = c(start.APOBEC3A, end.APOBEC3A),
    title = "APOBEC3A txrev downstream eQTL (GTEx Whole Blood)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none")
p.txrev.APOBEC3B.downstream =
  locusviz::plot_manhattan_panel(
    df.txrev.APOBEC3B.downstream,
    highlight_pos = highlight_pos.APOBEC3A,
    xlim = c(start.APOBEC3A, end.APOBEC3A),
    title = "APOBEC3B txrev downstream eQTL (GTEx Whole Blood)"
  ) +
  locusviz::get_default_theme(hide.xtitle = TRUE, legend.position = "none")

p.txrev.APOBEC3A_3B =
  list(
    p.txrev.APOBEC3A.downstream + labs(tag = "a"),
    p.txrev.APOBEC3B.downstream  + labs(tag = "b"),
    p.gene.APOBEC3A
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1,
                         heights = c(rep(1, 2), 0.2))

cowplot::save_plot(
  "figures/SFig15_APOBEC3A_3B.pdf",
  p.txrev.APOBEC3A_3B,
  base_height = 100,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig15_APOBEC3A_3B.png",
  p.txrev.APOBEC3A_3B,
  base_height = 100,
  base_width = 180,
  units = "mm",
  dpi = 300
)
################################################################################

# rs6979218 - chr7_100295525_C_G
p.gex.PILRB.forest =
  dplyr::filter(df.gex.PILRB, rsid == "chr7_100295525_C_G") %>%
  plot_forest(legend.position = c(1, 0),
              legend.justification = c(1, 0)) +
  labs(title = "PILRB eQTL (rs6979218)")

p.atac.PILRB.forest =
  dplyr::filter(df.atac.PILRB, rsid == "chr7_100295525_C_G") %>%
  plot_forest(legend.position = c(1, 0),
              legend.justification = c(1, 0)) +
  labs(title = "caQTL (rs6979218)")

p.link.PILRB.forest =
  dplyr::filter(df.open4gene,
                peak_id == "chr7-100295229-100295728" &
                  gene_id == "ENSG00000121716") %>%
  dplyr::rename(beta = hurdle_count_beta, se = hurdle_count_se) %>%
  dplyr::filter(cell_type %in% l1.cell_types.primary) %>%
  dplyr::mutate(sig = ifelse(sig, "FDR < 0.05", "FDR >= 0.05")) %>%
  plot_forest(legend.position = c(1, 0),
              legend.justification = c(1, 0)) +
  labs(title = "Link")

# rs12628403 - chr22_38962032_A_C
p.gex.APOBEC3A.forest =
  dplyr::filter(df.gex.APOBEC3A, rsid == "chr22_38962032_A_C") %>%
  plot_forest() +
  labs(title = "APOBEC3A eQTL (rs12628403)")

p.gex.CBX6.forest =
  dplyr::filter(df.gex.CBX6, rsid == "chr22_38962032_A_C") %>%
  plot_forest() +
  labs(title = "CBX6 eQTL (rs12628403)")

p.forest = list(
  p.gex.PILRB.forest,
  p.atac.PILRB.forest,
  p.link.PILRB.forest,
  p.gex.APOBEC3A.forest,
  p.gex.CBX6.forest
) %>%
  purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a")
p.forest

cowplot::save_plot(
  "figures/SFig14_PILRB_forest.pdf",
  p.forest,
  base_height = 120,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig14_PILRB_forest.png",
  p.forest,
  base_height = 120,
  base_width = 180,
  units = "mm",
  dpi = 300
)
