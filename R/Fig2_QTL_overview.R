library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.gex.n_samples <- rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/h5ad/integrated_gex_batch1_5.fgid.qc.n_samples.txt"
)
df.atac.n_samples <- rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/h5ad/integrated_atac_batch1_5.fgid.n_samples.txt"
)

#############################
if (!file.exists("data/df.gex.rds")) {
  df.gex = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/saige_qtl/step3/integrated_gex_batch1_5.fgid.qc.*.mean.inv.SAIGE.acat.txt.gz", function(df, remote_path) {
    dplyr::mutate(df, cell_type = parse_cell_type(remote_path))
  })
  saveRDS(df.gex, "data/df.gex.rds")
}

if (!file.exists("data/df.atac.rds")) {
  df.atac = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/atac_results/cis/integrated_atac_batch1_5.fgid.*.sum.inv.cis_qtl_acat.tsv.gz", function(df, remote_path) {
    dplyr::mutate(df, cell_type = parse_cell_type(remote_path))
  })
  saveRDS(df.atac, "data/df.atac.rds")
}

df.gex = readRDS("data/df.gex.rds")
df.atac = readRDS("data/df.atac.rds")

gr.atac =
  dplyr::filter(df.atac, qval < 0.05) %>%
  dplyr::distinct(phenotype_id) %>%
  dplyr::mutate(parse_peak_id(phenotype_id)) %>%
  GenomicRanges::makeGRangesFromDataFrame(
    seqnames.field = "peak_chrom",
    start.field = "peak_start",
    end.field = "peak_end"
  ) %>%
  GenomicRanges::reduce()

sum(GenomicRanges::width(gr.atac)) / 3e9
# 0.0472934

gr.gex =
  dplyr::distinct(df.gex, phenotype_id) %>%
  dplyr::left_join(df.features) %>%
  GenomicRanges::makeGRangesFromDataFrame(
    seqnames.field = "chrom",
    start.field = "start",
    end.field = "end"
  ) %>%
  GenomicRanges::reduce()
sum(GenomicRanges::width(gr.gex)) / 3e9
# 0.583146

dplyr::summarize(
  df.gex,
  n_tested = length(unique(phenotype_id)),
  n_sig = length(unique(phenotype_id[ACAT_q < 0.05])),
  frac_sig = n_sig / n_tested
)
#   n_tested n_sig  frac_sig
# 1    27294 20829 0.7631348

dplyr::summarize(
  df.atac,
  n_tested = length(unique(phenotype_id)),
  n_sig = length(unique(phenotype_id[qval < 0.05])),
  frac_sig = n_sig / n_tested
)
# n_tested  n_sig  frac_sig
# 1   297024 210584 0.7089797

##################
# Fig. 2a-c

count_sig_egenes = function(df, qval_col) {
  dplyr::group_by(df, phenotype_id) %>%
    dplyr::mutate(n_sig_cell_types = sum(!!as.symbol(qval_col) < 0.05)) %>%
    dplyr::group_by(cell_type, n_sig_cell_types) %>%
    dplyr::summarize(n = n(), .groups = "drop") %>%
    dplyr::filter(n_sig_cell_types > 0)
}

df.gex.egenes = dplyr::bind_rows(
  dplyr::filter(df.gex, filter_l1_cell_types(cell_type)) %>%
    count_sig_egenes(qval_col = "ACAT_q"),
  dplyr::filter(df.gex, !filter_l1_cell_types(cell_type)) %>%
    count_sig_egenes(qval_col = "ACAT_q")
)
df.gex.l1.egenes =
  dplyr::filter(df.gex.egenes, filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))
df.gex.l2.egenes =
  dplyr::filter(df.gex.egenes, !filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))

gex.l1.cell_type.order =
  dplyr::group_by(df.gex.l1.egenes, cell_type) %>%
  dplyr::summarize(n = sum(n), .groups = "drop") %>%
  dplyr::arrange(n) %>%
  dplyr::pull(cell_type)
gex.l2.cell_type.order =
  dplyr::group_by(df.gex.l2.egenes, cell_type) %>%
  dplyr::summarize(n = sum(n), .groups = "drop") %>%
  dplyr::arrange(n) %>%
  dplyr::pull(cell_type)

df.atac.egenes = dplyr::bind_rows(
  dplyr::filter(df.atac, filter_l1_cell_types(cell_type)) %>%
    count_sig_egenes(qval_col = "qval"),
  dplyr::filter(df.atac, !filter_l1_cell_types(cell_type)) %>%
    count_sig_egenes(qval_col = "qval")
)
df.atac.l1.egenes =
  dplyr::filter(df.atac.egenes, filter_l1_cell_types(cell_type))  %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))
df.atac.l2.egenes =
  dplyr::filter(df.atac.egenes, !filter_l1_cell_types(cell_type))  %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))

if (!file.exists("tables/ST7_gex_egenes.tsv")) {
  export_table(df.gex.egenes, "tables/ST7_gex_egenes.tsv", "ST7")
}
if (!file.exists("tables/ST11_atac_egenes.tsv")) {
  export_table(df.atac.egenes, "tables/ST11_atac_egenes.tsv", "ST11")
}

plot_n_egenes = function(df,
                         xlab,
                         cell_type.order,
                         labels = l1.labels,
                         fill.breaks = seq(2, 8, by = 2),
                         hide.ylab = FALSE) {
  dplyr::mutate(df, cell_type = factor(cell_type, levels = cell_type.order)) %>%
    ggplot(aes(n, cell_type)) +
    locusviz::or_missing(
      "PBMC" %in% cell_type.order,
      geom_hline(
        yintercept = length(cell_type.order) - 0.5,
        linewidth = 0.25,
        linetype = "dashed",
        color = "grey50"
      )
    ) +
    geom_col(aes(fill = n_sig_cell_types)) +
    scale_x_continuous(labels = scales::label_comma(), expand = expansion()) +
    scale_y_discrete(labels = labels) +
    scale_fill_viridis_c(breaks = fill.breaks, guide = guide_colorbar(title.position = "top")) +
    locusviz::get_default_theme(hide.ytitle = TRUE, hide.ylab = hide.ylab) +
    theme(
      legend.direction = "horizontal",
      legend.key.width = unit(3, "mm"),
      legend.title = element_text(margin = margin(b = 0.5, unit = "mm"))
    ) +
    labs(x = xlab, fill = "# cell types")
}

p.gex.l1.egenes =
  plot_n_egenes(
    df.gex.l1.egenes,
    xlab = "# cis-eGenes",
    cell_type.order = gex.l1.cell_type.order,
    labels = l1.labels
  ) +
  theme(legend.position = "none")
p.gex.l1.egenes

p.atac.l1.egenes =
  plot_n_egenes(
    df.atac.l1.egenes,
    xlab = "# cis-caPeaks",
    labels = l1.labels,
    cell_type.order = gex.l1.cell_type.order,
    hide.ylab = TRUE
  ) +
  theme(legend.position.inside = c(1, 0),
        legend.justification = c(1, 0)) +
  coord_cartesian(xlim = c(0, 200000))
p.atac.l1.egenes


p.gex.l2.egenes =
  plot_n_egenes(
    df.gex.l2.egenes,
    xlab = "# cis-eGenes",
    cell_type.order = gex.l2.cell_type.order,
    labels = l2.labels,
    fill.breaks = seq(5, 25, by = 5)
  ) +
  theme(legend.position = "none")
p.gex.l2.egenes

p.atac.l2.egenes =
  plot_n_egenes(
    df.atac.l2.egenes,
    xlab = "# cis-caPeaks",
    labels = l2.labels,
    fill.breaks = seq(5, 25, by = 5),
    cell_type.order = gex.l2.cell_type.order,
    hide.ylab = TRUE
  ) +
  theme(legend.position.inside = c(1, 0),
        legend.justification = c(1, 0)) +
  coord_cartesian(xlim = c(0, 150000))
p.atac.l2.egenes

df.power.all = dplyr::bind_rows(
  dplyr::group_by(df.gex, cell_type) %>%
    dplyr::summarize(n_tested = n(), n_sig = sum(ACAT_q < 0.05)) %>%
    dplyr::left_join(df.gex.n_samples) %>%
    dplyr::mutate(QTL = "eQTL"),
  dplyr::group_by(df.atac, cell_type) %>%
    dplyr::summarize(n_tested = n(), n_sig = sum(qval < 0.05)) %>%
    dplyr::left_join(df.atac.n_samples) %>%
    dplyr::mutate(QTL = "caQTL")
) %>%
  dplyr::select(QTL, tidyselect::everything())

if (!file.exists("tables/ST12_power_sig_features.tsv")) {
  export_table(df.power.all, "tables/ST12_power_sig_features.tsv", "ST12")
}


df.l1.power =
  dplyr::filter(df.power.all, filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = factor(remove_cell_type_prefix(cell_type), levels = l1.order))

df.l2.power =
  dplyr::filter(df.power.all, !filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = factor(remove_cell_type_prefix(cell_type)))

plot_qtl_power = function(df, legend.position = c(0, 1)) {
  ggplot(df, aes(n_nuclei / n_samples, n_sig)) +
    geom_point(aes(color = cell_type, shape = QTL)) +
    locusviz::get_default_theme(legend.position = legend.position,
                                legend.justification = legend.position) +
    theme(legend.title = element_text(margin = margin(b = 0.5, unit = "mm")),
          legend.spacing = unit(0, "mm")) +
    labs(x = "# nuclei per sample", color = "Cell type", shape = "QTL") +
    scale_color_manual(values = l1.colors, labels = l1.labels) +
    scale_shape_manual(values = qtl.shapes) +
    guides(color = guide_legend(order = 1))
}

p.l1.power =
  plot_qtl_power(df.l1.power) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_log()) +
  coord_cartesian(xlim = c(10, 10000), ylim = c(100, 200000)) +
  labs(y = "# cis-eGenes or caPeaks")

p.l1.power.gex =
  dplyr::filter(df.l1.power, QTL == "eQTL") %>%
  plot_qtl_power(legend.position = c(1, 0.1)) +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous(labels = scales::label_comma()) +
  coord_cartesian(xlim = c(0, 10000)) +
  labs(y = "# cis-eGenes") +
  guides(shape = "none")
p.l1.power.atac =
  dplyr::filter(df.l1.power, QTL == "caQTL") %>%
  plot_qtl_power() +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous(labels = scales::label_comma()) +
  coord_cartesian(xlim = c(0, 10000)) +
  labs(y = "# cis-caPeaks") +
  guides(shape = "none")

p.l1.power.tested =
  dplyr::mutate(df.l1.power, n_sig = n_tested) %>%
  plot_qtl_power() +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_log()) +
  coord_cartesian(xlim = c(10, 10000), ylim = c(100, 300000)) +
  labs(y = "# tested genes/peaks")

################################################################################
# SFig. 7

df.gex.trans = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/gex_results/mean/trans_nominal.sig.phenotypes.tsv.bgz"
)
df.atac.trans = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/atac_results/sum/trans_nominal.sig.phenotypes.tsv.bgz"
)

# Cross-mappability table (Saha & Battle 2018, hg38, GENCODE v26)
if (!file.exists("data/df.crossmap.rds")) {
  df.crossmap <- rgsutil::read_gsfile(
    "gs://expansion_areas/multiome/misc/cross_mappability/hg38_cross_mappability_strength.txt.gz",
    header = FALSE,
    col.names = c("gene_a", "gene_b", "crossmap_count"),
    data.table = TRUE
  ) %>%
    dplyr::mutate(
      gene_a = stringr::str_remove(gene_a, "\\..*$"),
      gene_b = stringr::str_remove(gene_b, "\\..*$")
    )
  df.crossmap <- df.crossmap[, .(crossmap_count = max(crossmap_count)), by = .(gene_a, gene_b)]
  
  saveRDS(df.crossmap, "data/df.crossmap.rds")
}
df.crossmap = readRDS("data/df.crossmap.rds")

flag_cross_map = function(df.gex.trans, df.crossmap) {
  CIS_WINDOW <- 1e6
  
  gr.features <-
    tidyr::drop_na(df.features) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  df.variants <-
    dplyr::distinct(df.gex.trans, top_variant_id) %>%
    tidyr::separate(
      top_variant_id,
      into = c("var_chrom", "var_pos", "var_ref", "var_alt"),
      sep = ":",
      convert = TRUE,
      remove = FALSE
    )
  
  gr.variants <- GenomicRanges::GRanges(
    seqnames = df.variants$var_chrom,
    ranges   = IRanges::IRanges(
      start = pmax(1L, df.variants$var_pos - CIS_WINDOW),
      end   = df.variants$var_pos + CIS_WINDOW
    ),
    top_variant_id = df.variants$top_variant_id
  )
  
  hits <- GenomicRanges::findOverlaps(gr.variants, gr.features)
  dt.var_genes <- data.table::data.table(top_variant_id = gr.variants$top_variant_id[S4Vectors::queryHits(hits)],
                                         cis_gene       = gr.features$phenotype_id[S4Vectors::subjectHits(hits)])
  
  # Symmetrize crossmap so one join covers both homology directions.
  dt.crossmap.sym <- unique(data.table::rbindlist(list(df.crossmap[, .(g1 = gene_a, g2 = gene_b)], df.crossmap[, .(g1 = gene_b, g2 = gene_a)])))
  
  dt.trans_pairs <- unique(data.table::as.data.table(df.gex.trans)[, .(phenotype_id, top_variant_id)])
  
  dt.flagged <- merge(dt.trans_pairs,
                      dt.var_genes,
                      by = "top_variant_id",
                      allow.cartesian = TRUE)[phenotype_id != cis_gene]
  
  dt.flagged <- merge(
    dt.flagged,
    dt.crossmap.sym,
    by.x = c("phenotype_id", "cis_gene"),
    by.y = c("g1", "g2")
  )
  dt.flagged <- unique(dt.flagged[, .(phenotype_id, top_variant_id)])
  
  dplyr::mutate(df.gex.trans,
      crossmap_flag = paste(phenotype_id, top_variant_id) %in%
        paste(dt.flagged$phenotype_id, dt.flagged$top_variant_id)
    )
}

df.gex.trans = flag_cross_map(df.gex.trans, df.crossmap)

cat(
  sprintf(
    "Cross-mappable trans-eQTL associations: %d / %d (%.1f%%)\n",
    sum(df.gex.trans$crossmap_flag),
    nrow(df.gex.trans),
    100 * mean(df.gex.trans$crossmap_flag)
  )
)
# Cross-mappable trans-eQTL associations: 25156 / 69625 (36.1%)

df.gex.trans.thresholds =
  dplyr::group_by(df.gex, cell_type) %>%
  dplyr::summarize(n_genes = n(),
                   threshold = 5e-8 / n_genes,
                   .groups = "drop")

df.atac.trans.thresholds =
  dplyr::group_by(df.atac, cell_type) %>%
  dplyr::summarize(n_genes = n(),
                   threshold = 5e-8 / n_genes,
                   .groups = "drop")

df.gex.trans.bonferroni =
  dplyr::left_join(df.gex.trans, df.gex.trans.thresholds, by = "cell_type") %>%
  dplyr::filter(!crossmap_flag, top_pval < threshold)

df.atac.trans.bonferroni =
  dplyr::left_join(df.atac.trans, df.atac.trans.thresholds, by = "cell_type") %>%
  dplyr::filter(top_pval < threshold)

df.gex.trans.egenes =
  dplyr::bind_rows(
    dplyr::mutate(df.gex.trans, sig = 0) %>%
      dplyr::filter(filter_l1_cell_types(cell_type)) %>%
      count_sig_egenes(qval_col = "sig"),
    dplyr::mutate(df.gex.trans, sig = 0) %>%
      dplyr::filter(!filter_l1_cell_types(cell_type)) %>%
      count_sig_egenes(qval_col = "sig")
  )

df.gex.trans.l1.egenes =
  dplyr::filter(df.gex.trans.egenes, filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))
df.gex.trans.l2.egenes =
  dplyr::filter(df.gex.trans.egenes, !filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))

df.atac.trans.egenes =
  dplyr::bind_rows(
    dplyr::mutate(df.atac.trans, sig = 0) %>%
      dplyr::filter(filter_l1_cell_types(cell_type)) %>%
      count_sig_egenes(qval_col = "sig"),
    plyr::mutate(df.atac.trans, sig = 0) %>%
      dplyr::filter(!filter_l1_cell_types(cell_type)) %>%
      count_sig_egenes(qval_col = "sig"),
  )
df.atac.trans.l1.egenes =
  dplyr::filter(df.atac.trans.egenes, filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))
df.atac.trans.l2.egenes =
  dplyr::filter(df.atac.trans.egenes, !filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))

p.gex.trans.l1.egenes =
  plot_n_egenes(
    df.gex.trans.l1.egenes,
    xlab = "# trans-eGenes",
    cell_type.order = gex.l1.cell_type.order,
    labels = l1.labels
  ) +
  theme(legend.position = "none")
p.gex.trans.l2.egenes =
  plot_n_egenes(
    df.gex.trans.l2.egenes,
    xlab = "# trans-eGenes",
    cell_type.order = gex.l2.cell_type.order,
    labels = l2.labels
  ) +
  theme(legend.position = "none")

p.atac.trans.l1.egenes =
  plot_n_egenes(
    df.atac.trans.l1.egenes,
    xlab = "# trans-caPeaks",
    labels = l1.labels,
    cell_type.order = gex.l1.cell_type.order,
    hide.ylab = TRUE
  ) +
  theme(legend.position.inside = c(1, 0),
        legend.justification = c(1, 0))
p.atac.trans.l2.egenes =
  plot_n_egenes(
    df.atac.trans.l2.egenes,
    xlab = "# trans-caPeaks",
    labels = l2.labels,
    cell_type.order = gex.l2.cell_type.order,
    hide.ylab = TRUE
  ) +
  theme(legend.position.inside = c(1, 0),
        legend.justification = c(1, 0))

p.trans =
  list(
    p.gex.trans.l1.egenes,
    p.atac.trans.l1.egenes,
    p.gex.trans.l2.egenes,
    p.atac.trans.l2.egenes
  ) %>% purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig7_trans_qtl.pdf",
  p.trans,
  base_height = 120,
  base_width = 120,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig7_trans_qtl.png",
  p.trans,
  base_height = 120,
  base_width = 120,
  units = "mm",
  dpi = 300
)

export_table(df.gex.trans.bonferroni,
             "tables/ST10_gex_trans_egenes.tsv",
             "ST10")
export_table(df.atac.trans.bonferroni,
             "tables/ST11_atac_trans_egenes.tsv",
             "ST11")


dplyr::distinct(df.gex.trans, phenotype_id) %>%
  dplyr::left_join(
    dplyr::filter(df.gex, ACAT_q < 0.05) %>%
      dplyr::distinct(phenotype_id) %>%
      dplyr::mutate(mode = "cis")
  ) %>%
  dplyr::summarize(n_trans = n(),
                   n_cis_frac = sum(mode == "cis", na.rm = TRUE) / n())
# n_trans n_cis_frac
# 1   20826  0.8517238

dplyr::distinct(df.atac.trans, phenotype_id) %>%
  dplyr::left_join(
    dplyr::filter(df.atac, qval < 0.05) %>%
      dplyr::distinct(phenotype_id) %>%
      dplyr::mutate(mode = "cis")
  ) %>%
  dplyr::summarize(n_trans = n(),
                   n_cis_frac = sum(mode == "cis", na.rm = TRUE) / n())
# n_trans n_cis_frac
# 1  215949  0.7602026

df.gex.trans.cs = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/saige_qtl/susie/full.with_trans.cs.tsv.bgz"
)
df.atac.trans.cs = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/atac_results/susie/full.with_trans.cs.tsv.bgz"
)

dplyr::left_join(df.gex.cred,
                 df.gex.trans.cs,
                 by = c("phenotype_id" = "phenotype", "cell_type", "cs")) %>%
  dplyr::filter(filter_l1_cell_types(cell_type)) %>%
  dplyr::summarize(
    n_with_trans = sum(!is.na(trans_phenotype_ids)),
    n_cs = n(),
    frac = n_with_trans / n()
  )
# n_with_trans  n_cs       frac
# 1         2704 86536 0.03124711

dplyr::left_join(df.atac.cred,
                 df.atac.trans.cs,
                 by = c("phenotype_id" = "phenotype", "cell_type", "cs")) %>%
  dplyr::filter(filter_l1_cell_types(cell_type)) %>%
  dplyr::summarize(
    n_with_trans = sum(!is.na(trans_phenotype_ids)),
    n_cs = n(),
    frac = n_with_trans / n()
  )
# n_with_trans   n_cs       frac
# 1        37916 636566 0.05956334
################################################################################
# R4C23: same correlation on constrained + highly-expressed subset
dplyr::group_by(df.power.all, QTL) %>%
  dplyr::summarize(cor(n_sig, n_nuclei / n_samples, method = "spearman"))

df.loeuf.v4 <- rgsutil::read_gsfile(
  "gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
) %>%
  dplyr::inner_join(df.features, by = c("gene_id" = "phenotype_id")) %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile)

df.gex.raw.mean <- rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/pseudobulk/integrated_gex_batch1_5.fgid.qc.predicted.celltype.l1.PBMC.mean.raw.bed.gz"
) %>%
  dplyr::transmute(gene_id, mean_expression = rowMeans(dplyr::select(., tidyselect::starts_with("FG")), na.rm = TRUE))

constrained_top_genes <-
  dplyr::inner_join(df.gex.raw.mean, df.loeuf.v4, by = "gene_id") %>%
  dplyr::filter(mean_expression >= quantile(mean_expression, 0.8),
                lof.oe_ci.upper_bin_decile == 0) %>%
  dplyr::pull(gene_id)

df.open4gene.sig <- readRDS("data/open4gene.sig.rds")
constrained_top_peaks <- dplyr::filter(df.open4gene.sig, gene_id %in% constrained_top_genes) %>%
  dplyr::pull(peak_id) %>% unique()

df.power.constrained <- dplyr::bind_rows(
  dplyr::filter(df.gex, phenotype_id %in% constrained_top_genes) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarize(n_tested = n(), n_sig = sum(ACAT_q < 0.05)) %>%
    dplyr::left_join(df.gex.n_samples) %>%
    dplyr::mutate(QTL = "eQTL"),
  dplyr::filter(df.atac, phenotype_id %in% constrained_top_peaks) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarize(n_tested = n(), n_sig = sum(qval < 0.05)) %>%
    dplyr::left_join(df.atac.n_samples) %>%
    dplyr::mutate(QTL = "caQTL")
)

dplyr::filter(df.power.constrained, filter_l1_cell_types(cell_type)) %>%
  dplyr::group_by(QTL) %>%
  dplyr::summarize(rho = cor(n_sig, n_nuclei / n_samples, method = "spearman"))

################################################################################
read_cred = function(df, path) {
  if (nrow(df) == 0) {
    print(path)
    return(NULL)
  }
  dplyr::filter(df, !low_purity) %>%
    dplyr::mutate(
      cell_type = stringr::str_replace(
        trait,
        "^.*\\.(predicted\\.celltype\\.l[12]\\.[^\\.]*)\\..*$",
        "\\1"
      )
    ) %>%
    dplyr::select(-trait) %>%
    dplyr::rename(phenotype_id = region)
}

if (!file.exists("data/df.atac.cred.rds")) {
  df.gex.cred = rgsutil::map_dfr_gsfiles(
    "gs://expansion_areas/multiome/batch1_5/saige_qtl/susie/full/*.SUSIE.cred.bgz",
    read_cred
  )
  saveRDS(df.gex.cred, "data/df.gex.cred.rds")
}

if (!file.exists("data/df.atac.cred.rds")) {
  df.atac.cred = rgsutil::map_dfr_gsfiles(
    "gs://expansion_areas/multiome/batch1_5/atac_results/susie/full/*.SUSIE.cred.bgz",
    read_cred
  )
  saveRDS(df.atac.cred, "data/df.atac.cred.rds")
}

read_high_pip = function(df, path) {
  if (nrow(df) == 0) {
    print(path)
    return(NULL)
  }
  dplyr::filter(df, prob > 0.5) %>%
    dplyr::mutate(
      cell_type = stringr::str_replace(
        trait,
        "^.*\\.(predicted\\.celltype\\.l[12]\\.[^\\.]*)\\..*$",
        "\\1"
      )
    ) %>%
    dplyr::select(
      -trait,-v,-(chromosome:allele2),-cs_specific_prob,-low_purity,-(mean_99:lbf_variable10)
    ) %>%
    dplyr::rename(phenotype_id = region, variant_id = rsid) %>%
    dplyr::select(cell_type, tidyselect::everything())
}

if (!file.exists("data/df.gex.high_pip.rds")) {
  df.gex.high_pip = rgsutil::map_dfr_gsfiles(
    "gs://expansion_areas/multiome/batch1_5/saige_qtl/susie/full/integrated_gex_batch1_5.fgid.qc.predicted.celltype.*.mean.inv.SAIGE.chr*.SUSIE.in_cs.snp.bgz",
    read_high_pip
  )
  saveRDS(df.gex.high_pip, "data/df.gex.high_pip.rds")
}

if (!file.exists("data/df.gex.high_pip.rds")) {
  df.atac.high_pip = rgsutil::map_dfr_gsfiles(
    "gs://expansion_areas/multiome/batch1_5/atac_results/susie/full/integrated_atac_batch1_5.fgid.predicted.celltype.*.sum.inv.chr*.SUSIE.in_cs.snp.bgz",
    read_high_pip
  )
  saveRDS(df.atac.high_pip, "data/df.atac.high_pip.rds")
}

df.gex.cred = readRDS("data/df.gex.cred.rds")
df.atac.cred = readRDS("data/df.atac.cred.rds")
df.gex.high_pip = readRDS("data/df.gex.high_pip.rds")
df.atac.high_pip = readRDS("data/df.atac.high_pip.rds")

dplyr::filter(df.gex.cred, filter_l1_cell_types(cell_type)) %>%
  dplyr::distinct(cell_type, phenotype_id, cs) %>%
  nrow()
# 86536
dplyr::filter(df.gex.cred, !filter_l1_cell_types(cell_type)) %>%
  dplyr::distinct(cell_type, phenotype_id, cs) %>%
  nrow()
# 89350

length(unique(df.gex.high_pip$variant_id))
# 16477

dplyr::filter(df.atac.cred, filter_l1_cell_types(cell_type)) %>%
  dplyr::distinct(cell_type, phenotype_id, cs) %>%
  nrow()
# 636566
dplyr::filter(df.atac.cred, !filter_l1_cell_types(cell_type)) %>%
  dplyr::distinct(cell_type, phenotype_id, cs) %>%
  nrow()
# 506183

length(unique(df.atac.high_pip$variant_id))
# 109163

length(unique(c(
  df.gex.high_pip$variant_id, df.atac.high_pip$variant_id
)))
# 119094

df.gex.n_cs =
  dplyr::group_by(df.gex.cred, cell_type, phenotype_id) %>%
  dplyr::summarize(cs_size = n(), .groups = "drop") %>%
  dplyr::group_by(cell_type, cs_size) %>%
  dplyr::summarize(count = n(), .groups = "drop") %>%
  dplyr::ungroup()

df.atac.n_cs =
  dplyr::group_by(df.atac.cred, cell_type, phenotype_id) %>%
  dplyr::summarize(cs_size = n(), .groups = "drop") %>%
  dplyr::group_by(cell_type, cs_size) %>%
  dplyr::summarize(count = n(), .groups = "drop") %>%
  dplyr::ungroup()

plot_n_cs = function(df,
                     xlab,
                     cell_type.order,
                     labels = l1.labels,
                     hide.ylab = FALSE) {
  dplyr::mutate(df, cell_type = factor(cell_type, levels = cell_type.order)) %>%
    ggplot(aes(count, cell_type, fill = cs_size)) +
    locusviz::or_missing(
      "PBMC" %in% cell_type.order,
      geom_hline(
        yintercept = length(cell_type.order) - 0.5,
        linewidth = 0.25,
        linetype = "dashed",
        color = "grey50"
      )
    ) +
    geom_col() +
    scale_x_continuous(labels = scales::label_comma(), expand = expansion()) +
    scale_y_discrete(labels = labels) +
    scale_fill_viridis_c(
      option = "plasma",
      breaks = seq(5),
      labels = c(seq(4), "5+"),
      limits = c(1, 5),
      oob = scales::squish,
      guide = guide_colorbar(title.position = "top")
    ) +
    locusviz::get_default_theme(hide.ytitle = TRUE, hide.ylab = hide.ylab) +
    theme(
      legend.direction = "horizontal",
      legend.key.width = unit(3, "mm"),
      legend.title = element_text(margin = margin(b = 0.5, unit = "mm"))
    ) +
    labs(x = xlab, fill = "# 95% CS")
}

p.gex.l1.n_cs =
  dplyr::filter(df.gex.n_cs, filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type)) %>%
  plot_n_cs(xlab = "# cis-eGenes", gex.l1.cell_type.order) +
  theme(legend.position = "none")
p.gex.l2.n_cs =
  dplyr::filter(df.gex.n_cs, !filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type)) %>%
  plot_n_cs(xlab = "# cis-eGenes", gex.l2.cell_type.order, labels = l2.labels) +
  theme(legend.position = "none")

p.atac.l1.n_cs =
  dplyr::filter(df.atac.n_cs, filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type)) %>%
  plot_n_cs(xlab = "# cis-caPeaks",
            gex.l1.cell_type.order,
            hide.ylab = TRUE) +
  theme(legend.position.inside = c(1, 0),
        legend.justification = c(1, 0))
p.atac.l2.n_cs =
  dplyr::filter(df.atac.n_cs, !filter_l1_cell_types(cell_type)) %>%
  dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type)) %>%
  plot_n_cs(
    xlab = "# cis-caPeaks",
    gex.l2.cell_type.order,
    labels = l2.labels,
    hide.ylab = TRUE
  ) +
  theme(legend.position.inside = c(1, 0),
        legend.justification = c(1, 0))

p.gex.l1.n_cs + p.atac.l1.n_cs

################################################################################
df.cluster.same = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/coloc/susie/clusters/same_feature_combined_clusters.tsv.gz"
)
df.cluster.all  = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/coloc/susie/clusters/all_feature_combined_clusters.tsv.gz"
)
df.cluster.all.pqtl  = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/coloc/susie/clusters/all_feature_pQTL_combined_clusters.tsv.gz"
)

length(unique(df.cluster.all.pqtl$cluster))
# 193506

df.qtl.upset =
  dplyr::group_by(df.cluster.all.pqtl, cluster) %>%
  dplyr::summarize(
    caQTL = any(QTL == "caQTL"),
    eQTL = any(QTL == "eQTL"),
    pQTL = any(QTL == "pQTL"),
    n_caQTL = length(unique(trait[QTL == "caQTL"])),
    n_eQTL = length(unique(trait[QTL == "eQTL"])),
    n_pQTL = length(unique(trait[QTL == "pQTL"])),
    n_features = length(unique(trait)),
    .groups = "drop"
  ) %>%
  dplyr::ungroup()

caQTL_eQTL_pQTL_clusters =
  dplyr::filter(df.qtl.upset, caQTL & eQTL & pQTL) %>%
  dplyr::pull(cluster)

df.eQTL_pQTL_same =
  dplyr::filter(
    df.cluster.all.pqtl,
    cluster %in% caQTL_eQTL_pQTL_clusters &
      QTL %in% c("eQTL", "pQTL")
  ) %>%
  dplyr::left_join(dplyr::select(df.features, phenotype_id, symbol),
                   by = c("trait" = "phenotype_id")) %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(trait %in% symbol | symbol %in% trait)

length(unique(df.eQTL_pQTL_same$cluster))
# 217

olink.genes = rgsutil::read_gsfile(
  "gs://finngen-production-library-green/omics/proteomics/release_2023_10_11/data/Olink/all_gene_info.txt"
)$geneName

expressed_genes = readRDS("data/expressed_genes.rds")
dplyr::filter(df.features, phenotype_id %in% expressed_genes) %>%
  dplyr::summarize(olink = sum(symbol %in% olink.genes),
                   frac = mean(symbol %in% olink.genes))

df.qtl.upset.rate =
  dplyr::group_by(df.qtl.upset, caQTL, eQTL, pQTL) %>%
  dplyr::summarize(n_cs = n(), frac = n() / nrow(df.qtl.upset), .groups = "drop")
df.qtl.upset.rate
# caQTL eQTL  pQTL       n     frac
# <lgl> <lgl> <lgl>  <int>    <dbl>
#   1 FALSE FALSE TRUE    1480 0.00765
# 2 FALSE TRUE  FALSE  13459 0.0696
# 3 FALSE TRUE  TRUE      78 0.000403
# 4 TRUE  FALSE FALSE 149384 0.772
# 5 TRUE  FALSE TRUE     586 0.00303
# 6 TRUE  TRUE  FALSE  28089 0.145
# 7 TRUE  TRUE  TRUE     430 0.00222

if (!file.exists("tables/ST13_qtl_upset_rate.tsv")) {
  export_table(df.qtl.upset.rate,
               "tables/ST13_qtl_upset_rate.tsv",
               "ST13")
}

dplyr::filter(df.qtl.upset, caQTL & eQTL) %>%
  dplyr::summarize(
    locusviz::median_ci(n_eQTL, colname = "median_eQTL"),
    locusviz::median_ci(n_caQTL, colname = "median_caQTL")
  )

dplyr::summarize(
  df.qtl.upset,
  n = n(),
  n_caQTL_eQTL = sum(caQTL & eQTL),
  frac_caQTL_eQTL = mean(caQTL & eQTL),
  n_caQTL_eQTL_pQTL = sum(caQTL & eQTL & pQTL),
  frac_caQTL_eQTL_pQTL = mean(caQTL & eQTL & pQTL)
)

p.qtl.upset =
  tidyr::pivot_longer(df.qtl.upset,
                      (caQTL:pQTL),
                      names_to = "set",
                      values_to = "value") %>%
  dplyr::filter(value) %>%
  locusviz::plot_upset(
    "cluster",
    "set",
    set_colors = qtl.colors,
    degree_colors = BuenColors::jdb_palette("brewer_blue")[c(3, 5, 7)],
    log10_scale = FALSE,
    return_list = TRUE
  )
p.qtl.upset[[1]] = p.qtl.upset[[1]] + labs(tag = "f", y = "No. molecular features")

dplyr::filter(df.cluster.same, QTL == "eQTL") %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarize(n_l1_cell_types = length(unique(cell_type[filter_l1_cell_types(cell_type)])),
                   n_l2_cell_types = length(unique(cell_type[!filter_l1_cell_types(cell_type)])))

compute_cooccur = function(df, QTL) {
  df.cluster =
    dplyr::filter(df, QTL == .env$QTL) %>%
    dplyr::mutate(cell_type_level = ifelse(filter_l1_cell_types(cell_type), "l1", "l2")) %>%
    dplyr::select(cluster, cell_type, cell_type_level)
  
  dplyr::inner_join(df.cluster,
                    df.cluster,
                    by = "cluster",
                    relationship = "many-to-many") %>%
    dplyr::filter(cell_type_level.x == "l1", cell_type_level.y == "l2") %>%
    dplyr::mutate(cell_type.x = remove_cell_type_prefix(cell_type.x),
                  cell_type.y = remove_cell_type_prefix(cell_type.y)) %>%
    dplyr::count(cell_type.x, cell_type.y) %>%
    dplyr::group_by(cell_type.x) %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::ungroup()
}

get_cooccur_order = function(df, method = "ward.D2") {
  mat <-
    dplyr::select(df, cell_type.x, cell_type.y, frac) %>%
    tidyr::pivot_wider(
      names_from = cell_type.y,
      values_from = frac,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("cell_type.x") %>%
    as.matrix()
  if (method == "cor.dist") {
    cor.dist.l1 <- as.dist(1 - cor(t(mat), use = "pairwise.complete"))
    hc.l1 <- hclust(cor.dist.l1, method = "average")
    cor.dist.l2 <- as.dist(1 - cor(mat, use = "pairwise.complete"))
    hc.l2 <- hclust(cor.dist.l2, method = "average")
  } else {
    hc.l1 <- hclust(dist(mat), method = method)
    hc.l2 <- hclust(dist(t(mat)), method = method)
  }
  
  return(list(
    cell_type.order.x = hc.l1$labels[hc.l1$order],
    cell_type.order.y = hc.l2$labels[hc.l2$order]
  ))
}

plot_cooccur = function(df,
                        cell_type.order.x,
                        cell_type.order.y,
                        hide.ylab = FALSE) {
  dplyr::mutate(
    df,
    cell_type.x = factor(cell_type.x, levels = cell_type.order.x),
    cell_type.y = factor(cell_type.y, levels = cell_type.order.y)
  ) %>%
    ggplot(aes(cell_type.x, cell_type.y, fill = frac)) +
    geom_tile() +
    locusviz::get_default_theme(hide.ylab = hide.ylab) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      legend.justification = "top",
      legend.direction = "horizontal",
      legend.key.width = unit(5, "mm"),
      legend.title = element_text(margin = margin(r = 4), vjust = 1),
      legend.margin = margin(b = -8)
    ) +
    scale_x_discrete(labels = l1.labels, expand = expansion()) +
    scale_y_discrete(labels = l2.labels, expand = expansion()) +
    scale_fill_viridis_c(option = "turbo", limit = c(0, 0.35)) +
    labs(x = "L1 cell types", y = "L2 cell types", fill = "% coloc")
}

df.gex.cooccur <- compute_cooccur(df.cluster.same, "eQTL")
df.atac.cooccur <- compute_cooccur(df.cluster.same, "caQTL")

if (!file.exists("tables/ST7_gex_cooccor.tsv")) {
  dplyr::rename(
    df.gex.cooccur,
    l1_cell_type = cell_type.x,
    l2_cell_type = cell_type.y,
    n_cs = n,
    frac_coloc = frac
  ) %>%
    export_table("tables/ST7_gex_cooccor.tsv", "ST7")
}

if (!file.exists("tables/ST9_atac_cooccor.tsv")) {
  dplyr::rename(
    df.atac.cooccur,
    l1_cell_type = cell_type.x,
    l2_cell_type = cell_type.y,
    n_cs = n,
    frac_coloc = frac
  ) %>%
    export_table("tables/ST9_atac_cooccor.tsv", "ST9")
}

gex.cooccur.order <- get_cooccur_order(df.gex.cooccur, method = "cor.dist")
atac.cooccur.order <- get_cooccur_order(df.atac.cooccur, method = "cor.dist")

p.gex.cooccur = plot_cooccur(
  df.gex.cooccur,
  cell_type.order.x = gex.cooccur.order$cell_type.order.x,
  cell_type.order.y = rev(gex.cooccur.order$cell_type.order.y),
)
p.atac.cooccur = plot_cooccur(
  df.atac.cooccur,
  cell_type.order.x = gex.cooccur.order$cell_type.order.x,
  cell_type.order.y = rev(gex.cooccur.order$cell_type.order.y),
  hide.ylab = TRUE
)


################################################################################
if (!file.exists("data/df.gex.tensorqtl.rds")) {
  df.gex.tensorqtl = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/gex_results/cis/integrated_gex_batch1_5.fgid.qc.*.mean.inv.cis_qtl_acat.tsv.gz", function(df, remote_path) {
    dplyr::mutate(df, cell_type = parse_cell_type(remote_path))
  })
  saveRDS(df.gex.tensorqtl, "data/df.gex.tensorqtl.rds")
}
df.gex.tensorqtl = readRDS("data/df.gex.tensorqtl.rds")

df.gex.cmp = dplyr::left_join(df.gex, df.gex.tensorqtl, by = c("cell_type", "phenotype_id")) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarize(n_saige_qtl = sum(ACAT_q < 0.05),
                   n_tensorqtl = sum(qval < 0.05))

with(df.gex.cmp, mean(n_saige_qtl / n_tensorqtl))
# 1.182499

if (!file.exists("tables/ST8_egenes_tensorqtl.tsv")) {
  export_table(df.gex.cmp, "tables/ST8_egenes_tensorqtl.tsv", "ST8")
}

p.gex.cmp =
  dplyr::mutate(df.gex.cmp,
                cell_type_level = ifelse(filter_l1_cell_types(cell_type), "l1", "l2")) %>%
  ggplot(aes(n_tensorqtl, n_saige_qtl)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = 'grey50'
  ) +
  geom_point(aes(color = cell_type, shape = cell_type_level)) +
  geom_smooth(
    method = "lm",
    formula = y ~ x - 1,
    color = "darkblue",
    show.legend = FALSE
  ) +
  locusviz::get_default_theme(legend.position = c(0, 1),
                              legend.justification = c(0, 1)) +
  scale_x_continuous(expand = expansion(), labels = scales::label_comma()) +
  scale_y_continuous(expand = expansion(), labels = scales::label_comma()) +
  scale_color_manual(values = cell_colors, guide = "none") +
  coord_cartesian(xlim = c(0, 20000), ylim = c(0, 20000)) +
  labs(x = "# eGenes (tensorQTL)", y = "# eGenes (SAIGE-QTL)", shape = "Cell type")
p.gex.cmp

################################################################################

plt.ext.qtl =
  list(
    p.gex.l2.egenes,
    p.gex.l2.n_cs,
    patchwork::free(p.gex.cooccur, side = "t"),
    p.atac.l2.egenes,
    p.atac.l2.n_cs,
    patchwork::free(p.atac.cooccur, side = "t"),
    p.l1.power.gex,
    p.l1.power.atac + theme(legend.position = "none"),
    p.l1.power.tested + theme(legend.position = "none")
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(byrow = FALSE) +
  patchwork::plot_annotation(tag_levels = "a")
plt.ext.qtl

cowplot::save_plot(
  "figures/ExtendedDataFig3_qtl.pdf",
  plt.ext.qtl,
  base_height = 170,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/ExtendedDataFig3_qtl.png",
  plt.ext.qtl,
  base_height = 170,
  base_width = 180,
  units = "mm",
  dpi = 300
)

################################################################################
if (!file.exists("data/df.gex.in_cs.rds")) {
  df.gex.in_cs = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/saige_qtl/susie/full/integrated_gex_batch1_5.fgid.qc.predicted.celltype.l1.PBMC.mean.inv.SAIGE.chr*.SUSIE.in_cs.snp.bgz", function(df, path) {
    dplyr::group_by(df, region) %>%
      dplyr::filter(!low_purity & max(prob) == prob) %>%
      dplyr::ungroup(region)
  })
  saveRDS(df.gex.in_cs, "data/df.gex.in_cs.rds")
}

if (!file.exists("data/df.atac.in_cs.rds")) {
  df.atac.in_cs = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/atac_results/susie/full/integrated_atac_batch1_5.fgid.predicted.celltype.l1.PBMC.sum.inv.chr*.SUSIE.in_cs.snp.bgz", function(df, path) {
    dplyr::group_by(df, region) %>%
      dplyr::filter(!low_purity & max(prob) == prob) %>%
      dplyr::ungroup(region)
  })
  saveRDS(df.atac.in_cs, "data/df.atac.in_cs.rds")
}

df.gex.in_cs = readRDS("data/df.gex.in_cs.rds")
df.atac.in_cs = readRDS("data/df.atac.in_cs.rds")

df.gex.distance =
  dplyr::ungroup(df.gex.in_cs) %>%
  dplyr::left_join(df.features, by = c("region" = "phenotype_id")) %>%
  dplyr::mutate(
    distance = ifelse(strand == "+", position - tss, tss - position),
    distance_within = dplyr::case_when(
      (start <= position) &
        (position <= end) ~ ifelse(strand == "+", position - start, end - position),
      TRUE ~ NA_integer_
    ),
    distance_to_feature = dplyr::case_when(
      position < start ~ ifelse(strand == "+", position - start, start - position),
      end < position ~ ifelse(strand == "+", position - end, end - position),
      TRUE ~ 0
    ),
    width = end - start + 1,
    normalized_distance = distance_within / width
  )

df.atac.distance =
  dplyr::mutate(df.atac.in_cs, parse_peak_id(region)) %>%
  dplyr::mutate(
    peak_pos = (peak_start + peak_start) %/% 2,
    distance = position - peak_pos,
    distance_within = dplyr::case_when(
      (peak_start <= position) &
        (position <= peak_end) ~ position - peak_start,
      TRUE ~ NA_integer_
    ),
    distance_to_feature = dplyr::case_when(
      position < peak_start ~ position - peak_start,
      peak_end < position ~ position - peak_end,
      TRUE ~ 0
    ),
    width = peak_end - peak_start + 1,
    normalized_distance = distance_within / width
  ) %>%
  dplyr::filter(width == 500)

df.distance = dplyr::bind_rows(
  dplyr::transmute(
    df.gex.distance,
    QTL = "eQTL",
    distance = distance,
    distance_to_feature = distance_to_feature,
    width = width,
    normalized_distance = normalized_distance
  ),
  dplyr::transmute(
    df.atac.distance,
    QTL = "caQTL",
    distance = distance,
    distance_to_feature = distance_to_feature,
    width = width,
    normalized_distance = normalized_distance
  )
)

dplyr::group_by(df.distance, QTL) %>%
  dplyr::summarize(median(distance, na.rm = TRUE))
#    QTL   `median(distance, na.rm = TRUE)`
#   <chr>                            <dbl>
# 1 caQTL                              249
# 2 eQTL                              5371

p.distance =
  dplyr::filter(df.distance, abs(distance) <= 1e5) %>%
  ggplot(aes(distance, color = QTL)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_line(stat = "density") +
  geom_text(
    aes(x, y, label = label),
    data = data.frame(
      x = c(-90000, 90000),
      y = 2e-5,
      label = c("5'", "3'")
    ),
    size = 2,
    color = "black"
  ) +
  labs(x = "Distance to TSS or peak summit (Kb)", y = expression("Density (" * "\u00D7" * 10^-5 * ")")) +
  scale_x_continuous(
    labels = scales::number_format(scale = 0.001),
    expand = expansion(),
    breaks = c(-1e5, -5e4, 0, 5e4, 1e5),
    limits = c(-1e5, 1e5)
  ) +
  scale_y_continuous(labels = scales::label_number(scale = 1e5, drop0trailing = TRUE)) +
  scale_color_manual(values = qtl.colors) +
  locusviz::get_default_theme()


distance_breaks = c(-Inf, seq(-1e5, 1e5, by = 10000), Inf)

distance_bin_labels = c("< -100",
                        seq(-100, -20, by = 10),
                        "> -10",
                        "0",
                        "< 10",
                        seq(20, 100, by = 10),
                        "> 100")

df.distance.bin = dplyr::mutate(
  df.distance,
  distance_bin = cut(distance_to_feature, distance_breaks),
  distance_bin_num = as.numeric(distance_bin),
  distance_bin = ifelse(!is.na(normalized_distance), mean(range(distance_bin_num)), distance_bin_num),
  distance_bin = factor(distance_bin, labels = distance_bin_labels)
) %>%
  dplyr::group_by(distance_bin, QTL) %>%
  dplyr::summarize(count = n(), .groups = "drop") %>%
  dplyr::group_by(QTL) %>%
  dplyr::mutate(frac = count / sum(count))


plot_distance_bin = function(df.distance.bin, hide.xtitle = FALSE) {
  ggplot(df.distance.bin, aes(factor(distance_bin), frac, fill = QTL)) +
    geom_col() +
    geom_text(aes(label = scales::percent(frac, accuracy = 1)),
              vjust = -0.5,
              size = 2) +
    scale_fill_manual(values = qtl.colors) +
    locusviz::get_default_theme(hide.xtitle = hide.xtitle) +
    theme(legend.title = element_blank()) +
    scale_x_discrete(limits = distance_bin_labels) +
    scale_y_continuous(expand = expansion(c(0, 0.1)), labels = scales::label_percent()) +
    labs(x = "Binned distance to gene body or peak (Kb)", y = "Fraction")
}

p.distance.bin.caqtl =
  dplyr::filter(df.distance.bin, QTL == "caQTL") %>%
  plot_distance_bin(hide.xtitle = TRUE)
p.distance.bin.eqtl =
  dplyr::filter(df.distance.bin, QTL == "eQTL") %>%
  plot_distance_bin() +
  theme(plot.margin = margin(t = 4))

p.distance.bin =
  p.distance.bin.caqtl + p.distance.bin.eqtl +
  patchwork::plot_layout(ncol = 1) +
  patchwork::plot_annotation(tag_levels = "a")
p.distance.bin

p.norm_distance =
  dplyr::filter(df.distance, !is.na(normalized_distance)) %>%
  ggplot(aes(normalized_distance, color = QTL)) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             color = "grey50") +
  geom_line(stat = "density") +
  geom_text(
    aes(
      x,
      y,
      label = label,
      hjust = hjust,
      color = QTL
    ),
    data = data.frame(
      x = c(rep(c(0.03, 0.97), 3), 0.5),
      y = c(rep(0.15, 2), c(2.2, 1), c(0.7, 0.45, 2)),
      label = c("5'", "3'", "TSS", "TTS", "Peak start", "end", "summit"),
      hjust = c(rep(c(0, 1), 3), 0.5),
      QTL = c(rep("black", 2), rep("eQTL", 2), rep("caQTL", 3))
    ),
    size = 2
  ) +
  labs(x = "Normalized distance within GB or peak", y = "Density") +
  scale_x_continuous(labels = scales::label_number(drop0trailing = TRUE),
                     expand = expansion()) +
  scale_color_manual(values = c(c("black" = "black"), qtl.colors)) +
  locusviz::get_default_theme(legend.position = "none")

p.distance
p.norm_distance

p.distance.inset =
  dplyr::filter(df.distance, abs(distance) <= 5000) %>%
  ggplot(aes(distance, color = QTL)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  stat_density(geom = "line") +
  labs(x = "Distance to TSS/peak summit (Kb)", y = "Density") +
  scale_x_continuous(
    labels = scales::label_number(scale = 0.001, drop0trailing = TRUE),
    expand = expansion(),
    limits = c(-5000, 5000)
  ) +
  scale_color_manual(values = qtl.colors) +
  locusviz::get_default_theme(
    legend.position = "none",
    hide.ylab = TRUE,
    hide.xtitle = TRUE
  ) +
  theme(axis.ticks.y = element_blank())

################################################################################
# enrichment

df.gex.max_pip = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/saige_qtl/susie/full.max_pip.annot.most_severe.txt.bgz"
)
df.atac.max_pip = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/atac_results/susie/full.max_pip.annot.most_severe.txt.bgz"
)

dplyr::filter(df.gex.max_pip,
              max_pip > 0.9 &
                consequence %in% c("pLoF", "Missense", "Synonymous")) %>%
  dplyr::group_by(consequence) %>%
  dplyr::summarize(count = n(), .groups = "drop")

dplyr::filter(df.gex.max_pip, max_pip > 0.9 &
                consequence == "Synonymous") %>%
  dplyr::group_by(most_severe) %>%
  dplyr::summarize(count = n())

dplyr::filter(df.atac.max_pip,
              max_pip > 0.9 &
                consequence %in% c("pLoF", "Missense", "Synonymous")) %>%
  dplyr::group_by(consequence) %>%
  dplyr::summarize(count = n())


with(dplyr::filter(df.gex.max_pip, max_pip > 0.9),
     table(consequence, encode4))

dplyr::filter(df.gex.max_pip, consequence == "Synonymous" &
                max_pip > 0.9)


df.enrichment.consequence =
  dplyr::bind_rows(
    locusviz:::compute_functional_enrichment(df.gex.max_pip, annot_levels = annot_levels) %>%
      dplyr::mutate(QTL = "eQTL"),
    locusviz:::compute_functional_enrichment(df.atac.max_pip, annot_levels = annot_levels) %>%
      dplyr::mutate(QTL = "caQTL")
  )

df.enrichment.encode4 =
  dplyr::bind_rows(
    locusviz:::compute_functional_enrichment(
      df.gex.max_pip,
      consequence_col = "encode4",
      annot_levels = encode4_levels
    )  %>%
      dplyr::mutate(QTL = "eQTL"),
    locusviz:::compute_functional_enrichment(
      df.atac.max_pip,
      consequence_col = "encode4",
      annot_levels = encode4_levels
    ) %>%
      dplyr::mutate(QTL = "caQTL")
  )

ca_link_levels = c("CA-Link+", "CA-Link-")
df.enrichment.ca_link =
  dplyr::bind_rows(
    locusviz:::compute_functional_enrichment(
      df.gex.max_pip,
      consequence_col = "finngen_ca_link",
      annot_levels = ca_link_levels
    ) %>%
      dplyr::mutate(QTL = "eQTL"),
    locusviz:::compute_functional_enrichment(
      df.atac.max_pip,
      consequence_col = "finngen_ca_link",
      annot_levels = ca_link_levels
    ) %>%
      dplyr::mutate(QTL = "caQTL")
  )

ca_link_mechanism_levels = c("Dual", "Switch/Rheostat", "Switch", "Rheostat", "CA-Link-")
df.enrichment.ca_link_mechanism =
  dplyr::bind_rows(
    locusviz:::compute_functional_enrichment(
      dplyr::mutate(
        df.gex.max_pip,
        finngen_ca_link_mechanism = forcats::fct_recode(
          finngen_ca_link_mechanism,
          "Dual" = "Switch & Rheostat",
          "Switch/Rheostat" = "Switch or Rheostat"
        )
      ),
      consequence_col = "finngen_ca_link_mechanism",
      annot_levels = ca_link_mechanism_levels
    ) %>%
      dplyr::mutate(QTL = "eQTL"),
    locusviz:::compute_functional_enrichment(
      dplyr::mutate(
        df.atac.max_pip,
        finngen_ca_link_mechanism = forcats::fct_recode(
          finngen_ca_link_mechanism,
          "Dual" = "Switch & Rheostat",
          "Switch/Rheostat" = "Switch or Rheostat"
        )
      ),
      consequence_col = "finngen_ca_link_mechanism",
      annot_levels = ca_link_mechanism_levels
    ) %>%
      dplyr::mutate(QTL = "caQTL")
  )

coding_cre_levels = stringr::str_c(rep(c("pLoF", "Missense", "Synonymous"), each =
                                         2), rep(c("-cCRE+", "-cCRE-"), 3))
df.enrichment.coding_cre =
  dplyr::bind_rows(
    dplyr::mutate(
      df.gex.max_pip,
      consequence_cre = dplyr::case_when(
        !is.na(encode4) ~ stringr::str_c(consequence, "-cCRE+"),
        TRUE ~ stringr::str_c(consequence, "-cCRE-"),
      )
    ) %>%
      locusviz:::compute_functional_enrichment(consequence_col = "consequence_cre", annot_levels = coding_cre_levels) %>%
      dplyr::mutate(QTL = "eQTL"),
    dplyr::mutate(
      df.atac.max_pip,
      consequence_cre = dplyr::case_when(
        !is.na(encode4) ~ stringr::str_c(consequence, "-cCRE+"),
        TRUE ~ stringr::str_c(consequence, "-cCRE-"),
      )
    ) %>%
      locusviz:::compute_functional_enrichment(consequence_col = "consequence_cre", annot_levels = coding_cre_levels) %>%
      dplyr::mutate(QTL = "caQTL")
    
  )

df.enrichments = rbind(
  df.enrichment.consequence %>% dplyr::mutate(annot = "Consequence"),
  df.enrichment.ca_link %>% dplyr::mutate(annot = "FG"),
  df.enrichment.encode4 %>% dplyr::mutate(annot = "ENCODE4")
) %>%
  dplyr::select(annot, QTL, tidyselect::everything())

if (!file.exists("tables/ST14_functional_enrichments.tsv")) {
  export_table(df.enrichments,
               "tables/ST14_functional_enrichments.tsv",
               "ST14")
}

if (!file.exists("tables/ST15_functional_enrichments_cre.tsv")) {
  dplyr::mutate(
    df.enrichment.coding_cre,
    annot = factor(
      stringr::str_split_fixed(consequence, "-", 2)[, 2],
      levels = c("cCRE+", "cCRE-")
    ),
    consequence = factor(stringr::str_split_fixed(consequence, "-", 2)[, 1], levels = annot_levels)
  ) %>%
    dplyr::select(annot, QTL, consequence, dplyr::everything()) %>%
    export_table("tables/ST15_functional_enrichments_cre.tsv",
                 "ST15",
                 save_googlesheet = FALSE)
}
if (!file.exists("tables/ST16_functional_enrichments_ca_link_mode.tsv")) {
  dplyr::mutate(
    df.enrichment.ca_link_mechanism,
    annot = factor(
      ifelse(consequence == "CA-Link-", "CA-Link-", "CA-Link+"),
      levels = c("CA-Link+", "CA-Link-")
    )
  ) %>%
    dplyr::rename(mode = consequence) %>%
    dplyr::select(annot, QTL, mode, dplyr::everything()) %>%
    export_table(
      "tables/ST16_functional_enrichments_ca_link_mode.tsv",
      "ST16"
    )
}

plot_enrichment = function(df,
                           xlim = NULL,
                           show.hline = TRUE) {
  pd = position_dodge(width = 0.75)
  ggplot(df, aes(enrichment, interaction(consequence, annot), color = QTL)) +
    locusviz::or_missing(show.hline,
                         geom_hline(
                           yintercept = c(length(encode4_levels) - 0.5, length(encode4_levels) + 1.5),
                           linetype = "dashed",
                           color = "grey50"
                         )) +
    geom_vline(xintercept = 1,
               linetype = "dashed",
               color = "grey50") +
    geom_errorbar(
      aes(xmin = lower, xmax = upper),
      width = 0,
      position = pd,
      show.legend = FALSE
    ) +
    geom_point(aes(shape = QTL), position = pd) +
    # scale_x_log10(
    #   expand = expansion(),
    #   breaks = c(0.1, 1, 10, 100),
    #   labels = scales::label_number(drop0trailing = TRUE)
    # ) +
    scale_x_continuous(
      expand = expansion(),
      breaks = c(0.1, 1, seq(5, 25, by = 5)),
      labels = scales::label_number(drop0trailing = TRUE)
    ) +
    scale_y_discrete(
      limits = rev,
      guide = legendry::guide_axis_nested(
        key = legendry::key_range_auto(sep = "\\."),
        levels_text = list(NULL, element_text(angle = 90, hjust = 0.5))
      )
    ) +
    scale_color_manual(values = qtl.colors) +
    scale_shape_manual(values = qtl.shapes) +
    locusviz::get_default_theme(
      hide.ytitle = TRUE,
      legend.position = c(1, 0),
      legend.justification = c(1, 0)
    ) +
    coord_cartesian(xlim = xlim) +
    labs(x = "Enrichment") +
    guides(color = guide_legend(reverse = TRUE),
           shape = guide_legend(reverse = TRUE))
}

p.enrichment =
  dplyr::mutate(
    df.enrichments,
    annot = factor(annot, levels = c("Consequence", "FG", "ENCODE4")),
    consequence = factor(
      consequence,
      levels = c(annot_levels, ca_link_levels, encode4_levels)
    ),
    QTL = factor(QTL, levels = rev(c("caQTL", "eQTL")))
  ) %>%
  plot_enrichment(xlim = c(0.1, 25))

p.enrichment

p.cre.enrichment =
  dplyr::mutate(
    df.enrichment.coding_cre,
    annot = factor(
      stringr::str_split_fixed(consequence, "-", 2)[, 2],
      levels = c("cCRE+", "cCRE-")
    ),
    consequence = factor(stringr::str_split_fixed(consequence, "-", 2)[, 1], levels = annot_levels),
    QTL = factor(QTL, levels = rev(c("caQTL", "eQTL")))
  ) %>%
  plot_enrichment(xlim = c(0.1, 20), show.hline = FALSE)
p.cre.enrichment

p.ca_link_mechanism.enrichment =
  dplyr::mutate(
    df.enrichment.ca_link_mechanism,
    annot = factor(
      ifelse(consequence == "CA-Link-", "", "CA-Link+"),
      levels = c("CA-Link+", "")
    ),
    consequence = consequence,
    QTL = factor(QTL, levels = rev(c("caQTL", "eQTL")))
  ) %>%
  plot_enrichment(xlim = c(0.1, 25), show.hline = FALSE)
p.ca_link_mechanism.enrichment

##################
## heritability
if (!file.exists("data/df.gex.hsq.rds")) {
  df.gex.hsq = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/mesc/saige_qtl/expscore_indiv/integrated_gex_batch1_5.fgid.qc.*.mesc.hsq.gz", function(df, path) {
    dplyr::mutate(df, cell_type = parse_cell_type(path), QTL = "eQTL")
  })
  saveRDS(df.gex.hsq, "data/df.gex.hsq.rds")
}
if (!file.exists("data/df.gex.hsq.rds")) {
  df.atac.hsq = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/mesc/tensorqtl/expscore_indiv/integrated_atac_batch1_5.fgid.*.mesc.hsq.gz", function(df, path) {
    dplyr::mutate(df, cell_type = parse_cell_type(path), QTL = "caQTL")
  })
  saveRDS(df.atac.hsq, "data/df.atac.hsq.rds")
}

df.gex.hsq = readRDS("data/df.gex.hsq.rds")
df.atac.hsq = readRDS("data/df.atac.hsq.rds")

df.hsq = dplyr::bind_rows(
  dplyr::left_join(
    df.gex.hsq,
    df.gex,
    by = c("Gene" = "phenotype_id", "cell_type" = "cell_type")
  ) %>%
    dplyr::filter(ACAT_q < 0.05),
  dplyr::left_join(
    df.atac.hsq,
    df.atac,
    by = c("Gene" = "phenotype_id", "cell_type" = "cell_type")
  ) %>%
    dplyr::filter(qval < 0.05)
) %>%
  tidyr::drop_na(h2cis)

dplyr::bind_rows(
  dplyr::transmute(
    df.gex.in_cs,
    var_explained = 2 * maf * (1 - maf) * beta ** 2,
    QTL = "eQTL"
  ),
  dplyr::transmute(
    df.atac.in_cs,
    var_explained = 2 * maf * (1 - maf) * beta ** 2,
    QTL = "caQTL"
  )
) %>%
  dplyr::group_by(QTL) %>%
  dplyr::summarize(mean(var_explained))
# 1 caQTL               0.0519
# 2 eQTL                0.00997


dplyr::filter(df.hsq, cell_type == "predicted.celltype.l1.PBMC") %>%
  dplyr::group_by(QTL) %>%
  dplyr::summarize(mean(h2cis))
# QTL   `mean(h2cis)`
# <chr>         <dbl>
#   1 caQTL         0.105
# 2 eQTL          0.130

if (!file.exists("tables/ST17_hsq_mean.tsv")) {
  df.hsq.mean.all =
    dplyr::group_by(df.hsq, cell_type, QTL) %>%
    dplyr::summarize(locusviz::mean_ci(h2cis))
  export_table(df.hsq.mean.all, "tables/ST17_hsq_mean.tsv", "ST17")
}

df.hsq.mean = dplyr::filter(df.hsq.mean.all, filter_l1_cell_types(cell_type))

p.hsq.pbmc =
  dplyr::filter(df.hsq, cell_type == "predicted.celltype.l1.PBMC") %>%
  dplyr::mutate(QTL = factor(QTL, levels = c("eQTL", "caQTL"))) %>%
  ggplot(aes(h2cis, QTL, fill = QTL)) +
  ggridges::geom_density_ridges(alpha = 0.8) +
  scale_fill_manual(values = qtl.colors) +
  locusviz::get_default_theme(hide.ytitle = TRUE, legend.position = "none") +
  theme(plot.title = element_text(hjust = 0, margin = margin(t = -8))) +
  labs(title = "PBMC", x = expression(italic(h)[cis]^2)) +
  scale_x_continuous(expand = expansion(),
                     labels = scales::label_number(drop0trailing = TRUE)) +
  scale_y_discrete(expand = expansion())
p.hsq.pbmc

df.line =
  dplyr::mutate(
    df.hsq.mean,
    cell_type = remove_cell_type_prefix(cell_type),
    cell_type = factor(cell_type, levels = gex.l1.cell_type.order),
    QTL = factor(QTL, levels = c("eQTL", "caQTL"))
  ) %>%
  dplyr::select(cell_type, QTL, mean) %>%
  tidyr::pivot_wider(names_from = QTL, values_from = mean)

pd = position_dodge(width = 1)
p.hsq.celltypes =
  dplyr::mutate(
    df.hsq.mean,
    cell_type = remove_cell_type_prefix(cell_type),
    cell_type = factor(cell_type, levels = gex.l1.cell_type.order),
    QTL = factor(QTL, levels = c("eQTL", "caQTL"))
  ) %>%
  ggplot(aes(mean, cell_type, color = QTL)) +
  geom_hline(
    yintercept = 8.5,
    linewidth = 0.25,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_segment(
    data = df.line,
    aes(
      x = eQTL,
      xend = caQTL,
      y = as.numeric(cell_type) - 0.25,
      yend = as.numeric(cell_type) + 0.25
    ),
    color = "grey50",
    linetype = "dotted",
    inherit.aes = FALSE
  ) +
  geom_errorbarh(aes(xmin = mean_lower, xmax = mean_upper),
                 height = 0,
                 position = pd) +
  geom_point(position = pd) +
  scale_color_manual(values = qtl.colors) +
  locusviz::get_default_theme(hide.ytitle = TRUE, legend.position = "none") +
  scale_x_continuous(expand = expansion(),
                     labels = scales::label_number(drop0trailing = TRUE)) +
  scale_y_discrete(labels = l1.labels) +
  labs(x = expression("Mean " * italic(h)[cis]^2)) +
  coord_cartesian(xlim = c(0, max(df.hsq.mean$mean_upper)))
p.hsq.celltypes

####################################

layout = "
ACE
BDF
GIJ
HIK"

plt =
  list(
    A = patchwork::free(p.gex.l1.egenes + labs(tag = "a"), side = "l"),
    B = patchwork::free(p.gex.l1.n_cs + labs(tag = "b"), side = "l"),
    C = patchwork::free(p.atac.l1.egenes + labs(tag = "c"), side = "l"),
    D = patchwork::free(p.atac.l1.n_cs + labs(tag = "d"), side = "l"),
    E = p.l1.power + labs(tag = "e"),
    F = (
      p.qtl.upset %>% purrr::reduce(`+`) + patchwork::plot_layout(ncol = 1, heights = c(1, 0.4))
    ),
    G = p.distance + labs(tag = "g"),
    H = p.norm_distance + labs(tag = "h"),
    I = patchwork::free(p.enrichment + labs(tag = "i")),
    J = p.hsq.pbmc + labs(tag = "j"),
    K = p.hsq.celltypes + labs(tag = "k")
    
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(design = layout, heights = c(1.2, 1.2, 0.6, 0.6))
plt

cowplot::save_plot(
  "figures/Fig2.overview.pdf",
  plt,
  base_height = 170,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/Fig2.overview.png",
  plt,
  base_height = 170,
  base_width = 180,
  units = "mm",
  dpi = 300
)

################################################################################
# rs776736, chr19_54660469_A_G
# tabix -h gs://finngen-production-library-green/omics/proteomics/release_2023_10_11/data/Olink/Finemap/Susie/all/LILRB4.SUSIE.snp.bgz chr19:53660469-55660469 | cat SUSIE.snp.header.txt - | bgzip -c > Olink_2023_10_11.LILRB4.SUSIE.snp.gz

start.LILRB4 <- 54660469 - 1e4
end.LILRB4 <- 54660469 + 1e4
highlight_pos.LILRB4 <- c(54660469)

df.olink.LILRB4 <-
  rgsutil:::fread_wrapper("data/Olink_2023_10_11.LILRB4.SUSIE.snp.gz") %>%
  dplyr::mutate(r2 = 0) %>%
  locusviz::preprocess(variant_col = "rsid",
                       pip_col = "prob",
                       cs_id_col = "cs") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.LILRB4) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)


df.gex.LILRB4 = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  "predicted.celltype.l1.Mono",
  "ENSG00000186818"
)

df.gex.LILRB4.Mono =
  dplyr::filter(df.gex.LILRB4, cell_type == "predicted.celltype.l1.Mono") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = position == highlight_pos.LILRB4) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

p.gex.LILRB4 =
  locusviz::plot_manhattan_panel(
    df.gex.LILRB4.Mono,
    highlight_pos = highlight_pos.LILRB4,
    xlim = c(start.LILRB4, end.LILRB4),
    title = "LILRB4 eQTL (Monocytes)"
  ) + 
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = c(1, 1.2)) +
  geom_text(
    aes(
      x = df.gex.LILRB4.Mono$position[df.gex.LILRB4.Mono$lead_variant],
      y = df.gex.LILRB4.Mono$nlog10p[df.gex.LILRB4.Mono$lead_variant],
      label = "rs776736"
    ),
    vjust = -1,
    size = 2.5,
    data = NULL
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.38), 0)) +
  theme(plot.margin = margin(t = 8))

p.olink.LILRB4 =
  locusviz::plot_manhattan_panel(
    df.olink.LILRB4,
    highlight_pos = highlight_pos.LILRB4,
    xlim = c(start.LILRB4, end.LILRB4),
    title = "LILRB4 pQTL (Olink)"
  ) +
  locusviz::get_default_theme(hide.xtitle = TRUE, legend.position = "none")

p.gene.LILRB4 =
  locusviz::plot_gene_panel(
    "chr19",
    start.LILRB4,
    end.LILRB4,
    genome_build = "hg38",
    highlight_pos = highlight_pos.LILRB4
  )
p.gex.LILRB4 + p.olink.LILRB4 + p.gene.LILRB4 + patchwork::plot_layout(ncol = 1)


################################################################################

layout = "
ABB
CCC
DE#
"

p.ext.qtl2 =
  list(
    p.gex.cmp + labs(tag = "a"),
    ((p.distance.bin.caqtl + labs(tag = "c")) + (p.distance.bin.eqtl + labs(tag = "d")) + patchwork::plot_layout(ncol = 1)
    ),
    (
      p.gex.LILRB4 + labs(tag = "b") + p.olink.LILRB4 + p.gene.LILRB4 + patchwork::plot_layout(ncol = 1, heights = c(1, 1, 0.1))
    ),
    patchwork::free(p.cre.enrichment + labs(tag = "e"), side = "l"),
    patchwork::free(p.ca_link_mechanism.enrichment + labs(tag = "f"), side = "l")
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(design = layout)

cowplot::save_plot(
  "figures/ExtendedDataFig4_qtl2.pdf",
  p.ext.qtl2,
  base_height = 170,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/ExtendedDataFig4_qtl2.png",
  p.ext.qtl2,
  base_height = 170,
  base_width = 180,
  units = "mm",
  dpi = 300
)
