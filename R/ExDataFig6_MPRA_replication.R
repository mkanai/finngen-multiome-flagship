library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.variant = readRDS("data/df.variant.rds") %>%
  dplyr::mutate(variant_id = stringr::str_replace_all(variant_id, "_", ":")) %>%
  dplyr::mutate(qtl_mechanism = qtl_pattern)

df.mpra = rgsutil:::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/MPRA_Siraj_etal/core_mpra.hg38.tsv.bgz"
)

if (!file.exists("data/df.mpra_per_ct.rds")) {
  df.mpra_per_ct <-
    dplyr::group_by(df.mpra, variant_hg38, cell_type) %>%
    dplyr::arrange(-emVar, -pmin(mean_RNA_alt, mean_RNA_ref)) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::ungroup()
  saveRDS(df.mpra_per_ct, "data/df.mpra_per_ct.rds")
}
df.mpra_per_ct = readRDS("data/df.mpra_per_ct.rds")

df.mpra_meta <- rgsutil:::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/MPRA_Siraj_etal/data/preprocess/mpra_meta.txt.gz"
)
df.mpra_k562_meta <-
  dplyr::filter(df.mpra_per_ct, cell_type == "K562") %>%
  dplyr::left_join(df.mpra_meta, by = c("variant", "cohort")) %>%
  dplyr::mutate(
    log2Skew_meta = ifelse(ref_alt_flip, -log2Skew_meta, log2Skew_meta),
    log2FC_meta = ifelse(ref_alt_flip, -log2FC_meta, log2FC_meta)
  )

df.mpra.eqtl.max_pip = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/MPRA_Siraj_etal/fg_eqtl.max_pip.tsv.bgz"
)
df.mpra.caqtl.max_pip = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/MPRA_Siraj_etal/fg_caqtl.max_pip.tsv.bgz"
)

pip_bin_breaks <- c(-Inf, 0.01, 0.1, 0.5, 0.9, 1.0)
df.mpra_k562_meta.eqtl = dplyr::left_join(df.mpra_k562_meta,
                                          df.mpra.eqtl.max_pip,
                                          by = c("variant_hg38" = "variant")) %>%
  dplyr::filter(type == "null_PIP10" | !is.na(eqtl.max_pip)) %>%
  dplyr::mutate(
    max_pip_bin = cut(eqtl.max_pip, pip_bin_breaks),
    max_pip_bin = factor(
      ifelse(is.na(max_pip_bin), "Null", as.character(max_pip_bin)),
      levels = c("Null", levels(max_pip_bin))
    ),
    max_pip_bin = forcats::fct_recode(max_pip_bin, "[0,0.01]" = "(-Inf,0.01]"),
    QTL = "eQTL"
  )

df.mpra_k562_meta.caqtl = dplyr::left_join(df.mpra_k562_meta,
                                           df.mpra.caqtl.max_pip,
                                           by = c("variant_hg38" = "variant")) %>%
  dplyr::filter(type == "null_PIP10" | !is.na(caqtl.max_pip)) %>%
  dplyr::mutate(
    max_pip_bin = cut(caqtl.max_pip, pip_bin_breaks),
    max_pip_bin = factor(
      ifelse(is.na(max_pip_bin), "Null", as.character(max_pip_bin)),
      levels = c("Null", levels(max_pip_bin))
    ),
    max_pip_bin = forcats::fct_recode(max_pip_bin, "[0,0.01]" = "(-Inf,0.01]"),
    QTL = "caQTL"
  )

dplyr::bind_rows(df.mpra_k562_meta.eqtl, df.mpra_k562_meta.caqtl) %>%
  dplyr::filter(max_pip_bin != "NULL") %>%
  dplyr::distinct(variant_hg38) %>%
  nrow()
# 172042
dplyr::bind_rows(df.mpra_k562_meta.eqtl, df.mpra_k562_meta.caqtl) %>%
  dplyr::filter(eqtl.max_pip > 0.5 | caqtl.max_pip > 0.5) %>%
  dplyr::distinct(variant_hg38) %>%
  nrow()
# 10428

df.max_pip_bin.emvar =
  dplyr::bind_rows(df.mpra_k562_meta.eqtl, df.mpra_k562_meta.caqtl) %>%
  dplyr::select(QTL, variant_hg38, max_pip_bin, emVar, emVar_meta) %>%
  tidyr::pivot_longer(-c(QTL, variant_hg38, max_pip_bin)) %>%
  dplyr::mutate(cell_type = ifelse(name == "emVar", "K562", "Meta")) %>%
  dplyr::group_split(QTL, max_pip_bin, cell_type) %>%
  purrr::map_dfr(function(data) {
    n_total <- nrow(data)
    n_emvar <- sum(data$value)
    
    # Return tibble with all results
    tibble::tibble(
      QTL = data$QTL[1],
      cell_type = data$cell_type[1],
      max_pip_bin = data$max_pip_bin[1],
      n_emvar = n_emvar,
      n_total = n_total,
      locusviz::binom_ci(n_emvar, n_total)
    ) %>%
      dplyr::rename(frac_emvar = frac,
                    lower = frac_lower,
                    upper = frac_upper)
  })


null_variants = dplyr::filter(df.mpra_k562_meta.eqtl, type == "null_PIP10") %>%
  dplyr::pull(variant_hg38) %>%
  unique()

df.atac.mpra.beta = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/MPRA_Siraj_etal/fg_caqtl.mpra_beta.tsv.bgz"
) %>%
  dplyr::mutate(pip = ifelse(variant %in% null_variants, -log10(p), pip)) %>%
  dplyr::group_by(variant) %>%
  dplyr::filter(pip == max(pip)) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(df.mpra_k562_meta.caqtl, by = c("variant" = "variant_hg38")) %>%
  dplyr::rename(variant_hg38 = variant)

df.gex.mpra.beta = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/MPRA_Siraj_etal/fg_eqtl.mpra_beta.tsv.bgz"
) %>%
  dplyr::mutate(pip = ifelse(variant %in% null_variants, -log10(p), pip)) %>%
  dplyr::group_by(variant) %>%
  dplyr::filter(pip == max(pip)) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(df.mpra_k562_meta.eqtl, by = c("variant" = "variant_hg38")) %>%
  dplyr::rename(variant_hg38 = variant)

df.max_pip_bin.log2skew_rho =
  dplyr::bind_rows(df.atac.mpra.beta, df.gex.mpra.beta) %>%
  dplyr::select(QTL,
                variant_hg38,
                max_pip_bin,
                beta,
                emVar,
                emVar_meta,
                log2Skew,
                log2Skew_meta) %>%
  tidyr::pivot_longer(-c(QTL, variant_hg38, max_pip_bin, beta, log2Skew, log2Skew_meta)) %>%
  dplyr::mutate(
    cell_type = ifelse(name == "emVar", "K562", "Meta"),
    log2Skew = ifelse(cell_type == "K562", log2Skew, log2Skew_meta)
  ) %>%
  dplyr::group_split(QTL, max_pip_bin, cell_type) %>%
  purrr::map_dfr(function(data) {
    ct = with(data, cor.test(rank(beta), rank(log2Skew)))
    ct_emvar = with(subset(data, value),
                    cor.test(rank(beta), rank(log2Skew)))
    
    # Return tibble with all results
    tibble::tibble(
      QTL = data$QTL[1],
      cell_type = data$cell_type[1],
      max_pip_bin = data$max_pip_bin[1],
      emVar = c(FALSE, TRUE),
      rho = c(ct$estimate, ct_emvar$estimate),
      lower = c(ct$conf.int[1], ct_emvar$conf.int[1]),
      upper = c(ct$conf.int[2], ct_emvar$conf.int[2])
    )
  })

df.cascade.emvar =
  dplyr::select(df.mpra_k562_meta, variant_hg38, emVar, emVar_meta) %>%
  tidyr::pivot_longer(-variant_hg38) %>%
  dplyr::mutate(cell_type = ifelse(name == "emVar", "K562", "Meta")) %>%
  dplyr::inner_join(df.variant, by = c("variant_hg38" = "variant_id")) %>%
  dplyr::group_split(qtl_mechanism_category, cell_type) %>%
  purrr::map_dfr(function(data) {
    n_total <- nrow(data)
    n_emvar <- sum(data$value)
    
    # Return tibble with all results
    tibble::tibble(
      cell_type = data$cell_type[1],
      qtl_mechanism_category = data$qtl_mechanism_category[1],
      n_emvar = n_emvar,
      n_total = n_total,
      locusviz::binom_ci(n_emvar, n_total)
    ) %>%
      dplyr::rename(frac_emvar = frac,
                    lower = frac_lower,
                    upper = frac_upper)
  })


if (!file.exists("data/df.gex.in_cs.mpra.tsv.gz")) {
  mpra_variants = stringr::str_replace_all(unique(df.mpra$variant_hg38), ":", "_")
  df.gex.in_cs.mpra = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/saige_qtl/susie/full/integrated_gex_batch1_5.fgid.qc.*.mean.inv.SAIGE.chr*.SUSIE.in_cs.snp.bgz", function(df, path) {
    if (nrow(df) == 0) {
      return(NULL)
    }
    dplyr::filter(df, rsid %in% mpra_variants) %>%
      dplyr::mutate(cell_type = parse_cell_type(trait)) %>%
      dplyr::select(cell_type, region, rsid, beta, se, prob, cs)
  })
  data.table::fwrite(
    df.gex.in_cs.mpra,
    "data/df.gex.in_cs.mpra.tsv.gz",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    na = "NA"
  )
}
df.gex.in_cs.mpra = data.table::fread("data/df.gex.in_cs.mpra.tsv.gz", data.table = FALSE)

if (!file.exists("data/df.atac.in_cs.mpra.tsv.gz")) {
  mpra_variants = stringr::str_replace_all(unique(df.mpra$variant_hg38), ":", "_")
  df.atac.in_cs.mpra = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/atac_results/susie/full/integrated_atac_batch1_5.fgid.*.sum.inv.chr*.SUSIE.in_cs.snp.bgz", function(df, path) {
    if (nrow(df) == 0) {
      return(NULL)
    }
    dplyr::filter(df, rsid %in% mpra_variants) %>%
      dplyr::mutate(cell_type = parse_cell_type(trait)) %>%
      dplyr::select(cell_type, region, rsid, beta, se, prob, cs)
  })
  data.table::fwrite(
    df.atac.in_cs.mpra,
    "data/df.atac.in_cs.mpra.tsv.gz",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    na = "NA"
  )
}
df.atac.in_cs.mpra = data.table::fread("data/df.atac.in_cs.mpra.tsv.gz", data.table = FALSE)

df.open4gene = readRDS("data/open4gene.rds")
df.open4gene.sig = readRDS("data/open4gene.sig.rds")

expressed_genes = readRDS("data/expressed_genes.rds")
df.loeuf.v4 = rgsutil::read_gsfile(
  "gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
) %>%
  dplyr::inner_join(df.features, by = c("gene_id" = "phenotype_id")) %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile) %>%
  dplyr::filter(gene_id %in% expressed_genes)

df.atac.in_cs.mpra.max_prob =
  dplyr::filter(df.atac.in_cs.mpra, prob > 0.5) %>%
  dplyr::inner_join(
    df.open4gene.sig,
    by = c("region" = "peak_id", "cell_type"),
    relationship = "many-to-many"
  ) %>%
  dplyr::left_join(df.loeuf.v4, by = "gene_id") %>%
  dplyr::group_by(rsid) %>%
  dplyr::filter(prob == max(prob)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(variant = stringr::str_replace_all(rsid, "_", ":"),
                QTL = "caQTL")

df.gex.in_cs.mpra.max_prob =
  dplyr::filter(df.gex.in_cs.mpra, prob > 0.5) %>%
  dplyr::inner_join(df.loeuf.v4, by = c("region" = "gene_id")) %>%
  dplyr::group_by(rsid) %>%
  dplyr::filter(prob == max(prob)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(variant = stringr::str_replace_all(rsid, "_", ":"),
                QTL = "eQTL")



df.loeuf.emvar =
  dplyr::bind_rows(df.gex.in_cs.mpra.max_prob, df.atac.in_cs.mpra.max_prob) %>%
  dplyr::mutate(lof.oe_ci.upper_bin_decile = lof.oe_ci.upper_bin_decile %/% 2) %>%
  dplyr::inner_join(df.mpra_k562_meta, by = c("variant" = "variant_hg38")) %>%
  dplyr::select(QTL,
                variant,
                lof.oe_ci.upper_bin_decile,
                emVar,
                emVar_meta,
                log2Skew,
                log2Skew_meta) %>%
  tidyr::pivot_longer(-c(
    QTL,
    variant,
    lof.oe_ci.upper_bin_decile,
    log2Skew,
    log2Skew_meta
  )) %>%
  dplyr::mutate(
    cell_type = ifelse(name == "emVar", "K562", "Meta"),
    log2Skew = ifelse(cell_type == "K562", log2Skew, log2Skew_meta)
  ) %>%
  dplyr::group_split(lof.oe_ci.upper_bin_decile, cell_type, QTL) %>%
  purrr::map_dfr(function(data) {
    n_total <- nrow(data)
    n_emvar <- sum(data$value)
    
    # Return tibble with all results
    tibble::tibble(
      QTL = data$QTL[1],
      cell_type = data$cell_type[1],
      lof.oe_ci.upper_bin_decile = data$lof.oe_ci.upper_bin_decile[1],
      n_emvar = n_emvar,
      n_total = n_total,
      locusviz::binom_ci(n_emvar, n_total),
      locusviz::median_ci(abs(data$log2Skew))
    )
  })


plot_max_pip_bin_emvar = function(df, xlim = NULL) {
  pd = position_dodge(width = 0.9)
  ggplot(df, aes(frac_emvar, max_pip_bin, color = QTL)) +
    geom_hline(yintercept = 1.5,
               linetype = "dashed",
               color = "grey50") +
    geom_point(position = pd) +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   height = 0,
                   position = pd) +
    locusviz::get_default_theme(legend.position = c(1, 0.2),
                                legend.justification = c(1, 0.2)) +
    scale_x_continuous(expand = expansion(), labels = scales::label_percent()) +
    scale_y_discrete(expand = expansion()) +
    scale_color_manual(values = qtl.colors) +
    coord_cartesian(xlim = xlim) +
    labs(x = "% emVar", y = "Max PIP bin", color = "QTL")
}

plot_max_pip_bin_emvar_rho = function(df, xlim = NULL) {
  pd = position_dodge(width = 0.9)
  ggplot(df, aes(rho, max_pip_bin, color = QTL)) +
    geom_hline(yintercept = 1.5,
               linetype = "dashed",
               color = "grey50") +
    geom_point(position = pd) +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   height = 0,
                   position = pd) +
    locusviz::get_default_theme(
      hide.ylab = TRUE,
      legend.position = c(1, 0.2),
                                legend.justification = c(1, 0.2)) +
    scale_x_continuous(expand = expansion()) +
    scale_y_discrete(expand = expansion()) +
    scale_color_manual(values = qtl.colors) +
    coord_cartesian(xlim = xlim) +
    labs(x = expression(paste("Allelic effect correlation (", italic(rho), ")")), y = "Max PIP bin", color = "QTL")
}

plot_cascade_embar = function(df, xlim = NULL) {
  dplyr::mutate(
    df,
    y = interaction(
      scales::comma(n_total),
      scales::comma(n_emvar),
      qtl_mechanism_category
    ),
    y = factor(y, levels = y)
  ) %>%
    dplyr::arrange(cell_type, y) %>%
    ggplot(aes(frac_emvar, y, color = qtl_mechanism_category)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
    locusviz::get_default_theme(hide.ytitle = TRUE, legend.position = "none") +
    scale_x_continuous(labels = scales::label_percent(drop0trailing = TRUE),
                       expand = expansion()) +
    scale_color_manual(values = qtl_mechanism_colors) +
    guides(y = legendry::guide_axis_nested(
      key = legendry::key_range_auto(sep = "\\."),
      levels_text = list(
        element_text(color = "grey50"),
        element_text(color = "grey50"),
        element_text()
      )
    )) +
    labs(x = "% emVar") +
    coord_cartesian(xlim = xlim)
}

plot_loeuf_emvar = function(df, ylim = NULL) {
  dplyr::mutate(df, x = (lof.oe_ci.upper_bin_decile + 0.5) / 5) %>%
    ggplot(aes(x, frac)) +
    geom_ribbon(aes(
      ymin = frac_lower,
      ymax = frac_upper,
      fill = QTL
    ), alpha = 0.1) +
    geom_point(aes(color = QTL)) +
    geom_line(aes(color = QTL)) +
    coord_cartesian(xlim = c(0, 1), ylim = ylim) +
    scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion()) +
    scale_color_manual(values = qtl.colors) +
    scale_fill_manual(values = qtl.colors) +
    locusviz::get_default_theme(legend.position = c(1, 0),
                                legend.justification = c(1, 0)) +
    theme(plot.title = element_text(hjust = 0.02)) +
    labs(x = "LOEUF quintile", y = "% emVar", color = "QTL")
}

plot_loeuf_log2skew = function(df, ylim = NULL) {
  dplyr::mutate(df, x = (lof.oe_ci.upper_bin_decile + 0.5) / 5) %>%
    ggplot(aes(x, median)) +
    geom_ribbon(aes(
      ymin = median_lower,
      ymax = median_upper,
      fill = QTL
    ), alpha = 0.1) +
    geom_point(aes(color = QTL)) +
    geom_line(aes(color = QTL)) +
    coord_cartesian(xlim = c(0, 1), ylim = ylim) +
    scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
    scale_y_continuous(expand = expansion()) +
    scale_color_manual(values = qtl.colors) +
    scale_fill_manual(values = qtl.colors) +
    locusviz::get_default_theme(legend.position = c(1, 0),
                                legend.justification = c(1, 0)) +
    labs(x = "LOEUF quintile",
         y = expression(paste(
           "Abs. allelic effect (|", Delta, log[2], "FC", "|)"
         )),
         color = "QTL")
}

plt =
  list(
    dplyr::filter(df.max_pip_bin.emvar, cell_type == "K562") %>%
      plot_max_pip_bin_emvar(xlim = c(0, 0.23)) + labs(title = "K562"),
    dplyr::filter(df.max_pip_bin.emvar, cell_type == "Meta") %>%
      plot_max_pip_bin_emvar(xlim = c(0, 0.5)) + labs(title = "Meta"),
    dplyr::filter(df.max_pip_bin.log2skew_rho, cell_type == "K562" & emVar) %>%
      plot_max_pip_bin_emvar_rho(xlim = c(-0.1, 0.65)),
    dplyr::filter(df.max_pip_bin.log2skew_rho, cell_type == "Meta" & emVar) %>%
      plot_max_pip_bin_emvar_rho(xlim = c(-0.1, 0.65)),
    dplyr::filter(df.cascade.emvar, cell_type == "K562") %>%
      plot_cascade_embar(xlim = c(0, 0.43)),
    dplyr::filter(df.cascade.emvar, cell_type == "Meta") %>%
      plot_cascade_embar(xlim = c(0, 0.78))
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 2, byrow = FALSE) +
  patchwork::plot_annotation(tag_levels = "a")
plt


cowplot::save_plot(
  "figures/ExtendedDataFig6_mpra_replication.pdf",
  plt,
  base_height = 100,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/ExtendedDataFig6_mpra_replication.png",
  plt,
  base_height = 100,
  base_width = 180,
  units = "mm",
  dpi = 300
)

################################################################################
# Alasoo et al

fastqtl_cols <- c(
  "phenotype_id",
  "n_snps",
  "beta1",
  "beta2",
  "dummy",
  "snp_id",
  "distance",
  "p_nominal",
  "beta",
  "p_perm",
  "p_beta"
)
condition_map <- c(
  naive = "Naive",
  IFNg = "IFNg (18h)",
  SL1344 = "Salmonella (5h)",
  IFNg_SL1344 = "IFNg (18h) + Salmonella (5h)"
)

df.alasoo.atac.meta <- rgsutil::read_gsfile("gs://expansion_areas/multiome/misc/Alasoo_2018/ATAC_peak_metadata.txt.gz") %>%
  dplyr::mutate(peak_id_hg19 = stringr::str_c("chr", chr, "_", start, "_", end))

chain <- rtracklayer::import.chain(system.file("extdata", "hg19ToHg38.over.chain", package = "locusviz"))

gr.alasoo.hg19 <- GenomicRanges::makeGRangesFromDataFrame(
  df.alasoo.atac.meta,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
) %>%
  GenomeInfoDb::`seqlevelsStyle<-`("UCSC")

gr.alasoo.hg38 <- rtracklayer::liftOver(gr.alasoo.hg19, chain) %>%
  .[lengths(.) == 1] %>%
  unlist()
df.alasoo.hg38 <-
  GenomicRanges::as.data.frame(gr.alasoo.hg38) %>%
  dplyr::mutate(peak_id_hg38 = paste0(seqnames, "-", start, "-", end))

gr.fg <-
  dplyr::distinct(df.open4gene, peak_id) %>%
  tidyr::separate(
    peak_id,
    into = c("chrom", "start", "end"),
    sep = "-",
    remove = FALSE,
    convert = TRUE
  ) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)


gr.overlaps <- GenomicRanges::findOverlaps(gr.alasoo.hg38, gr.fg)

df.alasoo_to_fg <- tibble::tibble(alasoo_peak_id = df.alasoo.hg38$gene_id[S4Vectors::queryHits(gr.overlaps)],
                                  peak_id = gr.fg$peak_id[S4Vectors::subjectHits(gr.overlaps)])


df.alasoo.gex = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/misc/Alasoo_2018/RNA_FastQTL_*_500kb_permuted.txt.gz", function(df, path) {
  dplyr::mutate(df,
                qval = qvalue::qvalue(p_beta)$qvalue,
                condition = condition_map[[stringr::str_extract(path, "FastQTL_(.*)_500kb", group = 1)]])
}, col.names = fastqtl_cols) %>%
  dplyr::mutate(condition = factor(condition, levels = condition_map))

df.alasoo.atac = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/misc/Alasoo_2018/ATAC_FastQTL_*_50kb_cqn_perm.txt.gz", function(df, path) {
  dplyr::mutate(df,
                qval = qvalue::qvalue(p_beta)$qvalue,
                condition = condition_map[[stringr::str_extract(path, "FastQTL_(.*)_50kb", group = 1)]])
}, col.names = fastqtl_cols)  %>%
  dplyr::mutate(condition = factor(condition, levels = condition_map)) %>%
  dplyr::inner_join(
    df.alasoo_to_fg,
    by = c("phenotype_id" = "alasoo_peak_id"),
    relationship = "many-to-many"
  )

df.alasoo.gex.loeuf =
  dplyr::inner_join(df.alasoo.gex, df.loeuf.v4, by = c("phenotype_id" = "gene_id")) %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile) %>%
  dplyr::mutate(lof.oe_ci.upper_bin_decile = lof.oe_ci.upper_bin_decile %/% 2) %>%
  dplyr::group_by(lof.oe_ci.upper_bin_decile, condition) %>%
  dplyr::summarize(locusviz::median_ci(abs(beta)), .groups = "drop")

df.alasoo.atac.max_abs_beta = 
  dplyr::inner_join(
    df.alasoo.atac,
    dplyr::distinct(df.open4gene.sig, peak_id, gene_id),
    by = "peak_id",
    relationship = "many-to-many"
  ) %>%
  dplyr::group_by(gene_id, condition) %>%
  dplyr::summarize(max_abs_beta = max(abs(beta)), .groups = "drop")

df.alasoo.atac.loeuf =
  dplyr::inner_join(df.alasoo.atac.max_abs_beta, df.loeuf.v4, by = "gene_id") %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile) %>%
  dplyr::mutate(lof.oe_ci.upper_bin_decile = lof.oe_ci.upper_bin_decile %/% 2) %>%
  dplyr::group_by(lof.oe_ci.upper_bin_decile, condition) %>%
  dplyr::summarize(locusviz::median_ci(max_abs_beta), .groups = "drop")

df.alasoo.gex.loeuf.rho =
  dplyr::inner_join(df.alasoo.gex, df.loeuf.v4, by = c("phenotype_id" = "gene_id")) %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile) %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(locusviz::spearman_ci(abs(beta), lof.oe_ci.upper), .groups = "drop")

df.alasoo.atac.loeuf.rho =
  dplyr::inner_join(df.alasoo.atac.max_abs_beta, df.loeuf.v4, by = "gene_id") %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile) %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(locusviz::spearman_ci(max_abs_beta, lof.oe_ci.upper), .groups = "drop")

condition.colors <- c(
  "Naive" = "grey50",
  "IFNg (18h)" = BuenColors::jdb_palette("corona")[3],
  "Salmonella (5h)" = BuenColors::jdb_palette("corona")[4],
  "IFNg (18h) + Salmonella (5h)" = BuenColors::jdb_palette("corona")[5]
)

p.alasoo.gex.loeuf=
  dplyr::mutate(df.alasoo.gex.loeuf, x = (lof.oe_ci.upper_bin_decile + 0.5) / 5) %>%
  ggplot(aes(x, median)) +
  geom_ribbon(aes(ymin = median_lower, ymax = median_upper, fill = condition),
              alpha = 0.1) +
  geom_point(aes(color = condition)) +
  geom_line(aes(color = condition)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  scale_color_manual(values = condition.colors) +
  scale_fill_manual(values = condition.colors, guide = "none") +
  locusviz::get_default_theme(legend.position = c(1, 0),
                              legend.justification = c(1, 0)) +
  theme(plot.title = element_text(hjust = 0.02)) +
  labs(x = "LOEUF quintile",
       y = expression(paste("Median |", italic(beta)[eQTL], "|")),
       color = "Stimulation",
       title = "Alasoo, et al. (2018)") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.3))

p.alasoo.atac.loeuf=
  dplyr::mutate(df.alasoo.atac.loeuf, x = (lof.oe_ci.upper_bin_decile + 0.5) / 5) %>%
  ggplot(aes(x, median)) +
  geom_ribbon(aes(ymin = median_lower, ymax = median_upper, fill = condition),
              alpha = 0.1) +
  geom_point(aes(color = condition)) +
  geom_line(aes(color = condition)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  scale_color_manual(values = condition.colors) +
  scale_fill_manual(values = condition.colors, guide = "none") +
  locusviz::get_default_theme(legend.position = "none") +
  labs(x = "LOEUF quintile",
       y = expression(paste("Median |", italic(beta)[caQTL], "|")),
       color = "Stimulation") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.85))



p.mpra.alasoo = list(
  p.alasoo.gex.loeuf + labs(tag = "a"),
  p.alasoo.atac.loeuf + labs(tag = "b"),
  dplyr::filter(df.loeuf.emvar, cell_type == "K562") %>%
    plot_loeuf_emvar(ylim = c(0, 0.3)) + labs(title = "K562", tag = "c"),
  dplyr::filter(df.loeuf.emvar, cell_type == "K562") %>%
    plot_loeuf_log2skew(ylim = c(0, 0.32)) + labs(tag = "e"),
  dplyr::filter(df.loeuf.emvar, cell_type == "Meta") %>%
    plot_loeuf_emvar(ylim = c(0, 0.6)) + labs(title = "Meta", tag = "d"),
  dplyr::filter(df.loeuf.emvar, cell_type == "Meta") %>%
    plot_loeuf_log2skew(ylim = c(0, 0.22)) + labs(tag = "f")
)%>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 2)
p.mpra.alasoo

cowplot::save_plot(
  "figures/ExtendedDataFig10_mpra_alasoo_loeuf.pdf",
  p.mpra.alasoo,
  base_height = 170,
  base_width = 120,
  units = "mm"
)

cowplot::save_plot(
  "figures/ExtendedDataFig10_mpra_alasoo_loeuf.png",
  p.mpra.alasoo,
  base_height = 170,
  base_width = 120,
  units = "mm",
  dpi = 300
)
