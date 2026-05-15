library(dplyr)
library(ggplot2)
source(here::here("R/const.R"))

df.gex.obs = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/h5ad/integrated_gex_batch1_5.fgid.qc.obs.tsv.gz"
) %>%
  dplyr::mutate(
    multiplexed = stringr::str_detect(pool_name, "^pool"),
    multiplexed = factor(multiplexed, levels = c("TRUE", "FALSE"))
  )

df.atac.obs = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/h5ad/per_chrom/integrated_atac_batch1_5.chr1.fgid.obs.tsv.gz"
) %>%
  dplyr::mutate(
    predicted.celltype.l1 = stringr::str_replace_all(predicted.celltype.l1, " ", "_"),
    predicted.celltype.l2 = stringr::str_replace_all(predicted.celltype.l2, " ", "_")
  ) %>%
  dplyr::mutate(
    multiplexed = stringr::str_detect(pool_name, "^pool"),
    multiplexed = factor(multiplexed, levels = c("TRUE", "FALSE"))
  )

df.gex.var = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/h5ad/integrated_gex_batch1_5.fgid.qc.var.tsv.gz"
)

df.gex.n_samples <- rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/h5ad/integrated_gex_batch1_5.fgid.qc.n_samples.txt"
)
df.atac.n_samples <- rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/h5ad/integrated_atac_batch1_5.fgid.n_samples.txt"
)

if (!file.exists("tables/ST1_multiome_overview.tsv")) {
  df.obs.out =
    dplyr::bind_rows(
      dplyr::mutate(df.gex.obs, Atlas = "snRNA"),
      dplyr::mutate(df.atac.obs, Atlas = "snATAC")
    ) %>%
    dplyr::group_by(FINNGENID) %>%
    dplyr::summarize(
      snrna_n_qc_nuclei = sum(Atlas == "snRNA"),
      snatac_n_qc_nuclei = sum(Atlas == "snATAC"),
      snmultiome_n_qc_nuclei = length(intersect(barcode[Atlas == "snRNA"], barcode[Atlas == "snATAC"])),
      snrna_median_total_counts = median(total_counts[Atlas == "snRNA"]),
      snrna_median_total_counts_wo_highly_expressed = median(total_counts_wo_highly_expressed[Atlas == "snRNA"]),
      snrna_median_n_genes = median(n_genes[Atlas == "snRNA"]),
      snrna_median_pct_counts_mito = median(pct_counts_mito[Atlas == "snRNA"]),
      snatac_median_total_counts = median(total_counts[Atlas == "snATAC"]),
      snatac_median_total_counts_wo_highly_expressed = median(total_counts_wo_highly_expressed[Atlas == "snATAC"]),
      snatac_median_n_peaks = median(n_genes[Atlas == "snATAC"]),
      snatac_median_TSS_enrichment = median(TSS.enrichment[Atlas == "snATAC"]),
      multiplexed = any(stringr::str_detect(pool_name, "^pool"))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(FINNGENID) %>%
    dplyr::rename(ID = FINNGENID) %>%
    dplyr::mutate(ID = seq_len(n()))
  
  export_table(df.obs.out,
               "tables/ST1_multiome_overview.tsv",
               save_googlesheet = FALSE)
}


# GEX OBS smaples
length(unique(df.gex.obs$FINNGENID))
# 1108
max(df.gex.n_samples$n_samples)
# 1103

# ATAC OBS samples
length(unique(df.atac.obs$FINNGENID))
# 1093
max(df.atac.n_samples$n_samples)
# 1088

dplyr::summarize(
  df.gex.obs,
  total_n_nuclei = n(),
  m_n_genes = median(n_genes),
  m_n_counts = median(n_counts)
)
#   total_n_nuclei m_n_genes m_n_counts
# 1       10612211      1056       1828

dplyr::group_by(df.gex.obs, multiplexed) %>%
  dplyr::summarize(
    total_n_nuclei = n(),
    m_n_genes = median(n_genes),
    m_n_counts = median(n_counts)
  )
# multiplexed total_n_nuclei m_n_genes m_n_counts
# <fct>                <int>     <int>      <int>
# 1 TRUE               7481200      1053       1824
# 2 FALSE              3131011      1070       1846

dplyr::group_by(df.gex.obs, FINNGENID) %>%
  dplyr::summarize(n_nuclei = n()) %>%
  dplyr::ungroup() %>%
  dplyr::summarize(m_n_nuclei = median(n_nuclei))
# 9506

dplyr::group_by(df.gex.obs, FINNGENID, multiplexed) %>%
  dplyr::summarize(n_nuclei = n()) %>%
  dplyr::group_by(multiplexed) %>%
  dplyr::summarize(m_n_nuclei = median(n_nuclei))
# multiplexed m_n_nuclei
# <fct>            <int>
#   1 TRUE              9571
# 2 FALSE             8933

dplyr::summarize(
  df.atac.obs,
  total_n_nuclei = n(),
  m_n_genes = median(n_genes),
  m_n_counts = median(n_counts),
  m_TSS = median(TSS.enrichment)
)
# total_n_nuclei m_n_genes m_n_counts    m_TSS
# 1        7100899      4174       4745 6.466506

dplyr::group_by(df.atac.obs, multiplexed) %>%
  dplyr::summarize(
    total_n_nuclei = n(),
    m_n_genes = median(n_genes),
    m_n_counts = median(n_counts),
    m_TSS = median(TSS.enrichment)
  )
# multiplexed total_n_nuclei m_n_genes m_n_counts m_TSS
# <fct>                <int>     <dbl>      <dbl> <dbl>
# 1 TRUE               5064006      4009       4543  6.52
# 2 FALSE              2036893      4742       5455  6.34

dplyr::group_by(df.atac.obs, FINNGENID) %>%
  dplyr::summarize(n_nuclei = n()) %>%
  dplyr::ungroup() %>%
  dplyr::summarize(m_n_nuclei = median(n_nuclei))
# 6,692

dplyr::group_by(df.atac.obs, FINNGENID, multiplexed) %>%
  dplyr::summarize(n_nuclei = n()) %>%
  dplyr::group_by(multiplexed) %>%
  dplyr::summarize(m_n_nuclei = median(n_nuclei))
# multiplexed m_n_nuclei
# <fct>            <dbl>
#   1 TRUE              6731
# 2 FALSE             6432

df.n_nuclei =
  dplyr::bind_rows(
    dplyr::group_by(df.gex.obs, FINNGENID) %>%
      dplyr::summarize(n_nuclei = n(), Atlas = "snRNA-Seq"),
    dplyr::group_by(df.atac.obs, FINNGENID) %>%
      dplyr::summarize(n_nuclei = n(), Atlas = "snATAC-Seq")
  )

p.n_nuclei =
  ggplot(df.n_nuclei, aes(n_nuclei, fill = Atlas)) +
  geom_histogram(bins = 50,
                 position = "identity",
                 alpha = 0.9) +
  geom_vline(
    aes(xintercept = m),
    linetype = "dashed",
    color = "grey50",
    data = dplyr::group_by(df.n_nuclei, Atlas) %>%
      dplyr::summarize(m = median(n_nuclei))
  ) +
  scale_fill_manual(values = BuenColors::jdb_palette("corona")) +
  scale_x_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  locusviz::get_default_theme() +
  theme(
    legend.position.inside = c(1, 1),
    legend.justification.inside = c(1, 1),
    legend.title = element_blank()
  ) +
  labs(x = "# nuclei", y = "# individuals")
p.n_nuclei

p.gex.n_genes =
  ggplot(df.gex.obs, aes(n_genes, fill = multiplexed)) +
  geom_histogram(bins = 50,
                 position = "identity",
                 alpha = 0.9) +
  geom_vline(
    aes(xintercept = m),
    linetype = "dashed",
    color = "grey50",
    data = dplyr::summarize(df.gex.obs, m = median(n_genes))
  ) +
  scale_x_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_fill_manual(values = as.character(locusviz::distinct_shades(qtl.colors["eQTL"], 2))) +
  locusviz::get_default_theme() +
  labs(x = "# genes", y = "# nuclei", fill = "Multiplexed")
p.gex.n_genes

p.gex.n_counts =
  ggplot(df.gex.obs, aes(n_counts, fill = multiplexed)) +
  geom_histogram(bins = 50,
                 position = "identity",
                 alpha = 0.9) +
  geom_vline(
    aes(xintercept = m),
    linetype = "dashed",
    color = "grey50",
    data = dplyr::summarize(df.gex.obs, m = median(n_counts))
  ) +
  scale_x_log10(labels = scales::label_comma(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_fill_manual(values = as.character(locusviz::distinct_shades(qtl.colors["eQTL"], 2))) +
  locusviz::get_default_theme() +
  labs(x = "# UMI", y = "# nuclei", fill = "Multiplexed")
p.gex.n_counts

p.gex.pct_counts_mito =
  ggplot(df.gex.obs, aes(pct_counts_mito / 100, fill = multiplexed)) +
  geom_histogram(bins = 50,
                 position = "identity",
                 alpha = 0.9) +
  geom_vline(
    aes(xintercept = m),
    linetype = "dashed",
    color = "grey50",
    data = dplyr::reframe(df.gex.obs, m = c(median(
      pct_counts_mito / 100
    ), 0.2))
  ) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_fill_manual(values = as.character(locusviz::distinct_shades(qtl.colors["eQTL"], 2))) +
  locusviz::get_default_theme() +
  labs(x = "% mito", y = "# nuclei", fill = "Multiplexed")
p.gex.pct_counts_mito

p.gex.var.mean_counts =
  ggplot(df.gex.var, aes(mean_counts)) +
  geom_histogram(bins = 50, fill = qtl.colors["caQTL"]) +
  geom_vline(
    aes(xintercept = m),
    linetype = "dashed",
    color = "grey50",
    data = dplyr::summarize(df.gex.var, m = median(mean_counts))
  ) +
  scale_x_log10(labels = scales::label_comma(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion()) +
  locusviz::get_default_theme() +
  labs(x = "Mean # UMI", y = "# genes")
p.gex.var.mean_counts

p.atac.n_counts =
  ggplot(df.atac.obs, aes(n_counts, fill = multiplexed)) +
  geom_histogram(bins = 50,
                 position = "identity",
                 alpha = 0.9) +
  geom_vline(
    aes(xintercept = m),
    linetype = "dashed",
    color = "grey50",
    data = dplyr::summarize(df.atac.obs, m = median(n_counts))
  ) +
  scale_x_log10(labels = scales::label_comma(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_fill_manual(values = as.character(locusviz::distinct_shades(qtl.colors["caQTL"], 2))) +
  locusviz::get_default_theme() +
  labs(x = "# fragments", y = "# nuclei", fill = "Multiplexed")
p.atac.n_counts

p.atac.tss =
  ggplot(df.atac.obs, aes(TSS.enrichment, fill = multiplexed)) +
  geom_histogram(bins = 50,
                 position = "identity",
                 alpha = 0.9) +
  geom_vline(
    aes(xintercept = m),
    linetype = "dashed",
    color = "grey50",
    data = dplyr::summarize(df.atac.obs, m = median(TSS.enrichment))
  ) +
  scale_x_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_fill_manual(values = as.character(locusviz::distinct_shades(qtl.colors["caQTL"], 2))) +
  locusviz::get_default_theme() +
  labs(x = "TSS enrichment", y = "# nuclei", fill = "Multiplexed")
p.atac.tss

calculate_celltype_freq <- function(df, cell_type_col, atlas_name, level_name) {
  dplyr::group_by(df, {{cell_type_col}}) %>%
    dplyr::summarize(n_nuclei = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      frac = n_nuclei / sum(n_nuclei),
      Atlas = atlas_name,
      cell_type = sprintf("predicted.celltype.%s.%s", level_name, {{cell_type_col}})
    ) %>%
    dplyr::select(Atlas, cell_type, tidyselect::everything()) %>%
    dplyr::select(-({{cell_type_col}}))
}

df.n_nuclei.cell_type = tidyr::crossing(tibble::tibble(
  df = list(df.gex.obs, df.atac.obs),
  atlas = c("snRNA-Seq", "snATAC-Seq")
),
tibble::tibble(
  col = c("predicted.celltype.l1", "predicted.celltype.l2"),
  level = c("l1", "l2")
)) %>%
  purrr::pmap_dfr(function(df, col, atlas, level) {
    calculate_celltype_freq(df, !!rlang::sym(col), atlas, level)
  })

if (!file.exists("tables/ST2_cell_type_freqs.tsv")) {
  export_table(df.n_nuclei.cell_type,
               "tables/ST2_cell_type_freqs.tsv",
               "ST2")
}

p.n_nuclei.cell_type =
  dplyr::bind_rows(
    dplyr::group_by(df.gex.obs, predicted.celltype.l2) %>%
      dplyr::summarize(n_nuclei = n(), Atlas = "snRNA-Seq") %>%
      dplyr::ungroup() %>%
      dplyr::rename(cell_type = predicted.celltype.l2),
    dplyr::group_by(df.atac.obs, predicted.celltype.l2) %>%
      dplyr::summarize(n_nuclei = n(), Atlas = "snATAC-Seq") %>%
      dplyr::ungroup() %>%
      dplyr::rename(cell_type = predicted.celltype.l2),
  ) %>%
  dplyr::filter(cell_type %in% remove_cell_type_prefix(l2.cell_types)) %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = cell_type[order(n_nuclei[which(Atlas == "snRNA-Seq")])])) %>%
  ggplot(aes(n_nuclei, cell_type, fill = Atlas)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = BuenColors::jdb_palette("corona")) +
  scale_x_continuous(expand = expansion(),
                     labels = scales::label_comma(scale = 1e-6)) +
  scale_y_discrete(labels = l2.labels) +
  locusviz::get_default_theme(legend.position = "none", hide.ytitle = TRUE) +
  labs(x = "# nuclei (M)")
p.n_nuclei.cell_type

################################################################################
df.consensus = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/atac/peaks/integrated_atac_batch1_5.iterative_merged.bed.gz",
  header = FALSE
) %>%
  dplyr::transmute(
    chrom = V1,
    start = V2,
    end = V3,
    width = end - start,
    merged_peak_names = stringr::str_remove_all(V4, "predicted\\.celltype\\.")
  )

if (!file.exists("./tables/ST3_consensus_peaks.tsv")) {
  dplyr::select(df.consensus, -merged_peak_names) %>%
    export_table("./tables/ST3_consensus_peaks.tsv", save_googlesheet = FALSE)
}

p.peak.width =
  ggplot(df.consensus, aes(width)) +
  geom_histogram(bins = 50, fill = qtl.colors["caQTL"]) +
  scale_x_log10(labels = scales::label_comma(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion()) +
  locusviz::get_default_theme() +
  labs(x = "Peak width", y = "# peaks")

p.peak.width

################################################################################
# open4gene

if (!file.exists("data/open4gene.rds")) {
  # rgsutil::read_gsfiles("gs://expansion_areas/multiome/batch1_5/multilinker/open4gene_offset_LSI/open4gene.*.results.tsv.gz", combine = "rows") %>%
  df.open4gene =
    rgsutil::read_gsfiles(
      "gs://expansion_areas/multiome/batch1_5/multilinker/open4gene_score/open4gene.*.results.tsv.gz",
      combine = "rows"
    ) %>%
    tidyr::drop_na(hurdle_zero_nlog10p, hurdle_count_nlog10p) %>%
    dplyr::filter(expr_cell_num >= 5000) %>%
    dplyr::left_join(df.features, by = c("gene_id" = "phenotype_id")) %>%
    dplyr::mutate(parse_peak_id(peak_id)) %>%
    dplyr::mutate(
      distance = (peak_start + peak_end) %/% 2 - tss,
      distance = ifelse(strand == "+", distance, -distance),
      nonzero_atac_prop = open_cell_num / total_cell_num,
      gene_global_start = locusviz::get_global_position(chrom, start, reference_genome = "GRCh38"),
      gene_global_end = locusviz::get_global_position(chrom, end, reference_genome = "GRCh38"),
      gene_global_pos = locusviz::get_global_position(chrom, (start + end) %/% 2, reference_genome = "GRCh38"),
      peak_global_start = locusviz::get_global_position(peak_chrom, peak_start, reference_genome = "GRCh38"),
      peak_global_end = locusviz::get_global_position(peak_chrom, peak_end, reference_genome = "GRCh38"),
      peak_global_pos = locusviz::get_global_position(peak_chrom, (peak_start + peak_end) %/% 2, reference_genome = "GRCh38")
    ) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::mutate(
      fasthurdle::joint_score_test(
        10 ** -hurdle_zero_nlog10p,
        10 ** -hurdle_count_nlog10p,
        alpha = 0.05
      )
    ) %>%
    dplyr::ungroup()
  
  saveRDS(df.open4gene, "data/open4gene.rds")
  saveRDS(dplyr::filter(df.open4gene, sig), "data/open4gene.sig.rds")
}

# df.open4gene = readRDS("data/open4gene.rds")
df.open4gene.sig = readRDS("data/open4gene.sig.rds")
link_peak_id = unique(df.open4gene.sig$peak_id)


export_peak_bed = function(df, path, overwrite = FALSE) {
  dplyr::distinct(df, peak_chrom, peak_start, peak_end) %>%
    dplyr::mutate(peak_chrom_num = ifelse(peak_chrom == "chrX", 23, as.numeric(
      stringr::str_remove(peak_chrom, "^chr")
    ))) %>%
    dplyr::arrange(peak_chrom_num, peak_start) %>%
    dplyr::select(-peak_chrom_num) %>%
    rgsutil::write_gsfile(path, col.names = FALSE, overwrite = overwrite)
}

export_peak_bed(
  df.open4gene.sig,
  "gs://expansion_areas/multiome/batch1_5/multilinker/open4gene_score/open4gene.sig.peak.bed",
  overwrite = FALSE
)


df.open4gene.sig.peak.mode =
  dplyr::group_by(df.open4gene.sig, peak_id, peak_chrom, peak_start, peak_end) %>%
  dplyr::summarize(
    n_mechanisms = length(unique(mode)),
    n_genes = length(unique(gene_id)),
    n_cell_types = length(unique(cell_type)),
    combined_modes = stringr::str_c(sort(unique(mode)), collapse = ","),
    mode = dplyr::case_when(
      any(mode == "dual") ~ "switch_and_rheostat",
      combined_modes %in%  c("rheostat", "rheostat,omnibus_only") ~ "rheostat",
      combined_modes %in% c("switch", "switch,omnibus_only") ~ "switch",
      combined_modes %in% c("switch,rheostat", "switch,rheostat,omnibus_only") ~ "switch_or_rheostat",
      TRUE ~ combined_modes
    ),
    .groups = "drop"
  )

if (!file.exists("./tables/ST5_peak_regulatory_mode.tsv")) {
  export_table(
    df.open4gene.sig.peak.mode,
    "./tables/ST5_peak_regulatory_mode.tsv",
    save_googlesheet = FALSE
  )
}

table(df.open4gene.sig.peak.mode$mode)
# omnibus_only            rheostat              switch switch_and_rheostat  switch_or_rheostat
#           26               16194                5418               78789                4353
table(df.open4gene.sig.peak.mode$mode) / nrow(df.open4gene.sig.peak.mode)
# omnibus_only            rheostat              switch switch_and_rheostat  switch_or_rheostat
#  0.000248139         0.154552395         0.051708341         0.751946936         0.041544188


dplyr::group_split(df.open4gene.sig.peak.mode, mode) %>%
  purrr::map(function(data) {
    mode = data$mode[1]
    export_peak_bed(
      data,
      sprintf(
        "gs://expansion_areas/multiome/batch1_5/multilinker/open4gene_score/open4gene.sig.%s.peak.bed",
        mode
      ),
      overwrite = FALSE
    )
  })

dplyr::distinct(df.open4gene.sig, peak_id, gene_id) %>%
  nrow()
# 593765
length(link_peak_id)
# 104780
length(link_peak_id) / nrow(df.consensus)
# 0.3513466
length(unique(df.open4gene.sig$gene_id))
# 12018

df.open4gene.n_genes =
  dplyr::group_by(df.open4gene.sig, peak_id) %>%
  dplyr::summarize(n_genes = length(unique(gene_id)), .groups = "drop")
df.open4gene.n_peaks =
  dplyr::group_by(df.open4gene.sig, gene_id) %>%
  dplyr::summarize(n_peaks = length(unique(peak_id)), .groups = "drop")

nrow(df.open4gene.sig)
# 1095949
dplyr::summarize(df.open4gene.n_genes,
                 median = median(n_genes),
                 max = max(n_genes))
# median   max
# <int> <int>
#   1      4    46

table(df.open4gene.sig$mode) / nrow(df.open4gene.sig)
#        dual          switch        rheostat    omnibus_only not_significant
# 0.419405465     0.216726326     0.356812224     0.007055985     0.000000000

dplyr::group_by(df.open4gene.sig, mode) %>%
  dplyr::summarize(median(abs(distance)))
# mode         `median(abs(distance))`
# <fct>                          <dbl>
#   1 dual                          254006
# 2 switch                        166827
# 3 rheostat                      444210
# 4 omnibus_only                  448073

dplyr::summarize(df.open4gene.n_peaks,
                 median = median(n_peaks),
                 max = max(n_peaks))
# median   max
# <int> <int>
#   1     36   395


dplyr::group_by(df.open4gene.sig, peak_id) %>%
  dplyr::summarize(n_genes = length(unique(gene_id))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_genes == max(n_genes))
# peak_id                 n_genes
# <chr>                     <int>
#   1 chr11-65150307-65151164      46
# 2 chr6-32437315-32437814       46

p.open4gene.n_genes =
  ggplot(df.open4gene.n_genes, aes(n_genes)) +
  geom_histogram(
    bins = max(df.open4gene.n_genes$n_genes),
    fill = BuenColors::jdb_palette("corona")[3]
  ) +
  geom_vline(
    aes(xintercept = m),
    linetype = "dashed",
    color = "grey50",
    data = dplyr::summarize(df.open4gene.n_genes, m = median(n_genes))
  ) +
  scale_x_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion()) +
  locusviz::get_default_theme() +
  labs(x = "# linked genes", y = "# peaks")
p.open4gene.n_genes

p.open4gene.n_peaks =
  ggplot(df.open4gene.n_peaks, aes(n_peaks)) +
  geom_histogram(bins = 50, fill = BuenColors::jdb_palette("corona")[4]) +
  geom_vline(
    aes(xintercept = m),
    linetype = "dashed",
    color = "grey50",
    data = dplyr::summarize(df.open4gene.n_peaks, m = median(n_peaks))
  ) +
  scale_x_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion()) +
  locusviz::get_default_theme() +
  labs(x = "# linked peaks", y = "# genes")

df.open4gene.sig.peaks.out = dplyr::group_by(df.open4gene.sig, peak_id) %>%
  dplyr::summarize(
    n_genes = length(unique(gene_id)),
    genes = stringr::str_c(sort(unique(symbol)), collapse = ","),
    sig_l1_cell_types = stringr::str_c(remove_cell_type_prefix(sort(unique(
      cell_type[filter_l1_cell_types(cell_type)]
    ))), collapse = ","),
    sig_l2_cell_types = stringr::str_c(remove_cell_type_prefix(sort(unique(
      cell_type[!filter_l1_cell_types(cell_type)]
    ))), collapse = ","),
    .groups = "drop"
  )

df.open4gene.sig.genes.out = dplyr::group_by(df.open4gene.sig, gene_id, symbol) %>%
  dplyr::summarize(
    n_sig_peaks = length(unique(peak_id)),
    genes = stringr::str_c(sort(unique(symbol)), collapse = ","),
    sig_l1_cell_types = stringr::str_c(remove_cell_type_prefix(sort(unique(
      cell_type[filter_l1_cell_types(cell_type)]
    ))), collapse = ","),
    sig_l2_cell_types = stringr::str_c(remove_cell_type_prefix(sort(unique(
      cell_type[!filter_l1_cell_types(cell_type)]
    ))), collapse = ","),
    .groups = "drop"
  )

if (!file.exists("tables/ST4_open4gene_sig.tsv")) {
  export_table(df.open4gene.sig.genes.out,
               "tables/ST4_open4gene_sig.tsv",
               save_googlesheet = FALSE)
}

p.link.mode.density =
  dplyr::filter(df.open4gene.sig, mode %in% c("rheostat", "dual", "switch")) %>%
  dplyr::mutate(mode = factor(mode, levels = link_mode_levels)) %>%
  ggplot(aes(distance, color = mode)) +
  geom_line(stat = "density") +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  scale_color_manual(values = link_mode_colors, labels = stringr::str_to_title) +
  scale_x_continuous(
    labels = scales::label_number(scale = 1e-6, drop0trailing = TRUE),
    expand = expansion()
  ) +
  scale_y_continuous(
    expand = expansion(),
    labels = scales::label_number(scale = 1e6, drop0trailing = TRUE)
  ) +
  coord_cartesian(xlim = c(-1e6, 1e6)) +
  locusviz::get_default_theme(hide.xlab = TRUE) +
  labs(y = expression("Density (" * "\u00D7" * 10^-6 * ")"),
       color = "Mode")
p.link.mode.density

df.link.mode.fraction <-
  dplyr::filter(df.open4gene.sig, mode %in% c("rheostat", "dual", "switch")) %>%
  dplyr::mutate(mode = factor(mode, levels = link_mode_levels)) %>%
  dplyr::mutate(distance_bin = cut(
    distance,
    breaks = seq(-1e6, 1e6, by = 20000),
    include.lowest = TRUE
  )) %>%
  dplyr::group_split(distance_bin) %>%
  purrr::map_dfr(function(data) {
    n_switch = sum(data$mode == "switch")
    n_rheostat = sum(data$mode == "rheostat")
    n_dual = sum(data$mode == "dual")
    n_total = nrow(data)
    
    tibble::tibble(
      distance = mean(data$distance),
      locusviz::binom_ci(n_switch, n_total, colname = "frac_switch"),
      locusviz::binom_ci(n_rheostat, n_total, colname = "frac_rheostat"),
      locusviz::binom_ci(n_dual, n_total, colname = "frac_dual")
    )
  }) %>%
  dplyr::rename_with(~ sub("(frac_(switch|rheostat|dual))$", "\\1_value", .x)) %>%
  tidyr::pivot_longer(
    cols = starts_with("frac_"),
    names_to = c("mode", ".value"),
    names_pattern = "frac_(switch|rheostat|dual)_(value|lower|upper)"
  ) %>%
  dplyr::mutate(mode = factor(mode, levels = link_mode_levels))


p.link.mode.fraction =
  ggplot(df.link.mode.fraction, aes(distance, value)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = mode), alpha = 0.1) +
  geom_line(aes(color = mode)) +
  scale_color_manual(values = link_mode_colors) +
  scale_fill_manual(values = link_mode_colors, guide = "none") +
  locusviz::get_default_theme(legend.position = "none") +
  scale_x_continuous(
    labels = scales::label_number(scale = 1e-6, drop0trailing = TRUE),
    expand = expansion()
  ) +
  scale_y_continuous(expand = expansion(), labels = scales::label_percent()) +
  coord_cartesian(xlim = c(-1e6, 1e6), ylim = c(0, 0.6)) +
  labs(x = "Distance to TSS (Mb)", y = "Fraction")
p.link.mode.fraction


############
# Alternating modes
df.link.switch_or_rheostat =
  dplyr::filter(df.open4gene.sig.peak.mode, mode == "switch_or_rheostat") %>%
  dplyr::select(peak_id) %>%
  dplyr::inner_join(df.open4gene.sig)

# IRF4, ENSG00000137265
dplyr::filter(df.link.switch_or_rheostat, gene_id == "ENSG00000137265") %>%
  dplyr::select(
    peak_id:cell_type,
    mode,
    hurdle_zero_beta,
    hurdle_zero_nlog10p,
    hurdle_count_beta,
    hurdle_count_nlog10p,
    distance,
    p_joint
  ) %>%
  View()

# TCF4, ENSG00000196628
dplyr::filter(df.link.switch_or_rheostat, gene_id == "ENSG00000196628") %>%
  dplyr::select(
    peak_id:cell_type,
    mode,
    hurdle_zero_beta,
    hurdle_zero_nlog10p,
    hurdle_count_beta,
    hurdle_count_nlog10p,
    distance,
    p_joint
  ) %>%
  View()

# CCDC50, ENSG00000152492
dplyr::filter(df.link.switch_or_rheostat, gene_id == "ENSG00000152492") %>%
  dplyr::select(
    peak_id:cell_type,
    mode,
    hurdle_zero_beta,
    hurdle_zero_nlog10p,
    hurdle_count_beta,
    hurdle_count_nlog10p,
    distance,
    p_joint
  ) %>%
  View()

plot_alt_mode_forest = function(df.open4gene.sig, gene_id, peak_id, title, legend.position = "none") {
  pd = position_dodge(width = 0.9)
  df = dplyr::filter(
    df.open4gene.sig,
    gene_id == .env$gene_id &
      peak_id == .env$peak_id &
      cell_type %in% c(
        "predicted.celltype.l1.B",
        "predicted.celltype.l1.Mono",
        "predicted.celltype.l1.DC"
      )
  ) %>%
    dplyr::select(peak_id,
                  gene_id,
                  cell_type,
                  tidyselect::starts_with("hurdle_")) %>%
    tidyr::pivot_longer(
      cols = tidyselect::starts_with("hurdle_"),
      names_to = c("model", ".value"),
      names_pattern = "hurdle_(zero|count)_(beta|se|nlog10p|spa)"
    ) %>%
    dplyr::mutate(
      cell_type = remove_cell_type_prefix(cell_type),
      cell_type = factor(cell_type, levels = c("B", "Mono", "DC")),
      model = factor(stringr::str_to_title(model), levels = rev(c("Zero", "Count")))
    )
  
  ggplot(df, aes(beta, cell_type, color = model)) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               color = "grey50") +
    geom_errorbar(
      aes(
        xmin = beta - qnorm(0.05 / 2, lower.tail = F) * se,
        xmax = beta + qnorm(0.05 / 2, lower.tail = F) * se
      ),
      width = 0,
      position = pd
    ) +
    geom_point(position = pd) +
    locusviz::get_default_theme(
      legend.position = legend.position,
      legend.justification = legend.position,
      hide.ytitle = TRUE
    ) +
    scale_color_manual(values = model.cols) +
    guides(color = guide_legend(reverse = TRUE)) +
    labs(x = sprintf("Beta (%s)", stringr::str_replace(peak_id, "-", ":")), color = "Model", title = title)
}

p.IRF4 = plot_alt_mode_forest(df.open4gene.sig,
                              "ENSG00000137265",
                              "chr6-246344-246843",
                              "IRF4",
                              legend.position = c(0, 0.9))
p.TCF4 = plot_alt_mode_forest(df.open4gene.sig,
                              "ENSG00000196628",
                              "chr18-55235635-55236760",
                              "TCF4")
p.CCDC50 = plot_alt_mode_forest(df.open4gene.sig,
                                "ENSG00000152492",
                                "chr3-192063473-192063972",
                                "CCDC50")

################################################################################
# ZNF_19p13.11
a =
  dplyr::left_join(df.open4gene.n_peaks,
                   df.features,
                   by = c("gene_id" = "phenotype_id")) %>%
  dplyr::arrange(dplyr::desc(n_peaks)) %>%
  dplyr::filter(chrom == "chr19" &
                  start > 18e6 & end <= 19e6)

b = dplyr::filter(df.open4gene.sig, gene_id %in% a$gene_id)
table(b$mode) / nrow(b)
# .      dual          switch        rheostat    omnibus_only not_significant
# 0.356898322     0.163300240     0.471071551     0.008729887     0.000000000

# TNFRSF10
a =
  dplyr::left_join(df.open4gene.n_peaks,
                   df.features,
                   by = c("gene_id" = "phenotype_id")) %>%
  dplyr::arrange(dplyr::desc(n_peaks)) %>%
  dplyr::filter(symbol %in% c("TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D"))

b = dplyr::filter(df.open4gene.sig, gene_id %in% a$gene_id)

length(unique(b$peak_id))
# 133
dplyr::group_by(b, peak_id) %>%
  dplyr::summarize(n_genes = length(unique(gene_id))) %>%
  dplyr::ungroup() %>%
  dplyr::summarize(mean(n_genes > 1))
# 0.639

dplyr::group_by(b, gene_id) %>%
  dplyr::summarize()

dplyr::filter(df.open4gene.sig,
              symbol %in% c("TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D")) %>%
  dplyr::distinct(peak_id) %>%
  nrow()
# 133

dplyr::left_join(df.open4gene.n_peaks,
                 df.features,
                 by = c("gene_id" = "phenotype_id")) %>%
  dplyr::arrange(dplyr::desc(n_peaks)) %>%
  dplyr::filter(stringr::str_detect(symbol, "TNFRSF"))

################################################################################

plot_open4gene_sig_count = function(df.open4gene.sig, l1 = TRUE) {
  df.open4gene.sig.count =
    dplyr::filter(df.open4gene.sig, filter_l1_cell_types(cell_type) == l1) %>%
    dplyr::group_by(cell_type, peak_id) %>%
    dplyr::summarize(n_genes = length(gene_id), .groups = "drop") %>%
    dplyr::group_by(cell_type, n_genes) %>%
    dplyr::summarize(n = length(unique(peak_id)), .groups = "drop") %>%
    dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))
  
  open4gene.order =
    dplyr::group_by(df.open4gene.sig.count, cell_type) %>%
    dplyr::summarize(count = sum(n)) %>%
    dplyr::arrange(dplyr::desc(count)) %>%
    dplyr::pull(cell_type)
  
  p.open4gene =
    dplyr::mutate(df.open4gene.sig.count, cell_type = factor(cell_type, levels = rev(open4gene.order))) %>%
    ggplot(aes(n, cell_type)) +
    locusviz::or_missing(
      l1,
      geom_hline(
        yintercept = length(open4gene.order) - 0.5,
        linewidth = 0.25,
        linetype = "dashed",
        color = "grey50"
      )
    ) +
    geom_col(aes(fill = n_genes)) +
    scale_x_continuous(labels = scales::label_comma(), expand = expansion()) +
    scale_y_discrete(labels = if (l1) {
      l1.labels
    } else {
      l2.labels
    }) +
    scale_fill_viridis_c(
      option = "magma",
      breaks = seq(2, 10, by = 2),
      labels = c(seq(2, 8, by = 2), "10+"),
      limits = c(1, 10),
      oob = scales::squish,
      guide = guide_colorbar(title.position = "top")
    ) +
    locusviz::get_default_theme(hide.ytitle = TRUE) +
    theme(
      legend.position.inside = c(1, 0),
      legend.justification = c(1, 0),
      legend.direction = "horizontal",
      legend.key.width = unit(3, "mm"),
      legend.title = element_text(margin = margin(b = 0.5, unit = "mm"))
    ) +
    labs(x = "# CA-Link+ peaks", fill = "# linked genes")
  p.open4gene
}

p.open4gene = plot_open4gene_sig_count(df.open4gene.sig)
p.open4gene.l2 = plot_open4gene_sig_count(df.open4gene.sig, l1 = FALSE)

plot_link_count_by_mode = function(df.open4gene.sig, l1 = TRUE) {
  df.link_count_by_mode =
    dplyr::filter(df.open4gene.sig, filter_l1_cell_types(cell_type) == l1) %>%
    dplyr::group_by(cell_type, mode) %>%
    dplyr::summarize(n = n(), .groups = "drop") %>%
    dplyr::mutate(cell_type = remove_cell_type_prefix(cell_type))
  
  open4gene.order =
    dplyr::group_by(df.link_count_by_mode, cell_type) %>%
    dplyr::summarize(count = sum(n)) %>%
    dplyr::arrange(dplyr::desc(count)) %>%
    dplyr::pull(cell_type)
  
  p.open4gene =
    dplyr::mutate(
      df.link_count_by_mode,
      cell_type = factor(cell_type, levels = rev(open4gene.order)),
      mode = stringr::str_replace(stringr::str_to_title(mode), "_", " "),
      mode = factor(mode, levels = c(
        "Switch", "Rheostat", "Dual", "Omnibus only"
      ))
    ) %>%
    ggplot(aes(n, cell_type)) +
    locusviz::or_missing(
      l1,
      geom_hline(
        yintercept = length(open4gene.order) - 0.5,
        linewidth = 0.25,
        linetype = "dashed",
        color = "grey50"
      )
    ) +
    geom_col(aes(fill = mode)) +
    scale_x_continuous(labels = scales::label_comma(), breaks = scales::breaks_pretty(n = 3), expand = expansion()) +
    scale_y_discrete(labels = if (l1) {
      l1.labels
    } else {
      l2.labels
    }) +
    scale_fill_manual(values = link_mode_colors) +
    locusviz::get_default_theme(
      hide.ytitle = TRUE,
      legend.position = c(1, 0),
      legend.justification = c(1, 0)
    ) +
    labs(x = "# peak-gene links", fill = "Mode")
  p.open4gene
}

p.link.count.mode = plot_link_count_by_mode(df.open4gene.sig)
p.link.count.mode.l2 = plot_link_count_by_mode(df.open4gene.sig, l1 = FALSE)

p.link.count.mode + p.link.count.mode.l2

p.open4gene.n_peaks.mode =
  dplyr::group_by(df.open4gene.sig.peak.mode, mode) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>%
  dplyr::mutate(
    mode = dplyr::case_when(
      mode == "switch_and_rheostat" ~ "Dual",
      mode == "switch_or_rheostat" ~ "Switch/Rheostat",
      mode == "omnibus_only" ~ "Omnibus only",
      TRUE ~ stringr::str_to_title(mode)
    ),
    mode = factor(mode, levels = rev(
      c("Dual", "Rheostat", "Switch", "Switch/Rheostat", "Omnibus only")
    ))
  ) %>%
  ggplot(aes(n, mode, fill = mode)) +
  geom_col() +
  geom_text(aes(label = scales::comma(n)),
            size = 2,
            hjust = -0.25) +
  scale_x_continuous(labels = scales::label_comma(), expand = expansion(c(0, 0.25))) +
  scale_y_discrete(expand = expansion()) +
  scale_fill_manual(values = link_mode_colors) +
  locusviz::get_default_theme(hide.ytitle = TRUE, legend.position = "none") +
  labs(x = "# peaks", fill = "Mode")

################################################################################

df.peak.encode4 = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/atac/peaks/integrated_atac_batch1_5.iterative_merged.encode4.intersect.bed.gz",
  header = FALSE
) %>%
  dplyr::transmute(
    chrom = V1,
    start = V2,
    end = V3,
    encode4 = V12
  ) %>%
  dplyr::group_by(chrom, start, end) %>%
  dplyr::summarize(
    peak_id = stringr::str_c(chrom[1], "-", start[1] + 1, "-", end[1]),
    encode4_detail = stringr::str_c(sort(unique(encode4)), collapse = ","),
    encode4 = ifelse(length(unique(encode4)) == 1, encode4_detail, "Mixed"),
    n_intersect = n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    encode4 = factor(encode4, levels = rev(encode4_levels)),
    link = ifelse(peak_id %in% link_peak_id, "CA-Link+", "CA-Link-")
  ) %>%
  dplyr::left_join(dplyr::select(df.open4gene.sig.peak.mode, peak_id, mode))

sort(table(df.peak.encode4$encode4) / nrow(df.consensus))
#       CA-TF          TF  CA-H3K4me3         PLS     CA-CTCF          CA        pELS       Mixed        dELS
# 0.005324856 0.006230216 0.014398573 0.016628440 0.022131016 0.023254332 0.089731209 0.113371157 0.636950748

p.encode4.link =
  dplyr::group_by(df.peak.encode4, encode4, link) %>%
  dplyr::summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(count, encode4, fill = link)) +
  geom_hline(yintercept = 1.5,
             linetype = "dashed",
             color = "grey50") +
  geom_col(position = position_dodge(width = 1)) +
  scale_x_continuous(labels = scales::label_comma(), expand = expansion()) +
  scale_y_discrete(limits = rev(encode4_levels)) +
  scale_fill_manual(values = link_colors) +
  locusviz::get_default_theme() +
  theme(legend.title = element_blank()) +
  labs(x = "# FinnGen CA peaks", y = "ENCODE4")

compute_link_enrichment = function(data) {
  data <-
    dplyr::count(data, link, encode4) %>%
    dplyr::group_by(link) %>%
    dplyr::mutate(total = sum(n), frac = n / total) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      link = ordered(link, levels = c("CA-Link+", "CA-Link-")),
      encode4 = ordered(encode4, levels = encode4_levels)
    ) %>%
    tidyr::pivot_wider(
      id_cols = "encode4",
      names_from = "link",
      values_from = c("n", "total")
    ) %>%
    dplyr::group_split(encode4) %>%
    purrr::map_dfr(function(data) {
      total_bottom <- "total_CA-Link-"
      n_bottom <-  "n_CA-Link-"
      total_top <- "total_CA-Link+"
      n_top <- "n_CA-Link+"
      null_result <- tibble::tibble(
        encode4 = data$encode4,
        enrichment = NA,
        lower = NA,
        upper = NA
      )
      
      if (!all(c(total_bottom, n_bottom, total_top, n_top) %in% colnames(data))) {
        return(null_result)
      }
      
      m <- with(data, matrix(
        c(
          get(total_bottom) - get(n_bottom),
          get(n_bottom),
          get(total_top) - get(n_top),
          get(n_top)
        ),
        nrow = 2,
        byrow = T
      ))
      
      if (!all(is.finite(m)) || !all(m > 0)) {
        return(null_result)
      }
      # print(epitools::riskratio(m, method = "boot"))
      measure <- epitools::riskratio(m, method = "boot")$measure
      tibble::tibble(
        encode4 = data$encode4,
        enrichment = measure[2, "estimate"],
        lower = measure[2, "lower"],
        upper = measure[2, "upper"],
        n_bottom = data[[n_bottom]],
        total_bottom = data[[total_bottom]],
        n_top = data[[n_top]],
        total_top = data[[total_top]]
      )
    })
  return(data)
}

df.encode4.enrichment =
  dplyr::bind_rows(
    compute_link_enrichment(df.peak.encode4) %>%
      dplyr::mutate(mode = "all"),
    purrr::map_dfr(c(
      "switch", "rheostat", "switch_and_rheostat"
    ), function(mode) {
      dplyr::filter(df.peak.encode4, mode == .env$mode |
                      is.na(mode)) %>%
        compute_link_enrichment() %>%
        dplyr::mutate(
          mode = .env$mode,
          mode = ifelse(mode == "switch_and_rheostat", "dual", mode)
        )
    })
  ) %>%
  dplyr::mutate(mode = factor(mode, levels = rev(link_mode_levels)))

if (!file.exists("tables/ST6_encode4_enrichment.tsv")) {
  dplyr::rename(
    df.encode4.enrichment,
    n_annot_without_link = n_bottom,
    n_without_link = n_bottom,
    n_annot_with_link = n_top,
    n_with_link = total_top
  ) %>%
    export_table("tables/ST6_encode4_enrichment.tsv", "ST6")
}

pd = position_dodge(width = 0.9)
p.encode4.enrichment.all =
  dplyr::filter(df.encode4.enrichment, mode == "all") %>%
  ggplot(aes(enrichment, encode4)) +
  geom_hline(yintercept = 1.5,
             linetype = "dashed",
             color = "grey50") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                 width = 0,
                 color = link_colors["CA-Link+"]) +
  geom_point(color = link_colors["CA-Link+"]) +
  scale_y_discrete(limits = rev(encode4_levels)) +
  locusviz::get_default_theme(hide.ylab = TRUE) +
  labs(x = "CA-Link+ enrichment")

p.encode4.enrichment =
  dplyr::filter(df.encode4.enrichment, mode != "all") %>%
  ggplot(aes(enrichment, encode4, color = mode)) +
  geom_hline(yintercept = 1.5,
             linetype = "dashed",
             color = "grey50") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                 width = 0,
                 position = pd) +
  geom_point(size = 1, position = pd) +
  scale_color_manual(values = link_mode_colors, labels = stringr::str_to_title) +
  scale_y_discrete(limits = rev(encode4_levels)) +
  guides(color = guide_legend(reverse = TRUE)) +
  locusviz::get_default_theme(legend.position = c(1, 0.15),
                              legend.justification = c(1, 0.15)) +
  theme(legend.title = element_blank()) +
  labs(x = "Enrichment", y = "ENCODE4")
p.encode4.enrichment

layout = "
AAB
CDF
CEF
"

plt =
  list(
    ggplot_spacer,
    p.n_nuclei + labs(tag = "d"),
    p.open4gene + labs(tag = "e"),
    p.link.mode.density + labs(tag = "f"),
    p.link.mode.fraction + labs(tag = "g"),
    patchwork::free(p.encode4.enrichment + labs(tag = "h"))
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(heights = c(1, 0.5, 0.5), design = layout)

cowplot::save_plot(
  "figures/Fig1_study_overview_bottom.pdf",
  plt,
  base_height = 121,
  base_width = 180,
  units = "mm"
)

################################################################################

p.ext.atlas =
  list(
    p.gex.n_genes,
    p.gex.n_counts,
    p.n_nuclei.cell_type ,
    p.atac.n_counts,
    p.atac.tss,
    p.peak.width
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig1_qc.pdf",
  p.ext.atlas,
  base_height = 120,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig1_qc.png",
  p.ext.atlas,
  base_height = 120,
  base_width = 180,
  units = "mm",
  dpi = 300
)

p.suppl.open4gene =
  list(p.open4gene.l2, p.open4gene.n_genes, p.open4gene.n_peaks) %>%
  purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig2_open4gene_l2.pdf",
  p.suppl.open4gene,
  base_height = 60,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig2_open4gene_l2.png",
  p.suppl.open4gene,
  base_height = 60,
  base_width = 180,
  units = "mm",
  dpi = 300
)

p.suppl.open4gene.mode =
  list(p.link.count.mode,
       p.link.count.mode.l2,
       p.open4gene.n_peaks.mode,
       p.IRF4, p.TCF4, p.CCDC50
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a")
p.suppl.open4gene.mode

cowplot::save_plot(
  "figures/SFig4_open4gene_mode.pdf",
  p.suppl.open4gene.mode,
  base_height = 120,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig4_open4gene_mode.png",
  p.suppl.open4gene.mode,
  base_height = 120,
  base_width = 180,
  units = "mm",
  dpi = 300
)

p.suppl.open4gene.encode =
  list(p.encode4.link,
       p.encode4.enrichment.all
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a")
p.suppl.open4gene.encode

cowplot::save_plot(
  "figures/SFig5_open4gene_encode_all.pdf",
  p.suppl.open4gene.encode,
  base_height = 60,
  base_width = 120,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig5_open4gene_encode_all.png",
  p.suppl.open4gene.encode,
  base_height = 60,
  base_width = 120,
  units = "mm",
  dpi = 300
)

