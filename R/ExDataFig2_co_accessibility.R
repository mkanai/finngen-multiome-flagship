library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.open4gene.sig = readRDS("data/open4gene.sig.rds")
# TNFRSF10
dplyr::filter(
  df.open4gene.sig,
  symbol %in% c("TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D") &
    cell_type == "predicted.celltype.l1.PBMC"
) %>%
  dplyr::distinct(peak_id) %>%
  rgsutil::write_gsfile("gs://expansion_areas/multiome/misc/co_accessibility/TNFRSF10.sig.peaks.tsv", overwrite = TRUE)

# IFI30
dplyr::filter(
  df.open4gene.sig,
  symbol %in% c("JUND", "UBA52", "IFI30") &
    cell_type == "predicted.celltype.l1.PBMC"
) %>%
  dplyr::distinct(peak_id) %>%
  rgsutil::write_gsfile("gs://expansion_areas/multiome/misc/co_accessibility/IFI30.sig.peaks.tsv", overwrite = TRUE)


################################################################################
TNFRSF10.sig.peaks <- rgsutil::read_gsfile("gs://expansion_areas/multiome/misc/co_accessibility/TNFRSF10.sig.peaks.tsv")[, 1]
IFI30.sig.peaks <- rgsutil::read_gsfile("gs://expansion_areas/multiome/misc/co_accessibility/IFI30.sig.peaks.tsv")[, 1]

df.peaks <- rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/pseudobulk/integrated_atac_batch1_5.fgid.predicted.celltype.l1.PBMC.sum.inv.bed.gz"
)
TNFRSF10.peak_idx = which(df.peaks$gene_id %in% TNFRSF10.sig.peaks)
mat.TNFRSF10.peaks =
  dplyr::mutate(df.peaks, idx = dplyr::row_number()) %>%
  dplyr::filter(min(TNFRSF10.peak_idx) <= idx &
                  idx <= max(TNFRSF10.peak_idx)) %>%
  tidyr::pivot_longer(tidyselect::starts_with("FG"), names_to = "FINNGENID") %>%
  dplyr::select(gene_id:value) %>%
  tidyr::pivot_wider(FINNGENID, names_from = "gene_id") %>%
  tibble::column_to_rownames("FINNGENID") %>%
  as.matrix()

IFI30.peak_idx = which(df.peaks$gene_id %in% IFI30.sig.peaks)
mat.IFI30.peaks =
  dplyr::mutate(df.peaks, idx = dplyr::row_number()) %>%
  dplyr::filter(min(IFI30.peak_idx) <= idx &
                  idx <= max(IFI30.peak_idx)) %>%
  tidyr::pivot_longer(tidyselect::starts_with("FG"), names_to = "FINNGENID") %>%
  dplyr::select(gene_id:value) %>%
  tidyr::pivot_wider(FINNGENID, names_from = "gene_id") %>%
  tibble::column_to_rownames("FINNGENID") %>%
  as.matrix()

compute_coaccessibility <- function(peak_matrix,
                                    cor_threshold = 0.3,
                                    fdr_threshold = 0.05) {
  # Compute correlations
  cor_result <- Hmisc::rcorr(peak_matrix, type = "pearson")
  
  # Get peak names
  peak_names <- rownames(cor_result$r)
  n_peaks <- length(peak_names)
  
  # Create indices for upper triangle
  indices <- which(upper.tri(cor_result$r), arr.ind = TRUE)
  
  # Extract values directly using matrix indexing
  results <- tibble::tibble(
    peak1 = peak_names[indices[, 1]],
    peak2 = peak_names[indices[, 2]],
    correlation = cor_result$r[indices],
    pvalue = cor_result$P[indices],
    n_samples = cor_result$n[indices]
  ) %>%
    dplyr::mutate(
      # FDR correction
      qvalue = p.adjust(pvalue, method = "BH"),
      # Status classification
      status = dplyr::case_when(
        is.na(qvalue) ~ NA_character_,
        qvalue >= fdr_threshold ~ "independent",
        correlation > cor_threshold ~ "co-accessible",
        correlation < -cor_threshold ~ "exclusive",
        TRUE ~ "weak"
      )
    )
  
  return(results)
}

df.coacc.TNFRSF10 = compute_coaccessibility(mat.TNFRSF10.peaks)
rgsutil::write_gsfile(
  df.coacc.TNFRSF10 ,
  "gs://expansion_areas/multiome/misc/co_accessibility/TNFRSF10.coaccessibility.tsv",
  overwrite = TRUE
)

df.coacc.IFI30 = compute_coaccessibility(mat.IFI30.peaks)
rgsutil::write_gsfile(
  df.coacc.IFI30,
  "gs://expansion_areas/multiome/misc/co_accessibility/IFI30.coaccessibility.tsv",
  overwrite = TRUE
)
################################################################################
triangular_coords <- function(x1, x2) {
  t1 <- pmin(x1, x2)
  t2 <- pmax(x1, x2)
  x <- (t1 + t2) / 2
  y <- (t2 - t1) / 2
  return(tibble::tibble(x = x, y = y))
}

read_coacc <- function(path) {
  rgsutil::read_gsfile(path) %>%
    dplyr::filter(status != "independent") %>%
    dplyr::mutate(parse_peak_id(peak1)) %>%
    dplyr::rename(peak_chrom1 = peak_chrom,
                  peak_start1 = peak_start,
                  peak_end1 = peak_end) %>%
    dplyr::mutate(parse_peak_id(peak2)) %>%
    dplyr::rename(peak_chrom2 = peak_chrom,
                  peak_start2 = peak_start,
                  peak_end2 = peak_end) %>%
    dplyr::mutate(
      peak_mid1 = (peak_start1 + peak_end1) / 2,
      peak_mid2 = (peak_start2 + peak_end2) / 2,
      triangular_coords(peak_mid1, peak_mid2)
    )
}



plot_coacc = function(df, bins = 64) {
  ggplot(df, aes(x, y)) +
    # ggrastr::rasterize(stat_summary_hex(aes(z = correlation), fun = median, bins = bins), dpi = 300) +
    ggrastr::rasterize(stat_summary_2d(aes(z = correlation), fun = median, bins = c(256, 256)), dpi = 300) +
    scale_fill_gradientn(colors = BuenColors::jdb_palette("brewer_yes"),
                         limits = c(-1, 1)) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(expand = expansion()) +
    locusviz::get_default_theme(
      hide.xlab = TRUE,
      hide.ylab = TRUE,
      legend.position = c(0, 1),
      legend.justification = c(0, 1)
    ) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(margin = margin(r = 4), vjust = 1),
      legend.direction = "horizontal",
      legend.key.width = unit(5, "mm")
    ) +
    labs(fill = expression(paste("Pearson's ", italic(r))))
}

plot_link_density = function(df, xlim = NULL, n = 512) {
  ggplot(df, aes(peak_mid, symbol)) +
    stat_density(
      aes(fill = after_stat(density)),
      geom = "raster",
      position = "identity",
      n = n
    ) +
    geom_text(
      aes(x = x, y = symbol, label = symbol),
      hjust = 1.1,
      size = 2,
      data = dplyr::reframe(df, x = min(peak_mid), symbol = unique(symbol))
    ) +
    scale_fill_viridis_c() +
    locusviz::get_default_theme(
      hide.xtitle = TRUE,
      hide.ytitle = TRUE,
      legend.position = "none"
    ) +
    theme(axis.line = element_blank(), axis.ticks.y = element_blank()) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_discrete(expand = expansion()) +
    coord_cartesian(xlim = xlim, clip = "off")
}


df.coacc.TNFRSF10 = read_coacc(
  "gs://expansion_areas/multiome/misc/co_accessibility/TNFRSF10.coaccessibility.tsv"
)
df.coacc.IFI30 = read_coacc(
  "gs://expansion_areas/multiome/misc/co_accessibility/IFI30.coaccessibility.tsv"
)

start.TNFRSF10 <- with(df.coacc.TNFRSF10, min(c(peak_mid1, peak_mid2)))
end.TNFRSF10 <- with(df.coacc.TNFRSF10, max(c(peak_mid1, peak_mid2)))
start.IFI30 <- with(df.coacc.IFI30, min(c(peak_mid1, peak_mid2)))
end.IFI30 <- with(df.coacc.IFI30, max(c(peak_mid1, peak_mid2)))


df.link.TNFRSF10 =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% l1.cell_types &
      symbol %in% c("TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D") &
      .env$start.TNFRSF10 < peak_end &
      peak_start < .env$end.TNFRSF10
  ) %>%
  dplyr::mutate(
    peak_mid = (peak_start + peak_end) / 2,
    start = ifelse(tss < peak_mid, tss, peak_mid),
    end = ifelse(tss < peak_mid, peak_mid, tss),
    group = peak_id,
    score = hurdle_count_beta,
    beta =  hurdle_count_beta
  )

df.link.IFI30 =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% l1.cell_types &
      symbol %in% c("JUND", "UBA52", "IFI30") &
      .env$start.IFI30 < peak_end &
      peak_start < .env$end.IFI30
  ) %>%
  dplyr::mutate(
    peak_mid = (peak_start + peak_end) / 2,
    start = ifelse(tss < peak_mid, tss, peak_mid),
    end = ifelse(tss < peak_mid, peak_mid, tss),
    group = peak_id,
    score = hurdle_count_beta,
    beta =  hurdle_count_beta
  )

################################################################################
p.coacc.IFI30 =
  dplyr::filter(df.coacc.IFI30, start.IFI30 <= x &
                  x <= end.IFI30) %>%
  plot_coacc() +
  coord_cartesian(xlim = c(start.IFI30, end.IFI30))

p.gene.IFI30 =
  locusviz::plot_gene_panel("chr19", start.IFI30, end.IFI30, genome_build = "hg38", ) +
  theme(axis.title.x = element_blank())

p.link.IFI30 =
  dplyr::mutate(df.link.IFI30, highlight = FALSE) %>%
  dplyr::filter(cell_type == "predicted.celltype.l1.PBMC") %>%
  plot_links(start.IFI30,
             end.IFI30,
             hide.xtitle = FALSE,
             alpha = c(`FALSE` = 0.5)) +
  labs(x = "Chromosome 19") +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

p.link.density.IFI30 =
  dplyr::distinct(df.link.IFI30, peak_mid, symbol) %>%
  dplyr::mutate(symbol = factor(symbol, levels = rev(c(
    "JUND", "UBA52", "IFI30"
  )))) %>%
  plot_link_density(xlim = c(start.IFI30, end.IFI30)) +
  theme(axis.text.y = element_blank())

p.coacc.TNFRSF10 =
  dplyr::filter(df.coacc.TNFRSF10, start.TNFRSF10 <= x &
                  x <= end.TNFRSF10) %>%
  plot_coacc() +
  coord_cartesian(xlim = c(start.TNFRSF10, end.TNFRSF10)) +
  theme(plot.margin = margin(t = 12))

p.gene.TNFRSF10 =
  locusviz::plot_gene_panel("chr8", start.TNFRSF10, end.TNFRSF10, genome_build = "hg38", ) +
  theme(axis.title.x = element_blank())

p.link.TNFRSF10 =
  dplyr::mutate(df.link.TNFRSF10, highlight = FALSE) %>%
  dplyr::filter(cell_type == "predicted.celltype.l1.PBMC") %>%
  plot_links(start.TNFRSF10,
             end.TNFRSF10,
             hide.xtitle = FALSE,
             alpha = c(`FALSE` = 1)) +
  labs(x = "Chromosome 8") +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

p.link.density.TNFRSF10 =
  dplyr::distinct(df.link.TNFRSF10, peak_mid, symbol) %>%
  dplyr::mutate(symbol = factor(symbol, levels = rev(
    c("TNFRSF10B", "TNFRSF10D", "TNFRSF10A", "TNFRSF10C")
  ))) %>%
  plot_link_density(xlim = c(start.TNFRSF10, end.TNFRSF10)) +
  theme(axis.text.y = element_blank())

p.coaccessibility =
  list(
    p.coacc.IFI30 + labs(tag = "a"),
    p.link.density.IFI30,
    p.gene.IFI30,
    p.link.IFI30,
    p.coacc.TNFRSF10 + labs(tag = "b"),
    p.link.density.TNFRSF10,
    p.gene.TNFRSF10,
    p.link.TNFRSF10
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1, heights = rep(c(0.5, 0.3, 0.2, 0.3), 2))
p.coaccessibility

cowplot::save_plot(
  "figures/ExtendedDataFig2_coacc.pdf",
  p.coaccessibility,
  base_height = 170,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/ExtendedDataFig2_coacc.png",
  p.coaccessibility,
  base_height = 170,
  base_width = 180,
  units = "mm",
  dpi = 300
)
