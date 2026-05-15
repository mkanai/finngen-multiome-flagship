library(ggplot2)

ggplot_spacer = ggplot() + locusviz::get_default_theme() + theme(plot.tag = element_blank())

cell_colors <- c(
  # B cell lineage (Purples - based on BuenColors B="#BA7FD0")
  "predicted.celltype.l1.B" = "#BA7FD0",
  # Same as L1
  "predicted.celltype.l2.B_naive" = "#BA7FD0",
  # Darker purple
  "predicted.celltype.l2.B_intermediate" = "#A65AC2",
  # Even darker purple
  "predicted.celltype.l2.B_memory" = "#8B3F99",
  # Lighter purple
  "predicted.celltype.l2.Plasmablast" = "#D4A5E3",
  
  # CD4 T cell lineage (Greens - to fill spectrum gap)
  # Green
  "predicted.celltype.l1.CD4_T" = "#2CA02C",
  # Same as L1
  "predicted.celltype.l2.CD4_Naive" = "#2CA02C",
  # Forest green
  "predicted.celltype.l2.CD4_TCM" = "#228B22",
  # Dark green
  "predicted.celltype.l2.CD4_TEM" = "#006400",
  # Light green
  "predicted.celltype.l2.CD4_CTL" = "#5CB85C",
  # Very light green
  "predicted.celltype.l2.CD4_Proliferating" = "#90EE90",
  # Very dark green
  "predicted.celltype.l2.Treg" = "#013220",
  
  # CD8 T cell lineage (Blues - based on BuenColors CD8="#001588")
  # Blue (more visible than dark navy)
  "predicted.celltype.l1.CD8_T" = "#1F77B4",
  # Same as L1
  "predicted.celltype.l2.CD8_Naive" = "#1F77B4",
  # Darker blue
  "predicted.celltype.l2.CD8_TCM" = "#17508F",
  # Even darker blue
  "predicted.celltype.l2.CD8_TEM" = "#0F3A6B",
  # Lighter blue
  "predicted.celltype.l2.CD8_Proliferating" = "#5FA2CE",
  
  # NK cells (Dark Purples/Magentas - based on BuenColors NK="#490C65")
  # Magenta (more visible)
  "predicted.celltype.l1.NK" = "#E377C2",
  # Same as L1
  "predicted.celltype.l2.NK" = "#E377C2",
  # Lighter pink
  "predicted.celltype.l2.NK_CD56bright" = "#F7B6D2",
  # Darker magenta
  "predicted.celltype.l2.NK_Proliferating" = "#C5589F",
  
  # Monocytes (Oranges - based on BuenColors mono="#FF5A00")
  # Orange
  "predicted.celltype.l1.Mono" = "#FF7F0E",
  # Same as L1
  "predicted.celltype.l2.CD14_Mono" = "#FF7F0E",
  # Lighter orange
  "predicted.celltype.l2.CD16_Mono" = "#FFBB78",
  
  # Dendritic cells (Reds - to fill spectrum gap)
  # Red
  "predicted.celltype.l1.DC" = "#D62728",
  # Darker red
  "predicted.celltype.l2.cDC1" = "#C71E1F",
  # Same as L1
  "predicted.celltype.l2.cDC2" = "#D62728",
  # Light red/pink
  "predicted.celltype.l2.pDC" = "#FF9896",
  # Salmon
  "predicted.celltype.l2.ASDC" = "#E58080",
  
  # Other T cells (Cyans - distinct from major T cells)
  # Cyan
  "predicted.celltype.l1.other_T" = "#17BECF",
  # Same as L1
  "predicted.celltype.l2.MAIT" = "#17BECF",
  # Darker cyan
  "predicted.celltype.l2.gdT" = "#0F94A8",
  # Lighter cyan
  "predicted.celltype.l2.dnT" = "#5FD3DD",
  
  # Progenitors and other cells
  # Yellow-green
  "predicted.celltype.l2.HSPC" = "#BCBD22",
  # Dark red
  "predicted.celltype.l2.Eryth" = "#8B0000",
  # Light salmon
  "predicted.celltype.l2.Platelet" = "#FFA07A",
  # Very light cyan
  "predicted.celltype.l2.ILC" = "#9EDAE5",
  
  # General categories
  # Gray
  "predicted.celltype.l1.PBMC" = "#7F7F7F",
  # Light gray
  "predicted.celltype.l1.other" = "#C7C7C7",
  # Dark gray
  "predicted.celltype.l2.Doublet" = "#4D4D4D"
)

l1.colors = rep(cell_colors[stringr::str_detect(names(cell_colors), "^predicted\\.celltype\\.l1\\.")], 2)
names(l1.colors)[seq_len(length(l1.colors) / 2)] = stringr::str_remove(names(l1.colors)[seq_len(length(l1.colors) / 2)], "^predicted\\.celltype\\.l1\\.")
l2.colors = rep(cell_colors[stringr::str_detect(names(cell_colors), "^predicted\\.celltype\\.l2\\.")], 2)
names(l2.colors)[seq_len(length(l2.colors) / 2)] = stringr::str_remove(names(l2.colors)[seq_len(length(l2.colors) / 2)], "^predicted\\.celltype\\.l2\\.")

l1.order = c("PBMC",
             "B",
             "CD4_T",
             "CD8_T",
             "Mono",
             "DC",
             "NK",
             "other",
             "other_T")

l1.cell_types = stringr::str_c("predicted.celltype.l1.", l1.order)
l1.cell_types.primary = setdiff(
  l1.cell_types,
  c(
    "predicted.celltype.l1.PBMC",
    "predicted.celltype.l1.other",
    "predicted.celltype.l1.other_T"
  )
)

l1.labels = c(
  "PBMC" = "PBMC",
  "CD4_T" = "CD4+ T",
  "CD8_T" = "CD8+ T",
  "Mono" = "Mono",
  "NK" = "NK",
  "B" = "B",
  "DC" = "DC",
  "other_T" = "Other T",
  "other" = "Other"
)

l2.cell_types = c(
  "predicted.celltype.l2.B_intermediate",
  "predicted.celltype.l2.B_memory",
  "predicted.celltype.l2.B_naive",
  "predicted.celltype.l2.CD14_Mono",
  "predicted.celltype.l2.CD16_Mono",
  "predicted.celltype.l2.CD4_CTL",
  "predicted.celltype.l2.CD4_Naive",
  "predicted.celltype.l2.CD4_TCM",
  "predicted.celltype.l2.CD4_TEM",
  "predicted.celltype.l2.CD8_Naive",
  "predicted.celltype.l2.CD8_TEM",
  "predicted.celltype.l2.HSPC",
  "predicted.celltype.l2.ILC",
  "predicted.celltype.l2.MAIT",
  "predicted.celltype.l2.NK",
  "predicted.celltype.l2.NK_CD56bright",
  "predicted.celltype.l2.NK_Proliferating",
  "predicted.celltype.l2.Plasmablast",
  "predicted.celltype.l2.Platelet",
  "predicted.celltype.l2.Treg",
  "predicted.celltype.l2.cDC2",
  "predicted.celltype.l2.dnT",
  "predicted.celltype.l2.gdT",
  "predicted.celltype.l2.pDC"
)

l2.labels = c(
  "ASDC" = "ASDC",
  "B_intermediate" = "B Intermediate",
  "B_memory" = "B Memory",
  "B_naive" = "B Naive",
  "CD14_Mono" = "CD14+ Mono",
  "CD16_Mono" = "CD16+ Mono",
  "CD4_CTL" = "CD4+ CTL",
  "CD4_Naive" = "CD4+ Naive",
  "CD4_Proliferating" = "CD4+ Proliferating",
  "CD4_TCM" = "CD4+ TCM",
  "CD4_TEM" = "CD4+ TEM",
  "CD8_Naive" = "CD8+ Naive",
  "CD8_Proliferating" = "CD8+ Proliferating",
  "CD8_TCM" = "CD8+ TCM",
  "CD8_TEM" = "CD8+ TEM",
  "HSPC" = "HSPC",
  "ILC" = "ILC",
  "MAIT" = "MAIT",
  "NK" = "NK",
  "NK_CD56bright" = "NK CD56bright",
  "NK_Proliferating" = "NK Proliferating",
  "Plasmablast" = "Plasmablast",
  "Platelet" = "Platelet",
  "Treg" = "Treg",
  "cDC1" = "cDC1",
  "cDC2" = "cDC2",
  "dnT" = "dnT",
  "gdT" = "gdT",
  "pDC" = "pDC"
)


annot_levels = c("UTR5",
                 "UTR3",
                 "Intron",
                 "pLoF",
                 "Missense",
                 "Synonymous",
                 "cCRE",
                 "Non-cCRE")
encode4_levels = c("PLS",
                   "pELS",
                   "dELS",
                   "CA-CTCF",
                   "CA-H3K4me3",
                   "CA-TF",
                   "CA",
                   "TF",
                   "Mixed")

variant_heterogeneity_colors = c(
  "shared_consistent" = as.character(shades::saturation(BuenColors::jdb_palette("corona")[5], 0.2)),
  "shared_heterogeneous" = BuenColors::jdb_palette("corona")[5],
  "shared_opposite" = as.character(shades::saturation(BuenColors::jdb_palette("corona")[5], 0.8)),
  "distinct_variants" = BuenColors::jdb_palette("corona")[7]
)
variant_heterogeneity_levels = names(variant_heterogeneity_colors)

cell_type_specificity_colors = c(
  "Cross-lineage shared" = BuenColors::jdb_palette("corona")[1],
  "Likely shared but underpowered" = as.character(shades::saturation(BuenColors::jdb_palette("corona")[1], 0.2)),
  "Lineage-specific" = BuenColors::jdb_palette("corona")[2],
  "T-cell-specific" = BuenColors::jdb_palette("corona")[3],
  "Single cell-type" = BuenColors::jdb_palette("corona")[4]
)
cell_type_specificity_levels = names(cell_type_specificity_colors)

qtl_mechanism_colors = c(
  "Local Cascade" = BuenColors::jdb_palette("corona")[10],
  "Positional Cascade" = BuenColors::jdb_palette("corona")[9],
  "Distal Cascade" = BuenColors::jdb_palette("corona")[13],
  "caQTL + eQTL (No Link)" = BuenColors::jdb_palette("corona")[12],
  "Only caQTL (With Link)" = BuenColors::jdb_palette("brewer_orange")[8],
  "Only caQTL (No Link)" = BuenColors::jdb_palette("brewer_orange")[6],
  "Only eQTL" = BuenColors::jdb_palette("brewer_green")[3],
  "Any pQTL" = "grey10",
  "Only pQTL" = "grey30",
  "No molQTL" = "grey50"
)
qtl_mechanism_levels = names(qtl_mechanism_colors)

link_mode_colors = c(
  "all" = "grey10",
  "rheostat" = "#26547c",
  "dual" = "#edae49",
  "switch" = "#ef476f",
  "Rheostat" = "#26547c",
  "Dual" = "#edae49",
  "Switch" = "#ef476f",
  "Switch/Rheostat" = "#267726",
  "omnibus_only" = "grey50",
  "Omnibus_only" = "grey50",
  "Omnibus only" = "grey50"
)
link_mode_levels = names(link_mode_colors)

model.cols = c(
  "zero" = BuenColors::jdb_palette("corona")[7],
  "count" = BuenColors::jdb_palette("corona")[9],
  "Zero" = BuenColors::jdb_palette("corona")[7],
  "Count" = BuenColors::jdb_palette("corona")[9]
)

peak_overlap_colors = c(
  "Peak overlap" = BuenColors::jdb_palette("brewer_red")[7],
  "No overlap" = "grey50"
)
caqtl_colors = c(
  "For overlapping peak" = BuenColors::jdb_palette("brewer_orange")[7],
  "For non-overlapping peak" = BuenColors::jdb_palette("brewer_orange")[5],
  "No detected caQTL" = "grey50"
)
link_colors = c(
  "For caQTL peak" = BuenColors::jdb_palette("brewer_blue")[7],
  "For overlapping peak (link)" = BuenColors::jdb_palette("brewer_blue")[5],
  "No link" = "grey50",
  "CA-Link+" = BuenColors::jdb_palette("brewer_blue")[7],
  "CA-Link-" = "grey50"
)
eqtl_colors = c(
  "For linked gene"  = BuenColors::jdb_palette("brewer_green")[7],
  "For non-linked gene"  = BuenColors::jdb_palette("brewer_green")[5],
  "eQTL"  = BuenColors::jdb_palette("brewer_green")[3],
  "No detected eQTL"  = "grey50"
)

gene_sankey_colors = c(cell_type_specificity_colors, variant_heterogeneity_colors)
variant_sankey_colors = c(
  peak_overlap_colors,
  caqtl_colors,
  link_colors,
  eqtl_colors,
  qtl_mechanism_colors,
  cell_type_specificity_colors
)

qtl.shapes = c("caQTL" = 16, "eQTL" = 17)
qtl.colors = c(
  "caQTL" = BuenColors::jdb_palette("corona")[1],
  "eQTL" = BuenColors::jdb_palette("corona")[2],
  "pQTL" = BuenColors::jdb_palette("corona")[3],
  "caQTL & eQTL" = BuenColors::jdb_palette("corona")[5],
  "eQTL & pQTL" = BuenColors::jdb_palette("corona")[4],
  "caQTL & pQTL" = BuenColors::jdb_palette("corona")[10],
  "caQTL & eQTL & pQTL" = BuenColors::jdb_palette("corona")[8]
)

immune_traits = c(
  # 255 - Hypothyroidism, strict autoimmune
  "E4_HYTHY_AI_STRICT",
  # 206 - Autoimmune diseases
  "AUTOIMMUNE",
  # 106 - Asthma (more control exclusions)
  "J10_ASTHMA_EXMORE",
  # 97 - Atopic dermatitis
  "L12_ATOPIC",
  # 96 - Inflammatory bowel disease, strict
  "K11_IBD_STRICT",
  # 75 - Dermatitis and eczema
  "L12_DERMATITISECZEMA",
  # 55 - Psoriasis
  "L12_PSORIASIS",
  # 54 - Ulcerative colitis (strict)
  "K11_UC_STRICT2",
  # 45 - Rheumatoid arthritis
  "M13_RHEUMA",
  # 17 - Multiple Sclerosis
  "G6_MS",
  # 16 - Crohn disease (strict)
  "K11_CD_STRICT2"
)

cardiometabolic_traits = c(
  # 424 - Hypertension
  "I9_HYPTENS",
  # 392 - Type 2 diabetes
  "T2D",
  # 310 - Statin medication
  "RX_STATIN",
  # 262 - Varicose veins
  "I9_VARICVE",
  # 177 - Atrial fibrillation and flutter
  "I9_AF",
  # 145 - Ischaemic heart disease, wide
  "I9_IHD",
  # 84 - Venous thromboembolism
  "I9_VTE",
  # 84 - Myocardial infarction, strict
  "I9_MI_STRICT"
)

cancer_traits = c(
  # 135 - Prostate cancer
  "C3_PROSTATE_EXALLC",
  # 119 - Malignant neoplasm of skin
  "C3_SKIN_EXALLC",
  # 119 - Basal cell carcinoma
  "C3_BASAL_CELL_CARCINOMA_EXALLC",
  # 105 - Other skin cancers
  "C3_OTHER_SKIN_EXALLC",
  # 98 - Breast cancer
  "C3_BREAST_EXALLC",
  # 85 - All malignant neoplasms
  "C3_CANCER_EXALLC"
)

blood_counts = c(
  # 1297
  "Platelets" = "3007461",
  # 939
  "Leukocytes" = "3010813",
  # 765
  "Erythrocytes" = "3026361",
  # 411
  "Hemoglobin" = "3000963",
  # 403
  "Hematocrit" = "3009542",
  # 305
  "Monocytes" = "3001604",
  # 241
  "Basophils" = "3006315",
  # 222
  "Eosinophils" = "3013115",
  # 211
  "Lymphocytes" = "3019198",
  # 177
  "Neutrophils" = "3013650"
)

lipids_glucose = c(
  # 921 - HDL cholesterol
  "HDLC" = "3023602",
  # 585 - Fasting triglycerides
  "FastTG" = "3048773",
  # 519 - Hemoglobin A1c
  "HbA1c" = "3004410",
  # 464 - Fasting glucose
  "FastGlucose" = "3018251",
  # 304 - Total cholesterol
  "TC" = "3019900",
  # 219 - Non-HDL cholesterol (derived)
  "nonHDL" = "nonHDL",
  # 189 - LDL cholesterol
  "LDLC" = "3001308",
  # 91 - Glucose
  "Glucose" = "3013826",
  # 77 - Triglycerides
  "TG" = "3025839"
)

parse_peak_id = function(peak_id) {
  x = stringr::str_split_fixed(peak_id, "-", 3)
  return(tibble::tibble(
    peak_chrom = x[, 1],
    peak_start = as.numeric(x[, 2]),
    peak_end = as.numeric(x[, 3])
  ))
}

parse_cell_type = function(path) {
  stringr::str_replace(path,
                       "^.*\\.(predicted\\.celltype\\.l[12]\\.[^\\.]*)\\..*$",
                       "\\1")
}

remove_cell_type_prefix = function(cell_type, keep_level = FALSE) {
  pattern = ifelse(keep_level,
                   "predicted\\.celltype\\.",
                   "predicted\\.celltype\\.l[12]\\.")
  return(stringr::str_replace_all(cell_type, pattern, ""))
}

filter_l1_cell_types = function(cell_type) {
  stringr::str_detect(cell_type, "predicted\\.celltype\\.l1\\.")
}

munge_zfile = function(zfile_pattern, cell_types, gene_id) {
  if (stringr::str_starts(gene_id, "ENSG")) {
    chrom = dplyr::filter(df.features, phenotype_id == .env$gene_id) %>%
      dplyr::pull(chrom)
  } else {
    chrom = stringr::str_split_fixed(gene_id, "-", 2)[1, 1]
  }
  
  purrr::map_dfr(cell_types, function(cell_type) {
    rgsutil::read_gsfile(sprintf(zfile_pattern, cell_type, cell_type, chrom[1], gene_id)) %>%
      dplyr::mutate(cell_type = .env$cell_type,
                    gene_id = .env$gene_id)
  })
}

import_peak_track = function(path_pattern, cell_types) {
  purrr::pmap_dfr(tidyr::crossing(i = seq(0, 2), cell_type = cell_types), function(i, cell_type) {
    # gs://expansion_areas/multiome/misc/qtl_boxplot/*.bw
    rtracklayer::import(sprintf(path_pattern, cell_type, i)) %>%
      GenomicRanges::as.data.frame() %>%
      dplyr::mutate(GT = i, cell_type = cell_type)
  }) %>%
    dplyr::mutate(GT = factor(GT), cell_type = remove_cell_type_prefix(cell_type)) %>%
    tidyr::pivot_longer(cols = c(start, end), values_to = "x")
}

munge_bezier_df = function(df) {
  # ref: https://github.com/stuart-lab/signac/blob/HEAD/R/visualization.R
  df$group =  seq_len(length.out = nrow(x = df))
  data.frame(
    x = c(df$start, (df$start + df$end) / 2, df$end),
    y = c(rep(x = 0, nrow(x = df)), -abs(df$beta), rep(x = 0, nrow(x = df))),
    # y = c(rep(x = 0, nrow(x = df)), rep(x = -1, nrow(x = df)), rep(x = 0, nrow(x = df))),
    group = rep(x = df$group, 3)
  ) %>%
    dplyr::left_join(dplyr::select(df, -start, -end), by = "group")
}

guide_rug <- function(x, ...) {
  legendry::primitive_ticks(key = legendry::key_manual(x), ...)
}

peak_background = function(xmin, xmax, flip_y = FALSE) {
  list(geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ifelse(flip_y, -Inf, 0),
      ymax = ifelse(flip_y, 0, Inf)
    ),
    fill = "grey90",
    data = data.frame(xmin = xmin, xmax = xmax)
  ))
}

plot_links = function(df,
                      start,
                      end,
                      highlight_pos = NULL,
                      peak.ranges = NULL,
                      xbreaks = ggplot2::waiver(),
                      hide.xtitle = TRUE,
                      background.layers = NULL,
                      alpha = c(`TRUE` = 1, `FALSE` = 0.2),
                      linewidth = c(`TRUE` = 1.5, `FALSE` = 0.5)) {
  if (!("highlight" %in% colnames(df))) {
    df = dplyr::mutate(df, highlight = FALSE)
  }
  max_abs_beta = max(abs(df$beta), na.rm = TRUE)

  munge_bezier_df(df) %>%
    ggplot() +
    background.layers +
    ggforce::geom_bezier(aes(
      x,
      y,
      group = group,
      color = beta,
      alpha = highlight,
      linewidth = highlight
    )) +
    locusviz::or_missing(
      !is.null(highlight_pos),
      geom_vline(
        xintercept = highlight_pos,
        linetype = "dashed",
        color = "grey50"
      )
    ) +
    scale_color_gradientn(
      colors = BuenColors::jdb_palette("solar_extra"),
      limits = c(-max_abs_beta, max_abs_beta)
    ) +
    scale_x_continuous(labels = scales::label_comma(), breaks = xbreaks) +
    scale_y_continuous(expand = expansion(), labels = ~ {
      -.x * 2
    }) +
    scale_alpha_manual(values = alpha, guide = "none") +
    scale_linewidth_manual(values = linewidth, guide = "none") +
    locusviz::get_default_theme(hide.xtitle = hide.xtitle) +
    theme(
      legend.position.inside = c(0, 0),
      legend.justification.inside = c(0, 0),
      legend.title = element_blank(),
      legend.direction = "horizontal",
      legend.key.width = unit(5, "mm")
    ) +
    coord_cartesian(xlim = c(start, end), clip = "off") +
    locusviz::or_missing(!is.null(peak.ranges), guides(x = legendry::compose_ontop(
      "axis",
      legendry::primitive_box(peak.ranges, theme = legendry::theme_guide(box = element_rect(
        fill = NA, color = NA
      )))
    ))) +
    guides(x.sec = guide_rug(unique(df$peak_mid), theme = theme(axis.ticks.length.x = unit(1, "mm")))) +
    labs(y = expression("|" ~ italic(beta)[Link] ~ "|"))
}

plot_peak = function(df.track,
                     cell_type,
                     start,
                     end,
                     highlight_pos,
                     fill_colors = NULL,
                     hide.xtitle = TRUE,
                     hide.ytitle = TRUE,
                     reverse.GT = FALSE,
                     rasterize = FALSE,
                     rasterize.dpi = 300,
                     x.breaks = waiver()) {
  if (is.null(fill_colors)) {
    fill_colors = locusviz::distinct_shades(l1.colors[cell_type])
  } else if (length(fill_colors) == 1) {
    fill_colors = locusviz::distinct_shades(fill_colors)
  } else if (length(fill_colors) != 3) {
    stop("fill_colors must specify one or three colors")
  }
  names(fill_colors) = seq(0, 2)
  
  df.track = dplyr::filter(df.track, cell_type == .env$cell_type) %>%
    dplyr::filter(.env$start <= x &
                    x <= .env$end)
  if (reverse.GT) {
    df.track = dplyr::mutate(df.track, GT = factor(GT, levels = seq(2, 0)))
  }
  
  rasterize_f <- ifelse(rasterize, function(p) {
    ggrastr::rasterize(p, dpi = rasterize.dpi)
  }, function(p) {
    p
  })
  
  ggplot() +
    rasterize_f(geom_area(aes(x, score, fill = GT), position = "identity", data = df.track)) +
    locusviz::or_missing(
      !is.null(highlight_pos),
      geom_vline(
        xintercept = highlight_pos,
        linetype = "dashed",
        color = "grey50"
      )
    ) +
    scale_x_continuous(
      breaks = x.breaks,
      expand = expansion(),
      labels = scales::label_comma()
    ) +
    scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    scale_fill_manual(breaks = seq(0, 2), values = fill_colors) +
    coord_cartesian(xlim = c(start, end)) +
    locusviz::get_default_theme(
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      hide.xtitle = hide.xtitle,
      hide.ytitle = hide.ytitle
    ) +
    theme(plot.title = element_text(hjust = 0, margin = margin()),
          legend.title = element_blank()) +
    labs(y = "Norm. CA")
}

plot_forest = function(df,
                       colors = l1.colors,
                       labels = l1.labels,
                       hide.xtitle = FALSE,
                       hide.ytitle = TRUE,
                       legend.position = "none",
                       legend.justification = NULL) {
  data =
    dplyr::mutate(
      df,
      cell_type = remove_cell_type_prefix(cell_type),
      cell_type = factor(cell_type, levels = cell_type[order(beta)])
    )
  if (!("sig" %in% colnames(data))) {
    data = dplyr::mutate(data, sig = factor(
      p < 5e-8,
      levels = c("FALSE", "TRUE"),
      labels = c("P >= 5e-8", "P < 5e-8")
    ))
  }
  ggplot(data, aes(beta, cell_type, color = cell_type)) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               color = "grey50") +
    geom_errorbar(aes(xmin = beta - se, xmax = beta + se), width = 0) +
    geom_point(aes(shape = sig)) +
    scale_color_manual(values = colors, guide = "none") +
    scale_shape_manual(values = c(
      "P >= 5e-8" = 4,
      "P < 5e-8" = 16,
      "FDR >= 0.05" = 4,
      "FDR < 0.05" = 16
    )) +
    scale_y_discrete(labels = labels) +
    locusviz::get_default_theme(
      legend.position = legend.position,
      legend.justification = legend.justification,
      hide.xtitle = hide.xtitle,
      hide.ytitle = hide.ytitle
    ) +
    theme(legend.title = element_blank()) +
    labs(x = "Effect size", shape = "Significance")
}

.google_ss_url <- NULL
if (file.exists(here::here("local/R/secrets.R"))) {
  source(here::here("local/R/secrets.R"))
}

export_table = function(df,
                        filename,
                        sheet = NULL,
                        range = "A3",
                        save_googlesheet = TRUE,
                        google_ss_url = .google_ss_url) {
  if (is.null(sheet) & save_googlesheet) {
    stop("sheet is null")
  }
  if (save_googlesheet && is.null(google_ss_url)) {
    stop("google_ss_url is NULL; set it in local/R/secrets.R or pass google_ss_url=")
  }

  data.table::fwrite(
    df,
    filename,
    quote = FALSE,
    row.names = FALSE,
    na = "NA",
    sep = "\t"
  )
  if (save_googlesheet) {
    googlesheets4::range_write(
      ss = google_ss_url,
      data = df,
      sheet = sheet,
      range = range
    )
  }
}

cache_panel <- function(name, expr, deps, dir = ".cache/") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  xfun::cache_rds(expr,
                  file = paste0(name, ".rds"),
                  dir = dir,
                  hash = deps)
}


df.tss = read.table(
  here::here("data/refdata-gex-GRCh38-2020-A.genes.tss.bed"),
  F,
  sep = "\t"
) %>%
  dplyr::transmute(phenotype_id = V4, tss = V3)

df.tss_full = read.table(
  here::here("data/refdata-gex-GRCh38-2020-A.genes.tss_full.bed"),
  F,
  sep = "\t"
) %>%
  dplyr::transmute(
    phenotype_id = V4,
    chrom = V1,
    start = V2 + 1,
    end = V3 + 1
  )

df.features = rgsutil:::fread_wrapper(here::here("data/features.gex.tsv.gz")) %>%
  dplyr::transmute(
    phenotype_id = V1,
    symbol = V2,
    chrom = V4,
    gene_type = V7
  ) %>%
  # Use collapsed model for start/end
  dplyr::left_join(df.tss_full) %>%
  dplyr::left_join(df.tss) %>%
  dplyr::mutate(strand = ifelse(abs(tss - start) < abs(tss - end), "+", "-"))

gene_symbols = df.features$symbol %>%
  purrr::set_names(df.features$phenotype_id)
