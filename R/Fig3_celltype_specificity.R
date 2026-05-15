library(dplyr)
library(ggplot2)
library(ggsankeyfier)
library(patchwork)
source(here::here("R/const.R"))

df.gene = data.table::fread("data/cascade_results/gene_categorization.tsv.gz",
                            data.table = FALSE) %>%
  dplyr::select(-variant_heterogeneity) %>%
  dplyr::rename(variant_heterogeneity = variant_heterogeneity_pattern) %>%
  dplyr::filter(cell_type_specificity != "No significance")

df.peak = data.table::fread("data/cascade_results/peak_categorization.tsv.gz",
                            data.table = FALSE) %>%
  dplyr::select(-variant_heterogeneity) %>%
  dplyr::rename(variant_heterogeneity = variant_heterogeneity_pattern) %>%
  dplyr::filter(cell_type_specificity != "No significance")

if (!file.exists("tables/ST18_egene_specificity.tsv")) {
  dplyr::left_join(df.gene, df.features, by = c("gene_id" = "phenotype_id")) %>%
    dplyr::select(gene_id, cell_type_specificity, n_significant_cts,n_tested_cts,variant_heterogeneity, top_variant) %>%
    export_table("tables/ST18_egene_specificity.tsv", save_googlesheet = FALSE)
}
if (!file.exists("tables/ST19_capeak_specificity.tsv")) {
  dplyr::select(df.peak, peak_id, cell_type_specificity, n_significant_cts,n_tested_cts,variant_heterogeneity, top_variant) %>%
    export_table("tables/ST19_capeak_specificity.tsv", save_googlesheet = FALSE)
}

plot_gene_sankey = function(df,
                            label.cex = 2,
                            width = 0.15) {
  pos_sankey <- ggsankeyfier::position_sankey(
    width = width,
    split_nodes = FALSE,
    align = "center",
    order = "as_is",
    v_space = "auto",
    h_space = "auto",
    direction = "backward"
  )
  pos_sankey_left <- ggsankeyfier::position_sankey(
    width = width,
    v_space = "auto",
    align = "center",
    order = "as_is",
    nudge_x = -0.05
  )
  pos_sankey_right <- ggsankeyfier::position_sankey(
    width = width,
    v_space = "auto",
    align = "center",
    order = "as_is",
    nudge_x = 0.05
  )
  
  data =
    dplyr::group_by(df, cell_type_specificity, variant_heterogeneity) %>%
    dplyr::summarize(count = n(), .groups = "drop") %>%
    dplyr::mutate(
      cell_type_specificity = factor(cell_type_specificity, levels = cell_type_specificity_levels),
      variant_heterogeneity = factor(variant_heterogeneity, levels = variant_heterogeneity_levels)
    ) %>%
    ggsankeyfier::pivot_stages_longer(
      stages_from = c("cell_type_specificity", "variant_heterogeneity"),
      values_from = "count"
    ) %>%
    dplyr::group_by(stage) %>%
    dplyr::mutate(frac_edge = count / sum(count)) %>%
    dplyr::group_by(edge_id) %>%
    dplyr::mutate(label_frac_edge = ifelse(frac_edge > 0.05 &
                                             all(!is.na(node)), frac_edge, NA)) %>%
    dplyr::group_by(node, stage) %>%
    dplyr::mutate(frac_node = sum(frac_edge)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!(connector == "to" & is.na(node)))
  
  p =
    ggplot(data,
           aes(
             x = stage,
             y = count,
             group = node,
             connector = connector,
             edge_id = edge_id
           )) +
    ggsankeyfier::geom_sankeyedge(aes(fill = node), position = pos_sankey) +
    ggsankeyfier::geom_sankeynode(aes(fill = node), position = pos_sankey) +
    geom_text(
      aes(
        label = scales::percent(frac_node, accuracy = 1),
        color = node
      ),
      stat = "sankeynode",
      position = pos_sankey,
      cex = label.cex,
      angle = 90
    ) +
    geom_text(
      aes(label = scales::percent(label_frac_edge, accuracy = 1)),
      stat = "sankeyedge",
      position = pos_sankey_right,
      cex = label.cex,
      angle = 90,
      color = "grey20"
    ) +
    # geom_text(
    #   aes(label = node),
    #   stat = "sankeynode",
    #   position = pos_sankey_text,
    #   cex = 2,
    #   angle = 90
    # ) +
    # geom_text(
    #   aes(label = node),
    #   stat = "sankeynode",
    #   position = pos_sankey_text2,
    #   cex = 2,
    #   angle = 90
    # ) +
    scale_fill_manual(values = gene_sankey_colors) +
    scale_color_manual(values =  purrr::set_names(
      ifelse(shades::lightness(gene_sankey_colors) > 30, "black", "white"),
      names(gene_sankey_colors)
    )) +
    scale_x_discrete(expand = expansion(add = 1)) +
    theme_void() +
    theme(legend.position = "none")
  return(p)
}

dplyr::group_by(df.gene, cell_type_specificity) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(frac = count / sum(count))
# # A tibble: 5 × 3
# cell_type_specificity          count   frac
# <chr>                          <int>  <dbl>
#   1 Cross-lineage shared           10594 0.509
# 2 Likely shared but underpowered  4681 0.225
# 3 Lineage-specific                1827 0.0877
# 4 Single cell-type                2859 0.137
# 5 T-cell-specific                  868 0.0417

dplyr::group_by(df.peak, cell_type_specificity) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(frac = count / sum(count))
# cell_type_specificity          count   frac
# <chr>                          <int>  <dbl>
#   1 Cross-lineage shared           68438 0.325
# 2 Likely shared but underpowered 93172 0.442
# 3 Lineage-specific               15445 0.0733
# 4 Single cell-type               28797 0.137
# 5 T-cell-specific                 4732 0.0225

dplyr::group_by(df.gene, variant_heterogeneity) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(frac = count / sum(count))
# variant_heterogeneity count   frac
# <chr>                 <int>  <dbl>
#   1 distinct_variants      2185 0.105
# 2 shared_consistent     10402 0.499
# 3 shared_heterogeneous   2718 0.130
# 4 shared_opposite         777 0.0373
# 5 NA                     4747 0.228

dplyr::group_by(df.peak, variant_heterogeneity) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(frac = count / sum(count))
# variant_heterogeneity count   frac
# <chr>                 <int>  <dbl>
#   1 distinct_variants     29257 0.139
# 2 shared_consistent     91237 0.433
# 3 shared_heterogeneous  27842 0.132
# 4 shared_opposite        3985 0.0189
# 5 NA                    58263 0.277

p.gene = plot_gene_sankey(df.gene)
p.peak = plot_gene_sankey(df.peak)

sanky_size = 3
cowplot::save_plot(
  "figures/Fig3_gene_sankey.pdf",
  p.gene,
  base_height = sanky_size,
  base_width = sanky_size
)

cowplot::save_plot(
  "figures/Fig3_peak_sankey.pdf",
  p.peak,
  base_height = sanky_size,
  base_width = sanky_size
)

################################################################################

# plink2 --pfile FG_EA5_batch1_5 --extract boxplot_variants.extract --export A --out boxplot_variants

# ENSG00000131196 - NFATC1

munge_inv_bed = function(bed_pattern, cell_types, gene_id) {
  inv = purrr::map_dfr(cell_types, function(cell_type) {
    rgsutil::read_gsfile(sprintf(bed_pattern, cell_type)) %>%
      dplyr::filter(gene_id == .env$gene_id) %>%
      dplyr::select(-(chr:end)) %>%
      tibble::column_to_rownames("gene_id") %>%
      t() %>%
      as.data.frame() %>%
      dplyr::mutate(cell_type = .env$cell_type) %>%
      tibble::rownames_to_column("IID")
  })
  return(inv)
}
hard_call_dosage = function(dosage) {
  dplyr::case_when(
    (0 <= dosage) & (dosage <= 0.1) ~ 0,
    (0.9 <= dosage) & (dosage <= 1.1) ~ 1,
    (1.9 <= dosage) ~ 2,
    TRUE ~ NA_real_
  )
}

df.gt = rgsutil::read_gsfile("gs://expansion_areas/multiome/misc/qtl_boxplot/boxplot_variants.raw")

# d-k
inv.IKZF1 =
  munge_inv_bed(
    "gs://expansion_areas/multiome/batch1_5/pseudobulk/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.bed.gz",
    c(
      "predicted.celltype.l2.CD4_TEM",
      "predicted.celltype.l2.CD8_TEM"
    ),
    "ENSG00000185811"
  ) %>%
  dplyr::left_join(df.gt, by = "IID") %>%
  dplyr::mutate(gt = hard_call_dosage(2 - chr7_50361744_T_C_T),
                cell_type = remove_cell_type_prefix(cell_type)) %>%
  tidyr::drop_na(gt)

inv.FOXO1 =
  munge_inv_bed(
    "gs://expansion_areas/multiome/batch1_5/pseudobulk/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.bed.gz",
    c(
      "predicted.celltype.l2.CD14_Mono",
      "predicted.celltype.l2.CD4_Naive"
    ),
    "ENSG00000150907"
  ) %>%
  dplyr::left_join(df.gt, by = "IID") %>%
  dplyr::mutate(gt = hard_call_dosage(2 - chr13_40595488_A_G_A),
                cell_type = remove_cell_type_prefix(cell_type)) %>%
  tidyr::drop_na(gt)

inv.NFATC1 =
  munge_inv_bed(
    "gs://expansion_areas/multiome/batch1_5/pseudobulk/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.bed.gz",
    c("predicted.celltype.l1.Mono", "predicted.celltype.l1.B"),
    "ENSG00000131196"
  ) %>%
  dplyr::left_join(df.gt, by = "IID") %>%
  dplyr::mutate(
    gt = hard_call_dosage(2 - chr18_79405991_G_C_G),
    cell_type = stringr::str_remove(cell_type, "predicted\\.celltype\\.l1\\.")
  ) %>%
  tidyr::drop_na(gt)

inv.ZBTB1 =
  munge_inv_bed(
    "gs://expansion_areas/multiome/batch1_5/pseudobulk/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.bed.gz",
    c(
      "predicted.celltype.l2.CD14_Mono",
      "predicted.celltype.l2.CD4_TCM"
    ),
    "ENSG00000126804"
  ) %>%
  dplyr::left_join(df.gt, by = "IID") %>%
  dplyr::mutate(gt = hard_call_dosage(2 - chr14_64505031_G_A_G),
                cell_type = remove_cell_type_prefix(cell_type)) %>%
  tidyr::drop_na(gt)

inv.ZNF804A =
  munge_inv_bed(
    "gs://expansion_areas/multiome/batch1_5/pseudobulk/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.bed.gz",
    c("predicted.celltype.l1.B", "predicted.celltype.l1.Mono"),
    "ENSG00000170396"
  ) %>%
  dplyr::left_join(df.gt, by = "IID") %>%
  dplyr::mutate(gt = hard_call_dosage(2 - chr2_184597314_T_C_T),
                cell_type = remove_cell_type_prefix(cell_type)) %>%
  tidyr::drop_na(gt)

inv.ZNF7 =
  munge_inv_bed(
    "gs://expansion_areas/multiome/batch1_5/pseudobulk/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.bed.gz",
    c(
      "predicted.celltype.l2.CD14_Mono",
      "predicted.celltype.l2.CD4_Naive"
    ),
    "ENSG00000147789"
  ) %>%
  dplyr::left_join(df.gt, by = "IID") %>%
  dplyr::mutate(gt = hard_call_dosage(2 - chr8_144851323_C_T_C),
                cell_type = remove_cell_type_prefix(cell_type)) %>%
  tidyr::drop_na(gt)

plot_boxplot = function(df,
                        gene_id,
                        hide.ylab = FALSE,
                        cell_type.colors = l1.colors,
                        cell_type.labels = l1.labels) {
  ggplot(df, aes(gt, !!as.symbol(gene_id))) +
    geom_boxplot(
      aes(group = interaction(gt, cell_type), fill = cell_type),
      position = position_dodge(width = 0.8),
      linewidth = 0.25,
      outlier.size = 0.25
    ) +
    geom_smooth(
      aes(group = cell_type),
      method = lm,
      linewidth = 0.5,
      color = "grey80",
      show.legend = FALSE
    ) +
    locusviz::get_default_theme(
      hide.ylab = hide.ylab,
      legend.position = c(1, 1.25),
      legend.justification = c(1, 1.25)
    ) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(
        hjust = 0,
        margin = margin(),
        face = "italic"
      ),
      strip.text = element_blank(),
      strip.background = element_blank()
    ) +
    scale_x_continuous(breaks = seq(0, 2)) +
    scale_fill_manual(values = cell_type.colors, labels = cell_type.labels) +
    scale_color_manual(values = cell_type.colors, labels = cell_type.labels) +
    coord_cartesian(ylim = c(-3, 3)) +
    facet_grid(~ cell_type) +
    guides(fill = guide_legend(override.aes = list(color = NA, linetype = 0)))
}

# d-k
p.IKZF1 =
  plot_boxplot(
    inv.IKZF1,
    "ENSG00000185811",
    hide.ylab = FALSE,
    cell_type.colors = l2.colors,
    cell_type.labels = l2.labels
  ) +
  labs(x = "rs6976036 (ALT)", y = "Normalized expression", title = "IKZF1")

p.FOXO1 =
  plot_boxplot(
    inv.FOXO1,
    "ENSG00000150907",
    hide.ylab = TRUE,
    cell_type.colors = l2.colors,
    cell_type.labels = l2.labels
  ) +
  labs(x = "rs17061503 (ALT)", y = "Normalized expression", title = "FOXO1")

p.NFATC1 =
  plot_boxplot(inv.NFATC1, "ENSG00000131196", hide.ylab = TRUE) +
  labs(x = "rs68734 (ALT)", y = "Normalized expression", title = "NFATC1")

p.ZBTB1 =
  plot_boxplot(
    inv.ZBTB1,
    "ENSG00000126804",
    hide.ylab = FALSE,
    cell_type.colors = l2.colors,
    cell_type.labels = l2.labels
  ) +
  labs(x = "rs9323451 (ALT)", y = "Normalized expression", title = "ZBTB1")

p.ZNF804A =
  plot_boxplot(
    inv.ZNF804A,
    "ENSG00000170396",
    hide.ylab = TRUE,
    cell_type.colors = l1.colors,
    cell_type.labels = l1.labels
  ) +
  labs(x = "rs10497655 (ALT)", y = "Normalized expression", title = "ZNF804A")

p.ZNF7 =
  plot_boxplot(
    inv.ZNF7,
    "ENSG00000147789",
    hide.ylab = TRUE,
    cell_type.colors = l2.colors,
    cell_type.labels = l2.labels
  ) +
  labs(x = "rs1209879 (ALT)", y = "Normalized expression", title = "ZNF7")


plot_locuszoom_distinct = function(df.z,
                                   gene_id,
                                   highlight_pos,
                                   title = NULL,
                                   window = 50000) {
  x = dplyr::filter(df.features, phenotype_id == .env$gene_id)
  chrom = x$chrom[1]
  start = x$start[1] - 50000
  end = x$end[1] + 50000
  
  
  p.manhattan =
    ggplot(df.z, aes(position, nlog10p, color = cell_type)) +
    geom_vline(
      xintercept = highlight_pos,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = -log10(5e-8),
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.5
    ) +
    ggrastr::rasterize(geom_point(
      size = 1.5,
      alpha = 0.25,
      show.legend = TRUE
    ), dpi = 300) +
    geom_point(
      data = dplyr::filter(df.z, lead_variant),
      aes(x = position, y = nlog10p),
      shape = 18,
      size = 3
    ) +
    labs(
      x = sprintf("Chromosome %s", stringr::str_remove(chrom, "^chr")),
      y = expression(paste(-log[10], "(", italic(P), ")")),
      title = title
    ) +
    scale_x_continuous(labels = scales::label_number(scale = 1e-6, drop0trailing = TRUE)) +
    scale_y_continuous(expand = expansion(c(0, 0.1), 0)) +
    scale_color_manual(values = l1.colors, labels = l1.labels) +
    locusviz::get_default_theme(
      hide.xtitle = TRUE,
      legend.position = c(0, 1),
      legend.justification = c(0, 1)
    ) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(
        hjust = 0,
        margin = margin(),
        face = "italic"
      )
    ) +
    coord_cartesian(xlim = c(start, end))
  
  p.gene =
    locusviz::plot_gene_panel(chrom,
                              start,
                              end,
                              genome_build = "hg38",
                              highlight_pos = highlight_pos) +
    labs(x = sprintf("Chromosome %s (Mb)", stringr::str_remove(chrom, "^chr")))
  return(list(p.manhattan, p.gene))
}

# NFATC2
df.z = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/saige_qtl/step2_cis/zfiles/%s/integrated_gex_batch1_5.fgid.qc.%s.mean.inv.SAIGE.%s.%s.z",
  c("predicted.celltype.l1.B", "predicted.celltype.l1.CD4_T"),
  "ENSG00000101096"
) %>%
  dplyr::mutate(
    nlog10p = pchisq((beta / se) ** 2, 1, log.p = TRUE, lower.tail =
                       FALSE) / -log(10),
    lead_variant = dplyr::case_when(
      rsid == "chr20_51549451_T_A" &
        cell_type == "predicted.celltype.l1.CD4_T" ~ TRUE,
      rsid == "chr20_51460945_T_C" &
        cell_type == "predicted.celltype.l1.B" ~ TRUE,
      TRUE ~ FALSE
    ),
    cell_type = stringr::str_remove(cell_type, "predicted\\.celltype\\.l1\\.")
  )


p.NFATC2.locuszoom = plot_locuszoom_distinct(
  df.z,
  "ENSG00000101096",
  highlight_pos = c(51549451, 51460945),
  title = "NFATC2"
)
p.NFATC2.locuszoom[[1]] = p.NFATC2.locuszoom[[1]] + labs(tag = "g")

################################################################################

df.open4gene.sig = readRDS("data/open4gene.sig.rds")

start.APOBEC3A = 38734725
end.APOBEC3A = 39320389
highlight_pos.APOBEC3A = c(38962032)
peaks.APOBEC3A = c("chr22-38926990-38927724", "chr22-39003064-39003563")

cell_types.APOBEC3A = setdiff(
  l1.cell_types,
  c(
    "predicted.celltype.l1.PBMC",
    "predicted.celltype.l1.other",
    "predicted.celltype.l1.other_T"
  )
)
df.atac.APOBEC3A.i = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/atac_results/cis_nominal/zfiles/%s/integrated_atac_batch1_5.fgid.%s.sum.inv.%s.%s.z",
  cell_types.APOBEC3A,
  "chr22-38926990-38927724"
)
df.atac.APOBEC3A.ii = munge_zfile(
  "gs://fg-temp-saige-qtl2/batch1_5/atac_results/cis_nominal/zfiles/%s/integrated_atac_batch1_5.fgid.%s.sum.inv.%s.%s.z",
  cell_types.APOBEC3A,
  "chr22-39003064-39003563"
)

p.atac.APOBEC3A.i.forest =
  dplyr::filter(df.atac.APOBEC3A.i, rsid == "chr22_38962032_A_C") %>%
  plot_forest(hide.xtitle = TRUE)

p.atac.APOBEC3A.ii.forest =
  dplyr::filter(df.atac.APOBEC3A.ii, rsid == "chr22_38962032_A_C") %>%
  plot_forest(hide.xtitle = TRUE)

df.atac.APOBEC3A.i.CD4_T =
  dplyr::filter(df.atac.APOBEC3A.i,
                cell_type == "predicted.celltype.l1.CD4_T") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = ifelse(variant == "chr22_38962032_A_C", TRUE, FALSE)) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.atac.APOBEC3A.i.Mono =
  dplyr::filter(df.atac.APOBEC3A.i, cell_type == "predicted.celltype.l1.Mono") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = ifelse(variant == "chr22_38962032_A_C", TRUE, FALSE)) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.atac.APOBEC3A.ii.CD4_T =
  dplyr::filter(df.atac.APOBEC3A.ii,
                cell_type == "predicted.celltype.l1.CD4_T") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = ifelse(variant == "chr22_38962032_A_C", TRUE, FALSE)) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

df.atac.APOBEC3A.ii.Mono =
  dplyr::filter(df.atac.APOBEC3A.ii,
                cell_type == "predicted.celltype.l1.Mono") %>%
  dplyr::mutate(r2 = 0, pip = NA, cs_id = NA) %>%
  locusviz::preprocess(variant_col = "rsid") %>%
  dplyr::mutate(lead_variant = ifelse(variant == "chr22_38962032_A_C", TRUE, FALSE)) %>%
  locusviz::annotate_r2(reference_panel = "sisu42", window = 2e6)

peak.APOBEC3A.start = c(38926990, 39003064)
peak.APOBEC3A.end = c(38927724, 39003563)

background.layers = list(geom_rect(
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf
  ),
  fill = "grey90",
  data = data.frame(xmin = peak.APOBEC3A.start, xmax = peak.APOBEC3A.end)
))

p.atac.APOBEC3A.i.CD4_T =
  locusviz::plot_manhattan_panel(
    df.atac.APOBEC3A.i.CD4_T,
    highlight_pos = highlight_pos.APOBEC3A,
    background.layers = background.layers,
    xlim = c(start.APOBEC3A, end.APOBEC3A),
    title = "chr22:38926990-38927724 caQTL (CD4+ T)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none")

p.atac.APOBEC3A.ii.CD4_T =
  locusviz::plot_manhattan_panel(
    df.atac.APOBEC3A.ii.CD4_T,
    highlight_pos = highlight_pos.APOBEC3A,
    background.layers = background.layers,
    xlim = c(start.APOBEC3A, end.APOBEC3A),
    title = "chr22:39003064-39003563 caQTL (CD4+ T)"
  ) +
  locusviz::get_default_theme(hide.xlab = TRUE, legend.position = "none")

p.gene.APOBEC3A =
  locusviz::plot_gene_panel(
    "chr22",
    start.APOBEC3A,
    end.APOBEC3A,
    genome_build = "hg38",
    highlight_pos = highlight_pos.APOBEC3A
  )

df.link.APOBEC3A =
  dplyr::filter(
    df.open4gene.sig,
    cell_type %in% cell_types.APOBEC3A &
      peak_id %in% peaks.APOBEC3A &
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

df.atac.APOBEC3A.tracks = import_peak_track(
  "data/fragment.%s.chr22_38962032_A_C.%d.bw",
  cell_types = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.CD4_T")
)

peak.ranges.APOBEC3A = legendry::key_range_manual(start = peak.APOBEC3A.start,
                                                  end = peak.APOBEC3A.end,
                                                  name = c("i", "ii"))

list(
  plot_peak(
    df.atac.APOBEC3A.tracks,
    start = 38926990 - 250,
    end = 38927724 + 250,
    highlight_pos = highlight_pos.APOBEC3A,
    cell_type = "Mono",
    hide.xtitle = TRUE,
    hide.ytitle = FALSE
  ) +
    guides(x.sec = legendry::primitive_bracket(peak.ranges.APOBEC3A)),
  plot_peak(
    df.atac.APOBEC3A.tracks,
    start = 38926990 - 250,
    end = 38927724 + 250,
    highlight_pos = highlight_pos.APOBEC3A,
    cell_type = "CD4_T",
    hide.xtitle = TRUE,
    hide.ytitle = FALSE,
    reverse.GT = TRUE
  ) +
    guides(x.sec = legendry::primitive_bracket(peak.ranges.APOBEC3A)),
  plot_peak(
    df.atac.APOBEC3A.tracks,
    start = 39003064 - 250,
    end = 39003563 + 250,
    highlight_pos = highlight_pos.APOBEC3A,
    cell_type = "Mono",
    hide.xtitle = TRUE,
    hide.ytitle = FALSE,
    reverse.GT = TRUE
  ) +
    guides(x.sec = legendry::primitive_bracket(peak.ranges.APOBEC3A)),
  plot_peak(
    df.atac.APOBEC3A.tracks,
    start = 39003064 - 250,
    end = 39003563 + 250,
    highlight_pos = highlight_pos.APOBEC3A,
    cell_type = "CD4_T",
    hide.xtitle = TRUE,
    hide.ytitle = FALSE
  ) +
    guides(x.sec = legendry::primitive_bracket(peak.ranges.APOBEC3A))
) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 1)


p.atac.APOBEC3A.tracks.CD4_T =
  plot_peak(
    df.atac.APOBEC3A.tracks,
    start = start.APOBEC3A,
    end = end.APOBEC3A,
    highlight_pos = highlight_pos.APOBEC3A,
    cell_type = "CD4_T",
    hide.xtitle = TRUE,
    hide.ytitle = FALSE,
    rasterize = TRUE
  )

p.atac.APOBEC3A.tracks.Mono =
  plot_peak(
    df.atac.APOBEC3A.tracks,
    start = start.APOBEC3A,
    end = end.APOBEC3A,
    highlight_pos = highlight_pos.APOBEC3A,
    cell_type = "Mono",
    hide.xtitle = TRUE,
    hide.ytitle = FALSE,
    rasterize = TRUE
  )


p.link.APOBEC3A =
  plot_links(
    df.link.APOBEC3A,
    start.APOBEC3A,
    end.APOBEC3A,
    highlight_pos.APOBEC3A,
    peak.ranges.APOBEC3A,
    # xbreaks = c(27300000, 27350000, 27450000),
    background.layers = background.layers
  ) +
  labs(x = "Chromosome 22", title = "Peak-gene links")

list(# p.atac.APOBEC3A.i.CD4_T,
  # p.atac.APOBEC3A.ii.CD4_T,
  # p.atac.APOBEC3A.tracks.CD4_T,
  # p.atac.APOBEC3A.tracks.Mono,
  p.gene.APOBEC3A, p.link.APOBEC3A) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1)


#################################

# chr5-56327396-56328237
window = 500
start = 56327396 - window
# end = 56328237 + window
end = 56329186 + window
# chr5_56329186_G_A
highlight_pos = 56329186



df.track = import_peak_track(
  "data/fragment.%s.chr5_56329186_G_A.%d.bw",
  cell_types = c(
    "predicted.celltype.l1.B",
    "predicted.celltype.l1.CD4_T",
    "predicted.celltype.l1.Mono"
  )
) %>%
  dplyr::filter(start - 100 <= x & x <= end + 100)

dplyr::filter(
  df.open4gene.sig,
  peak_id %in% c("chr5-56327396-56328237", "chr5-56328850-56329349")
) %>% View()

plot_peak_fig3 = function(df.track,
                          start = NULL,
                          end = NULL,
                          highlight_pos = NULL) {
  ggplot() +
    geom_area(aes(x, score, fill = GT), position = "identity", data = df.track) +
    locusviz::or_missing(
      !is.null(highlight_pos),
      geom_vline(
        xintercept = highlight_pos,
        linetype = "dashed",
        color = "grey50"
      )
    ) +
    locusviz::or_missing(!is.null(start) &
                           !is.null(end), coord_cartesian(xlim = c(start, end)))     +
    scale_x_continuous(
      expand = expansion(),
      labels = scales::label_number(scale = 1e-6, drop0trailing = TRUE)
    ) +
    scale_y_continuous(expand = expansion(c(0, 0.05)),
                       labels = scales::label_number(drop0trailing = TRUE)) +
    locusviz::get_default_theme(legend.position = c(0, 1),
                                legend.justification = c(0, 1)) +
    theme(plot.title = element_text(hjust = 0, margin = margin()),
          legend.title = element_blank()) +
    labs(x = "Chromosome 5 (Mb)", y = "Norm. CA", fill = "rs158488")
}


peak_start = c(56327396, 56328850)
peak_end = c(56328237, 56329349)

peak.ranges <- legendry::key_range_manual(start = peak_start,
                                          end   = peak_end,
                                          name  = c("i", "ii"))

p.peak.Mono =
  dplyr::filter(df.track, cell_type == "Mono") %>%
  plot_peak_fig3(start, end, highlight_pos) +
  scale_fill_manual(values = rev(locusviz::distinct_shades(l1.colors["Mono"]))) +
  labs(title = "Mono", tag = "l") +
  guides(x.sec = legendry::primitive_bracket(peak.ranges))

p.peak.B =
  ggplot() +
  geom_area(
    aes(x, score, fill = GT),
    position = "identity",
    data = dplyr::filter(df.track, cell_type == "B" &
                           x < 56329186 - 600)
  ) +
  geom_area(
    aes(x, score, fill = GT),
    position = "identity",
    data = dplyr::filter(df.track, cell_type == "B" &
                           x >= 56329186 - 650) %>%
      dplyr::mutate(GT = factor(GT, levels = seq(2, 0)))
  ) +
  geom_vline(xintercept = highlight_pos,
             linetype = "dashed",
             color = "grey50") +
  coord_cartesian(xlim = c(start, end)) +
  
  scale_fill_manual(values = rev(locusviz::distinct_shades(l1.colors["B"]))) +
  scale_x_continuous(
    expand = expansion(),
    labels = scales::label_number(scale = 1e-6, drop0trailing = TRUE)
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.05)),
                     labels = scales::label_number(drop0trailing = TRUE)) +
  locusviz::get_default_theme() +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0, margin = margin()),
    legend.title = element_blank(),
    legend.position.inside = c(0, 1),
    legend.justification.inside = c(0, 1)
  ) +
  labs(
    x = "Chromosome 5 (Mb)",
    fill = "ALT",
    title = "B",
    tag = "m"
  ) +
  guides(x.sec = legendry::primitive_bracket(peak.ranges))

p.peak.CD4_T =
  dplyr::filter(df.track, cell_type == "CD4_T") %>%
  dplyr::mutate(GT = factor(GT, levels = seq(2, 0))) %>%
  plot_peak_fig3(start, end, highlight_pos) +
  scale_fill_manual(values = locusviz::distinct_shades(l1.colors["CD4_T"]), breaks = seq(0, 2)) +
  theme(axis.title.y = element_blank()) +
  labs(title = "CD4+ T", tag = 'n') +
  guides(x.sec = legendry::primitive_bracket(peak.ranges))

df.peak.z =
  dplyr::bind_rows(
    munge_zfile(
      "gs://fg-temp-saige-qtl2/batch1_5/atac_results/cis_nominal/zfiles/%s/integrated_atac_batch1_5.fgid.%s.sum.inv.%s.%s.z",
      l1.cell_types.primary,
      "chr5-56327396-56328237"
    ),
    munge_zfile(
      "gs://fg-temp-saige-qtl2/batch1_5/atac_results/cis_nominal/zfiles/%s/integrated_atac_batch1_5.fgid.%s.sum.inv.%s.%s.z",
      setdiff(l1.cell_types.primary, "predicted.celltype.l1.DC"),
      "chr5-56328850-56329349"
    )
  ) %>%
  dplyr::filter(rsid == "chr5_56329186_G_A")

pd = position_dodge(width = 0.5)
p.peak.forest =
  dplyr::mutate(
    df.peak.z,
    cell_type = stringr::str_remove(cell_type, "predicted\\.celltype\\.l1\\."),
    cell_type = factor(cell_type, levels = rev(
      c("CD4_T", "CD8_T", "B", "NK", "DC", "Mono")
    )),
    prefix = ifelse(gene_id == "chr5-56327396-56328237", "i. ", "ii. "),
    gene_id = stringr::str_c(prefix, stringr::str_replace(gene_id, "-", ":"))
  ) %>%
  ggplot(aes(beta, cell_type, color = gene_id, shape = gene_id)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = beta - se, xmax = beta + se),
    height = 0,
    position = pd,
    show.legend = FALSE
  ) +
  geom_point(position = pd) +
  locusviz::get_default_theme(
    hide.ytitle = TRUE,
    legend.position = c(0, 1.2),
    legend.justification = c(0, 1.2)
  ) +
  theme(legend.title = element_blank()) +
  scale_y_discrete(labels = l1.labels) +
  scale_color_manual(values = BuenColors::jdb_palette("corona")) +
  labs(
    x = "Effect size (rs158488)",
    shape = "Peak",
    color = "Peak",
    tag = "o"
  )
p.peak.forest

layout = "
ABCD
ABCE
FGHI
JKLM
"
p.box = c(
  list(
    p.IKZF1 + labs(tag = "d"),
    patchwork::free(p.FOXO1 + labs(tag = "e"), side = "l"),
    patchwork::free(p.NFATC1 + labs(tag = "f"), side = "l")
  ),
  purrr::map(p.NFATC2.locuszoom, function(p) {
    patchwork::free(p, side = "b")
  }),
  list(
    p.ZBTB1 + labs(tag = "h"),
    patchwork::free(p.ZNF804A + labs(tag = "i"), side = "l"),
    patchwork::free(p.ZNF7 + labs(tag = "j"), side = "l"),
    ggplot_spacer + labs(tag = "k"),
    p.peak.Mono,
    patchwork::free(p.peak.B, side = "l"),
    patchwork::free(p.peak.CD4_T, side = "l"),
    p.peak.forest
  )
) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(
    design = layout,
    widths = c(rep(0.6, 3), 1),
    height = c(1, 0.01, 1.01, 1.1)
  )
p.box

cowplot::save_plot(
  "figures/Fig3_boxplot.pdf",
  p.box,
  base_height = 130,
  base_width = 180,
  units = "mm"
)
################################################################################

df.motif = data.table::fread(
  "https://raw.githubusercontent.com/kundajelab/MotifCompendium/refs/heads/main/pipeline/data/MotifCompendium-Database-Human.metadata.tsv"
)
df.features.TF = dplyr::filter(df.features, symbol %in% unique(unlist(stringr::str_split(df.motif$TF, ","))))
df.features.ZNF = dplyr::filter(df.features.TF,
                                stringr::str_detect(symbol, "^ZNF|^ZBTB|^ZSCAN|^ZKSCAN|^ZFP"))

df.all.ZNF =
  rgsutil::read_gsfile(
    "gs://expansion_areas/multiome/batch1_5/coloc/susie/clusters/all_feature_combined_clusters.tsv.gz"
  ) %>%
  dplyr::filter(trait %in% df.features.ZNF$phenotype_id)
ZNF.cluster.113164 = dplyr::filter(df.all.ZNF, cluster == 113164) %>%
  dplyr::distinct(trait) %>%
  dplyr::pull()

df.gex.in_cs.ZNF = data.table::fread("data/integrated_gex_batch1_5.fgid.qc.in_cs.ZNF.tsv.gz",
                                     data.table = FALSE)
df.ZNF =
  dplyr::inner_join(df.gex.in_cs.ZNF,
                    df.features.ZNF,
                    by = c("region" = "phenotype_id")) %>%
  dplyr::filter(region %in% ZNF.cluster.113164 &
                  rsid == "chr19_21587341_C_A") %>%
  dplyr::mutate(cell_type = parse_cell_type(trait)) %>%
  dplyr::filter(
    filter_l1_cell_types(cell_type) &
      !cell_type %in% c(
        "predicted.celltype.l1.PBMC",
        "predicted.celltype.l1.other_T"
      )
  ) %>%
  dplyr::transmute(symbol = symbol, cell_type = l1.labels[remove_cell_type_prefix(cell_type)], beta)

mat =
  tidyr::pivot_wider(
    df.ZNF,
    id_cols = symbol,
    names_from = cell_type,
    values_from = beta
  ) %>%
  tibble::column_to_rownames("symbol") %>%
  as.matrix()
mat[is.na(mat)] = 0

row_clust <- hclust(dist(mat))
col_clust <- hclust(dist(t(mat)))

# Reorder matrix based on clustering
mat_clustered <- mat[row_clust$order, col_clust$order]

# Convert to long format for ggplot2
mat_long <- mat_clustered %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variant_id") %>%
  tidyr::pivot_longer(cols = -variant_id,
                      names_to = "cell_type",
                      values_to = "beta") %>%
  dplyr::mutate(
    variant_id = factor(variant_id, levels = rownames(mat_clustered)),
    cell_type = factor(cell_type, levels = colnames(mat_clustered))
  )

# Create color palette
col_fun = circlize::colorRamp2(c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)), BuenColors::jdb_palette("brewer_yes")[c(1, 5, 9)])

# Create heatmap with clustering
ht <- ComplexHeatmap::Heatmap(
  mat,
  name = "rs55778393 (ALT)",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 6, fontface = "italic"),
  row_dend_width = unit(4, "mm"),
  column_names_gp = grid::gpar(fontsize = 6),
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_dend_height = unit(4, "mm"),
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontsize = 7),
    labels_gp = grid::gpar(fontsize = 6),
    title_position = "leftcenter",
    direction = "horizontal",
    grid_height = unit(1, "mm"),
    grid_width = unit(2, "mm")
  )
)
ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")

pdf("figures/Fig3_ZNF_heatmap.pdf",
    width = 2.6,
    height = 1.55)
ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
dev.off()

