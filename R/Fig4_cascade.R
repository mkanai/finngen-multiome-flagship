library(dplyr)
library(ggplot2)
library(patchwork)
library(ggsankeyfier)
source(here::here("R/const.R"))

if (!file.exists("data/df.variant.rds")) {
  df.variant =
    data.table::fread("data/cascade_results/variant_categorization.tsv.gz", data.table = FALSE) %>%
    dplyr::mutate(
      cell_type_specificity = factor(
        cell_type_specificity,
        levels = rev(cell_type_specificity_levels)
      ),
      qtl_mechanism_category = factor(qtl_mechanism_category, levels = rev(qtl_mechanism_levels)) ,
      peak_overlap = factor(
        ifelse(peak_overlap, "Peak overlap", "No overlap"),
        levels = rev(c("Peak overlap", "No overlap"))
      ),
      caqtl = factor(caqtl, levels = rev(
        c(
          "For overlapping peak",
          "For non-overlapping peak",
          "No detected caQTL"
        )
      )),
      peak_gene_link = factor(
        ifelse(
          peak_gene_link == "For overlapping peak",
          "For overlapping peak (link)",
          peak_gene_link
        ),
        levels = rev(c(
          "For caQTL peak", "For overlapping peak (link)", "No link"
        ))
      ),
      eqtl = factor(eqtl, levels = rev(
        c(
          "For linked gene",
          "For non-linked gene",
          "No detected eQTL",
          "eQTL"
        )
      ))
    )
  
  saveRDS(df.variant, "data/df.variant.rds")
}

df.variant = readRDS("data/df.variant.rds")

if (!file.exists("tables/ST20_variant_categorization.tsv")) {
  df.variant.out =
    dplyr::select(df.variant, -qtl_pattern) %>%
    dplyr::mutate(locusviz::parse_variant(variant_id, sep = "_")) %>%
    dplyr::arrange(chromosome, position) %>%
    dplyr::select(-(chromosome:alt)) %>%
    dplyr::mutate(
      best_cell_types = stringr::str_replace_all(best_cell_types, "predicted\\.celltype\\.l[12]\\.", ""),
      significant_cts = stringr::str_replace_all(significant_cts, "predicted\\.celltype\\.l[12]\\.", ""),
      cascade_gene_symbols = purrr::map_chr(stringr::str_split(cascade_peak_genes, ","), ~
                                              {
                                                stringr::str_c(gene_symbols[unique(stringr::str_split_fixed(., "->", 2)[, 2])], collapse = ",")
                                              })
    )
  
  dplyr::select(
    df.variant.out,
    -c(
      significant_cts,
      gene_affected_cell_types,
      peak_affected_cell_types,
      associated_genes,
      associated_peaks
    )
  ) %>%
    export_table("tables/ST20_variant_categorization.tsv",
                 save_googlesheet = FALSE)
}

plot_variant_sankey = function(df,
                               label.cex = 2,
                               plot.label = FALSE) {
  pos_sankey <- ggsankeyfier::position_sankey(
    width = 0.15,
    split_nodes = FALSE,
    align = "center",
    order = "as_is",
    v_space = "auto",
    h_space = "auto",
    direction = "forward"
  )
  pos_sankey_left <- ggsankeyfier::position_sankey(
    v_space = "auto",
    align = "center",
    order = "as_is",
    nudge_x = -0.05
  )
  pos_sankey_right <- ggsankeyfier::position_sankey(
    v_space = "auto",
    align = "center",
    order = "as_is",
    nudge_x = 0.05
  )
  pos_sankey_far_right <- ggsankeyfier::position_sankey(
    v_space = "auto",
    align = "center",
    order = "as_is",
    nudge_x = 0.1
  )
  
  data = dplyr::group_by(
    df,
    peak_overlap,
    caqtl,
    peak_gene_link,
    eqtl,
    qtl_mechanism_category,
    cell_type_specificity
  ) %>%
    dplyr::summarize(count = n()) %>%
    ggsankeyfier::pivot_stages_longer(
      stages_from = c(
        "peak_overlap",
        "caqtl",
        "peak_gene_link",
        "eqtl",
        "qtl_mechanism_category",
        "cell_type_specificity"
      ),
      values_from = "count"
    ) %>%
    dplyr::group_by(stage) %>%
    dplyr::mutate(
      frac_edge = count / sum(count),
      count_from = ifelse(connector == "from", count, NA),
      frac_edge_from = count_from / sum(count_from, na.rm = TRUE),
      count_stage = sum(count)
    ) %>%
    dplyr::group_by(edge_id) %>%
    dplyr::mutate(label_frac_edge_from = ifelse(frac_edge_from > 0.05 &
                                                  all(!is.na(node)), frac_edge_from, NA)) %>%
    dplyr::group_by(node, stage) %>%
    dplyr::mutate(
      count_node = sum(count),
      frac_node = count_node / dplyr::first(count_stage),
      frac_edge_to_node = count / count_node
    ) %>%
    dplyr::ungroup()
  
  dplyr::distinct(data, stage, node, count) %>%
    dplyr::group_by(stage, node) %>%
    dplyr::summarize(count = sum(count)) %>%
    dplyr::group_by(stage) %>%
    dplyr::mutate(frac = count / sum(count)) %>%
    print()
  
  dplyr::select(
    data,
    -count_from,
    -frac_edge_from,
    -count_stage,
    -label_frac_edge_from,
    -count_node
  ) %>%
    as.data.frame() %>%
    print()
  
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
    # geom_text(
    #   aes(label = scales::percent(frac_node, accuracy = 1)),
    #   stat = "sankeynode",
    #   position = pos_sankey,
    #   cex = label.cex
    # ) +
    locusviz::or_missing(
      plot.label,
      geom_text(
        aes(
          label = sprintf("%s (%s)", node, scales::percent(frac_node, accuracy = 1)),
          color = node
        ),
        stat = "sankeynode",
        position = pos_sankey_far_right,
        cex = label.cex,
        hjust = 0
      )
    )     +
    locusviz::or_missing(
      plot.label,
      geom_text(
        aes(label = scales::percent(label_frac_edge_from, accuracy = 1)),
        stat = "sankeyedge",
        position = pos_sankey_right,
        cex = label.cex * 0.8,
        color = "grey50",
        hjust = 0
      )
    ) +
    geom_text(
      aes(
        label = scales::percent(frac_node, accuracy = 1),
        color = node
      ),
      stat = "sankeynode",
      position = pos_sankey,
      # hjust = 0.25,
      cex = label.cex
    ) +
    scale_fill_manual(values = variant_sankey_colors) +
    scale_color_manual(values =  purrr::set_names(
      ifelse(
        shades::lightness(variant_sankey_colors) > 50,
        "black",
        "white"
      ),
      names(variant_sankey_colors)
    )) +
    # scale_x_discrete(expand = expansion(add = 1)) +
    # scale_y_discrete(expand = expansion(mult = 0.3)) +
    theme_void() +
    theme(legend.position = "none")
}

p.variant = plot_variant_sankey(df.variant, label.cex = 2, plot.label = FALSE)
p.variant

cowplot::save_plot(
  "figures/Fig4_variant_cascade.pdf",
  p.variant,
  base_height = 2.4,
  base_width = 7.2
)

################################################################################
# Extended

p.variant.gex = dplyr::filter(df.variant, !qtl_mechanism_category %in% c("Only caQTL (No Link)", "Only caQTL (With Link)")) %>%
  plot_variant_sankey(label.cex = 2, plot.label = FALSE)

cowplot::save_plot(
  "figures/ExtendedDataFig5_variant_cascade_gex.pdf",
  p.variant.gex,
  base_height = 2.4,
  base_width = 7.2
)

table(df.variant$qtl_mechanism_category)
with(df.variant,
     table(cell_type_specificity, qtl_mechanism_category))

table(df.variant$peak_overlap) / nrow(df.variant)
dplyr::group_by(df.variant, peak_overlap, caqtl) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::group_by(peak_overlap) %>%
  dplyr::mutate(count / sum(count), count / nrow(df.variant))

table(df.variant$peak_gene_link) / nrow(df.variant)
dplyr::mutate(df.variant, peak_gene_link = peak_gene_link %in% c("For caQTL peak", "For overlapping peak (link)")) %>%
  dplyr::group_by(peak_gene_link, caqtl) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::group_by(peak_gene_link) %>%
  dplyr::mutate(count / sum(count), count / nrow(df.variant))


p.variant.spec.heatmap =
  dplyr::mutate(
    df.variant,
    cell_type_specificity = forcats::fct_recode(cell_type_specificity, "Likely shared" = "Likely shared but underpowered")
  ) %>%
  dplyr::group_by(qtl_mechanism_category, cell_type_specificity) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::group_by(qtl_mechanism_category) %>%
  dplyr::mutate(frac = count / sum(count)) %>%
  ggplot(aes(cell_type_specificity, qtl_mechanism_category, fill = frac)) +
  geom_tile() +
  geom_text(aes(
    label = scales::percent(frac, accuracy = 1),
    color = frac < 0.2
  ), size = 2) +
  scale_fill_viridis_c(labels = scales::label_percent()) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black"),
                     guide = "none") +
  scale_x_discrete(limits = rev) +
  locusviz::get_default_theme(
    legend.position = "none",
    legend.justification = "right",
    hide.xtitle = TRUE,
    hide.ytitle = TRUE
  ) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.title = element_blank())

################################################################################
highlight_pos.BACH2 = c(90267049, 90296024)
start.BACH2 = highlight_pos.BACH2[1] - 5100
end.BACH2 = highlight_pos.BACH2[2] + 5000

peak.BACH2.start = c(90266796, 90295895)
peak.BACH2.end = c(90267690, 90296663)
peak.ranges.BACH2 = legendry::key_range_manual(
  start = peak.BACH2.start,
  end = peak.BACH2.end,
  name = c("Peak i (Enhancer)", "Peak ii (Promoter)")
)

df.open4gene.sig = readRDS("data/open4gene.sig.rds")
dplyr::filter(df.open4gene.sig, symbol == "BACH2" & peak_id == "chr6-90266796-90267690" & cell_type == "predicted.celltype.l1.CD4_T") %>% View()
dplyr::filter(df.open4gene.sig, symbol == "BACH2" & peak_id == "chr6-90295895-90296663" & cell_type == "predicted.celltype.l1.CD4_T") %>% View()

df.atac.BACH2.tracks = import_peak_track(
  "data/fragment.%s.chr6_90267049_G_A.%d.bw",
  cell_types = c("predicted.celltype.l1.CD4_T")
)

p.BACH2.rs72494581 = plot_peak(
  df.atac.BACH2.tracks,
  start = 90267049 - 400,
  end = 90267049 + 750,
  highlight_pos = highlight_pos.BACH2,
  cell_type = "CD4_T",
  hide.xtitle = TRUE,
  hide.ytitle = FALSE,
  x.breaks = c(90267000, 90267500)
) +
  geom_text(aes(x = 90267049, y = 2.1, label = "rs72494581"),
            size = 2,
            hjust = -0.1) +
  guides(x.sec = legendry::primitive_bracket(peak.ranges.BACH2))

p.BACH2.rs6908626 = plot_peak(
  df.atac.BACH2.tracks,
  start = 90296024 - 250,
  end = 90296024 + 750,
  highlight_pos = highlight_pos.BACH2,
  cell_type = "CD4_T",
  hide.xtitle = TRUE,
  hide.ytitle = TRUE,
  reverse.GT = TRUE,
  x.breaks = c(90296000, 90296500)
) +
  geom_text(aes(x = 90296024, y = 3.4, label = "rs6908626"),
            size = 2,
            hjust = -0.1) +
  guides(x.sec = legendry::primitive_bracket(peak.ranges.BACH2))

p.gene.BACH2 =
  locusviz::plot_gene_panel("chr6",
                            start.BACH2,
                            end.BACH2,
                            genome_build = "hg38",
                            highlight_pos = highlight_pos.BACH2) +
  scale_x_continuous(labels = scales::label_comma()) +
  theme(axis.text.x = element_text(size = 6.4),
        plot.margin = margin(t = 6))

layout = "
ABD
CCD
"

p.ext.cascade =
  list(
    p.BACH2.rs72494581,
    p.BACH2.rs6908626,
    p.gene.BACH2,
    patchwork::free(p.variant.spec.heatmap, side = "tb")
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(design = layout, heights = c(1, 0.05))

cowplot::save_plot(
  "figures/ExtendedDataFig5_BACH2.pdf",
  p.ext.cascade,
  base_height = 33,
  base_width = 180,
  units = "mm"
)
