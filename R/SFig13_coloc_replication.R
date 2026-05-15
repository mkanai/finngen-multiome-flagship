library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.trait = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/metadata/phenotypes/FG_R12_endpoints.KANTA.summary.txt"
) %>%
  dplyr::filter(H2_Z > 2 & num_gw_significant > 10) %>%
  dplyr::select(-sumstats, -ldsc_munged) %>%
  dplyr::rename(
    ldsc_h2 = H2,
    ldsc_h2_se = SE,
    ldsc_h2_z = H2_Z,
    ldsc_int = INT,
    ldsc_int_se = INT_SE,
    ldsc_ratio = RATIO,
    ldsc_ratio_se = RATIO_SE
  )

df.coloc.raw = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/coloc/susie/*.eQTL.colocQC.tsv.gz", function(df, path) {
  dplyr::filter(
    df,
    low_purity1 == 0 &
      low_purity2 == 0 &
      cs_overlap > 0 &
      pmin(probmass_1, probmass_2) > 0.5
  ) %>%
    dplyr::mutate(
      cell_type = stringr::str_replace(
        dataset2,
        "^.*\\.(predicted\\.celltype\\.l[12]\\.[^\\.]*)\\..*$",
        "\\1"
      ),
      QTL = stringr::str_replace(dataset2, "^.*_([^_]*QTL)$", "\\1"),
      QTL = ifelse(dataset2 == "FIN-R12-Olink", "pQTL", QTL)
    )
}) %>%
  dplyr::filter(dataset1 %in% c("FinnGen-R12--GWAS", "FinnGen-KANTA--GWAS", "FinnGen-R12")) %>%
  dplyr::filter(trait1 %in% df.trait$phenocode)

df.coloc.tk.l1 = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/TenK10K/coloc_100kb/TenK10K_coloc.l1.tsv.bgz"
)
df.coloc.tk.l2 = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/TenK10K/coloc_100kb/TenK10K_coloc.l2.tsv.bgz"
)
df.coloc.tk.tested.l2 = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/replication/TenK10K/coloc_100kb/TenK10K_coloc.l2.tested.tsv.bgz"
)

df.coloc.tk.tested = dplyr::bind_rows(
  dplyr::transmute(
    df.coloc.tk.tested.l2,
    gene,
    finngen_phenocode,
    cell_type = predicted_celltype_l1,
    tested = TRUE
  ),
  dplyr::transmute(
    df.coloc.tk.tested.l2,
    gene,
    finngen_phenocode,
    cell_type = predicted_celltype_l2,
    tested = TRUE
  )
)  %>%
  dplyr::distinct()

df.coloc.tk =
  dplyr::bind_rows(
    dplyr::transmute(
      df.coloc.tk.l1,
      gene,
      finngen_phenocode,
      tenk10k_phenotype,
      cell_type = predicted_celltype_l1,
      PP.H4 = max_PP_H4
    ),
    dplyr::transmute(
      df.coloc.tk.l2,
      gene,
      finngen_phenocode,
      tenk10k_phenotype,
      cell_type = predicted_celltype_l2,
      PP.H4 = PP.H4.abf
    )
  ) %>%
  dplyr::right_join(df.coloc.tk.tested,
                    by = c("gene", "finngen_phenocode", "cell_type"))



matching_traits = dplyr::distinct(df.coloc.tk, finngen_phenocode) %>%
  dplyr::pull(finngen_phenocode)

pip_bin_breaks <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

df.coloc.tk.merged =
  dplyr::filter(df.coloc.raw, trait1 %in% matching_traits &
                  QTL == "eQTL") %>%
  dplyr::group_by(trait1, trait2, cell_type) %>%
  dplyr::summarize(
    max_PP_H4 = max(PP.H4.abf),
    max_PP_H4_bin = cut(max_PP_H4, pip_bin_breaks),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    df.coloc.tk,
    by = c(
      "trait1" = "finngen_phenocode",
      "trait2" = "gene",
      "cell_type" = "cell_type"
    )
  )


df.coloc.tk.rate =
  dplyr::group_by(df.coloc.tk.merged, max_PP_H4_bin) %>%
  dplyr::summarize(
    locusviz::binom_ci(sum(PP.H4 > 0.5, na.rm = TRUE), sum(tested, na.rm = TRUE), colname = "frac_05"),
    locusviz::binom_ci(sum(PP.H4 > 0.8, na.rm = TRUE), sum(tested, na.rm = TRUE), colname = "frac_08"),
    .groups = "drop"
  ) %>%
  dplyr::rename_with( ~ sub("^(frac_([0-9]+))$", "\\1_value", .x)) %>%
  tidyr::pivot_longer(
    cols = tidyselect::starts_with("frac_"),
    names_to = c("threshold", ".value"),
    names_pattern = "frac_([0-9]+)_(value|lower|upper)"
  ) %>%
  dplyr::rename(frac = value) %>%
  dplyr::mutate(
    threshold = forcats::fct_recode(threshold, "PP.H4 > 0.5" = "05", "PP.H4 > 0.8" = "08"),
    threshold = forcats::fct_rev(threshold)
  )

pd = position_dodge(width = 0.9)
sig.colors = c(
  "PP.H4 > 0.8" = BuenColors::jdb_palette("brewer_heat")[7],
  "PP.H4 > 0.5" = BuenColors::jdb_palette("brewer_heat")[5],
  "PP.H4 <= 0.5" = BuenColors::jdb_palette("brewer_heat")[3],
  "FinnGen only" = BuenColors::jdb_palette("calma_morado")[5],
  "Missing trait" = "grey50"
)

p.coloc.tk =
  ggplot(df.coloc.tk.rate, aes(frac, max_PP_H4_bin, color = threshold)) +
  geom_hline(yintercept = 3.5,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                width = 0,
                position = pd) +
  geom_point(position = pd) +
  locusviz::get_default_theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    hide.ylab = TRUE
  ) +
  theme(plot.title = element_text(hjust = 0.02)) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_color_manual(values = sig.colors) +
  guides(color = guide_legend(reverse = TRUE)) +
  coord_cartesian(xlim = c(0, 0.51)) +
  labs(x = "% replication", y = "FinnGen PP.H4 bin", color = "TenK10K")
p.coloc.tk


################################################################################
immune_tissue_ids <- c(
  # Monocyte / Macrophage
  "CL_0002057",
  # monocyte
  "CL_0002396",
  # CD16+ monocyte
  "CL_0000860",
  # classical monocyte
  "CL_0000875",
  # non-classical monocyte
  "CL_0000235",
  # macrophage
  "CL_0000583",
  # alveolar macrophage
  
  # Dendritic cells
  "CL_0000451",
  # dendritic cell (cDC/cDC2)
  "CL_0000990",
  # conventional dendritic cell (DC2)
  "CL_0000784",
  # plasmacytoid dendritic cell
  
  # NK cells
  "CL_0000623",
  # NK cell
  "CL_0000938",
  # CD56+ NK cell
  
  # Neutrophil
  "CL_0000775",
  # neutrophil
  
  # Platelet / Megakaryocyte
  "CL_0000233",
  # platelet
  "CL_0000556",
  # megakaryocyte
  "CL_0000553",
  # megakaryocyte progenitor cell
  
  # HSPC
  "CL_0008001",
  # hematopoietic precursor cell
  
  # B cells
  "CL_0000236",
  # B cell
  "CL_0000788",
  # naive B cell
  "CL_0000787",
  # memory B cell
  "CL_0000980",
  # plasmablast
  
  # CD4+ T cells
  "CL_0000624",
  # CD4+ T cell
  "CL_0000904",
  # CD4+ TCM cell
  "CL_0000905",
  # CD4+ TEM cell
  "CL_0000934",
  # CD4+ CTL cell
  "CL_0000897",
  # CD4+ memory T cell
  "CL_0000545",
  # Th1 cell
  "CL_0000546",
  # Th2 cell
  "CL_0000899",
  # Th17 cell
  "CL_0002038",
  # Tfh cell
  
  # Treg
  "CL_0002677",
  # Treg memory
  "CL_0002678",
  # Treg naive
  "CL_0000815",
  # regulatory T cell
  
  # CD8+ T cells
  "CL_0000625",
  # CD8+ T cell
  "CL_0000907",
  # CD8+ TCM cell
  "CL_0000913",
  # CD8+ TEM cell
  "CL_0001062",
  # Tem/Temra cytotoxic T cell
  
  # Other T / lymphoid
  "CL_0000084",
  # T cell (generic)
  "CL_0000940",
  # MAIT cell
  "CL_0000798",
  # gdT cell
  "CL_0002489",
  # dnT cell
  
  # Blood (whole)
  "UBERON_0000178" # blood
)


if (!file.exists("data/df.coloc.eqtl_catalogue.rds")) {
  df.eq.metadata = dplyr::bind_rows(
    data.table::fread(
      "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/refs/heads/master/data_tables/dataset_metadata_r7.tsv",
      data.table = FALSE
    ) %>%
      dplyr::mutate(release = "r7"),
    data.table::fread(
      "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/refs/heads/master/data_tables/dataset_metadata_r8_beta.tsv",
      data.table = FALSE
    ) %>%
      dplyr::mutate(release = "r8_beta")
  ) %>%
    dplyr::group_by(dataset_id) %>%
    dplyr::filter(n() == 1 | release == "r8_beta") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(is_immune = tissue_id %in% immune_tissue_ids)
  saveRDS(df.eq.metadata, "data/df.eq.metadata.rds")

  # n_egenes > 100 & eQTL (not aptamer pQTL)
  df.eq.permuted = rgsutil::read_gsfile(
    "gs://expansion_areas/multiome/misc/eQTL_Catalogue/eQTL_Catalogue.r8_beta.ge.permuted.tsv.gz"
  ) %>%
    tidyr::drop_na(molecular_trait_id)
  df.eq.n_egenes = dplyr::group_by(df.eq.permuted, dataset_id) %>%
    dplyr::summarize(n_egenes = sum(qval < 0.05), .groups = "drop")
  valid_dataset_ids = df.eq.n_egenes %>%
    dplyr::filter(n_egenes > 100) %>%
    dplyr::left_join(df.eq.metadata, by = "dataset_id") %>%
    dplyr::filter(quant_method != "aptamer") %>%
    dplyr::pull(dataset_id)
  saveRDS(valid_dataset_ids, "data/valid_dataset_ids.rds")

  df.coloc.eqtl_catalogue = rgsutil::read_gsfile(
    "gs://expansion_areas/multiome/batch1_5/coloc/susie/eQTL_Catalogue/FinnGen-R12-KANTA.eQTL_Catalogue.colocQC.tsv.gz"
  ) %>%
    dplyr::filter(
      low_purity1 == 0 &
        cs_overlap > 0 &
        pmin(probmass_1, probmass_2) > 0.5 &
        stringr::str_starts(dataset2, "eQTL_Catalogue_")
    ) %>%
    dplyr::mutate(dataset_id = sub("^eQTL_Catalogue_", "", dataset2)) %>%
    dplyr::filter(dataset_id %in% valid_dataset_ids) %>%
    dplyr::left_join(df.eq.metadata, by = "dataset_id")

  saveRDS(df.coloc.eqtl_catalogue,
          "data/df.coloc.eqtl_catalogue.rds")
}
df.coloc.eqtl_catalogue = readRDS("data/df.coloc.eqtl_catalogue.rds")

if (!file.exists("tables/ST28_eQTL_catalogue_metadata.tsv")) {
  df.eq.metadata = readRDS("data/df.eq.metadata.rds")
  valid_dataset_ids = readRDS("data/valid_dataset_ids.rds")
  
  dplyr::filter(df.eq.metadata, dataset_id %in% valid_dataset_ids) %>%
    export_table("tables/ST28_eQTL_catalogue_metadata.tsv", "ST28")
}

df.coloc.eqtl_catalogue.max_PP_H4 =
  dplyr::group_by(df.coloc.eqtl_catalogue, trait1, trait2, study_type, is_immune) %>%
  dplyr::summarize(
    tested = TRUE,
    max_PP_H4_eq = max(PP.H4.abf, na.rm = TRUE),
    .groups = "drop"
  )


df.coloc.eq.merged =
  dplyr::filter(df.coloc.raw, QTL == "eQTL") %>%
  dplyr::group_by(trait1, trait2) %>%
  dplyr::summarize(
    max_PP_H4 = max(PP.H4.abf),
    max_PP_H4_bin = cut(max_PP_H4, pip_bin_breaks),
    .groups = "drop"
  ) %>%
  dplyr::left_join(df.coloc.eqtl_catalogue.max_PP_H4, by = c("trait1", "trait2"))

study_type.colors = c(
  "Bulk (immune)" = BuenColors::jdb_palette("brewer_green")[7],
  "Bulk (non-immune)" = BuenColors::jdb_palette("brewer_green")[5],
  "Single-cell (immune)" = BuenColors::jdb_palette("brewer_purple")[7],
  "Single-cell (non-immune)" = BuenColors::jdb_palette("brewer_purple")[5]
)

df.coloc.eq.rate =
  dplyr::group_by(df.coloc.eq.merged, max_PP_H4_bin, study_type, is_immune) %>%
  dplyr::summarize(locusviz::binom_ci(sum(max_PP_H4_eq > 0.8, na.rm = TRUE), sum(tested, na.rm = TRUE)), .groups = "drop") %>%
  tidyr::drop_na() %>%
  dplyr::mutate(
    study = paste(
      stringr::str_to_sentence(study_type),
      ifelse(is_immune, "(immune)", "(non-immune)")
    ),
    study = factor(study, levels = rev(names(study_type.colors)))
  )


p.coloc.eq =
  ggplot(df.coloc.eq.rate, aes(frac, max_PP_H4_bin, color = study)) +
  geom_hline(yintercept = 3.5,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbar(aes(xmin = frac_lower, xmax = frac_upper),
                width = 0,
                position = pd) +
  geom_point(position = pd) +
  locusviz::get_default_theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1)
  ) +
  theme(plot.title = element_text(hjust = 0.02)) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_color_manual(values = study_type.colors) +
  guides(color = guide_legend(reverse = TRUE),
         shape = guide_legend(reverse = TRUE)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "% replication (PP.H4 > 0.8)", y = "FinnGen PP.H4 bin", color = "eQTL Catalogue")

df.coloc.tk.count =
  dplyr::filter(df.coloc.tk.merged, max_PP_H4 > 0.8) %>%
  dplyr::group_by(trait1, trait2) %>%
  dplyr::summarize(PP.H4 = if (all(is.na(PP.H4))) {
    NA_real_
  } else{
    max(PP.H4, na.rm = TRUE)
  }, .groups = "drop") %>%
  dplyr::mutate(
    type = dplyr::case_when(
      PP.H4 > 0.8 ~ "PP.H4 > 0.8",
      PP.H4 > 0.5 ~ "PP.H4 > 0.5",
      PP.H4 <= 0.5 ~ "PP.H4 <= 0.5",
      TRUE ~ "FinnGen only"
    )
  ) %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(
    dataset = "",
    study = "TenK10K",
    count = n(),
    .groups = "drop"
  )

n_fg_coloc =
  dplyr::filter(df.coloc.raw, PP.H4.abf > 0.8) %>%
  dplyr::distinct(trait1, trait2) %>%
  nrow()


df.coloc.tk.count.missing =
  tibble::tibble(
    type = "Missing trait",
    dataset = "",
    study = "TenK10K",
    count = n_fg_coloc - sum(df.coloc.tk.count$count)
  )


df.coloc.eq.count =
  dplyr::filter(df.coloc.eq.merged, max_PP_H4 > 0.8)  %>%
  dplyr::mutate(
    type = dplyr::case_when(
      max_PP_H4_eq > 0.8 ~ "PP.H4 > 0.8",
      max_PP_H4_eq > 0.5 ~ "PP.H4 > 0.5",
      max_PP_H4_eq <= 0.5 ~ "PP.H4 <= 0.5",
      TRUE ~ NA_character_
    )
  ) %>%
  tidyr::drop_na(type) %>%
  dplyr::group_by(study_type, is_immune, type) %>%
  dplyr::summarize(
    dataset = "eQTL Catalogue",
    study = paste(
      stringr::str_to_sentence(study_type[1]),
      ifelse(is_immune[1], "(immune)", "(non-immune)")
    ),
    count = n(),
    .groups = "drop"
  )

df.coloc.eq.count.fg_only =
  dplyr::group_by(df.coloc.eq.count, study) %>%
  dplyr::summarize(
    count = n_fg_coloc - sum(count),
    type = "FinnGen only",
    dataset = "eQTL Catalogue",
    .groups = "drop"
  )


df.coloc.count =
  dplyr::bind_rows(
    df.coloc.eq.count,
    df.coloc.eq.count.fg_only,
    df.coloc.tk.count,
    df.coloc.tk.count.missing
  ) %>%
  dplyr::mutate(
    type = factor(type, levels = rev(names(sig.colors))),
    dataset = factor(dataset, levels = c("eQTL Catalogue", "")),
    study = factor(study, levels = c("TenK10K", names(study_type.colors))),
    y = interaction(study, dataset)
  )

p.coloc.count =
  ggplot(df.coloc.count, aes(count, y, fill = type)) +
  geom_col() +
  scale_x_continuous(expand = expansion(), labels = scales::label_comma()) +
  scale_y_discrete(
    limits = rev,
    guide = legendry::guide_axis_nested(
      key = legendry::key_range_auto(sep = "\\."),
      levels_text = list(NULL, element_text(angle = 90, hjust = 0.5))
    )
  ) +
  scale_fill_manual(values = sig.colors) +
  locusviz::get_default_theme(
    legend.position = "right",
    legend.justification = "right",
    hide.ytitle = TRUE
  ) +
  theme(
    legend.margin = margin(l = -4),
    plot.tag = element_text(margin = margin(r = -8))
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "# FinnGen coloc (PP.H4 > 0.8)", fill = "Replication")


plt =
  p.coloc.count + p.coloc.eq + p.coloc.tk + patchwork::plot_layout(nrow = 1, widths = c(0.7, 1, 1)) + patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig13_coloc_replication.pdf",
  plt,
  base_height = 60,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig13_coloc_replication.png",
  plt,
  base_height = 60,
  base_width = 180,
  dpi = 300,
  units = "mm"
)

################################################################################
# ST29 - coloc novel
df.finngen_best =
  dplyr::filter(df.coloc.raw, QTL == "eQTL") %>%
  dplyr::group_by(trait1, trait2) %>%
  dplyr::slice_max(PP.H4.abf, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    trait1, trait2,
    finngen_best_cell_type = cell_type,
    finngen_max_PP_H4      = PP.H4.abf
  )

df.eq_best =
  dplyr::mutate(df.coloc.eqtl_catalogue, category = paste0(
    ifelse(study_type == "bulk", "bulk", "sc"),
    "_",
    ifelse(is_immune, "imm", "non")
  )) %>%
  dplyr::group_by(trait1, trait2, category) %>%
  dplyr::slice_max(PP.H4.abf, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(trait1, trait2, category,
                max_PP_H4 = PP.H4.abf,
                best_dataset = dataset_id,
                best_study = study_id) %>%
  tidyr::pivot_wider(
    names_from = category,
    values_from = c(max_PP_H4, best_dataset, best_study),
    names_glue = "{category}_{.value}"
  )

df.coloc.novel =
  dplyr::filter(df.finngen_best, finngen_max_PP_H4 > 0.8) %>%
  dplyr::left_join(df.eq_best, by = c("trait1", "trait2")) %>%
  dplyr::mutate(
    bulk_imm_replicated = !is.na(bulk_imm_max_PP_H4) & bulk_imm_max_PP_H4 > 0.8,
    bulk_non_replicated = !is.na(bulk_non_max_PP_H4) & bulk_non_max_PP_H4 > 0.8,
    sc_imm_replicated   = !is.na(sc_imm_max_PP_H4)   & sc_imm_max_PP_H4   > 0.8,
    sc_non_replicated   = !is.na(sc_non_max_PP_H4)   & sc_non_max_PP_H4   > 0.8,
    any_replicated      = bulk_imm_replicated | bulk_non_replicated |
      sc_imm_replicated   | sc_non_replicated,
    novel               = !any_replicated
  ) %>%
  dplyr::left_join(
     dplyr::distinct(df.features, phenotype_id, symbol),
    by = c("trait2" = "phenotype_id")
  ) %>%
  dplyr::select(
    trait1, trait2, symbol,
    finngen_best_cell_type, finngen_max_PP_H4,
    tidyselect::any_of("cs_category"),
    tidyselect::starts_with("bulk_imm_"),
    tidyselect::starts_with("bulk_non_"),
    tidyselect::starts_with("sc_imm_"),
    tidyselect::starts_with("sc_non_"),
    any_replicated, novel
  )

cat(sprintf("Rows: %d  Replicated: %d  Novel: %d\n",
            nrow(df.coloc.novel), sum(df.coloc.novel$any_replicated), sum(df.coloc.novel$novel)))

if (!file.exists("tables/ST29_coloc_novel.tsv")) {
  export_table(df.coloc.novel, "tables/ST29_coloc_novel.tsv", save_googlesheet = FALSE)
}
