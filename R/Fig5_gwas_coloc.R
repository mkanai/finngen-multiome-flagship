library(dplyr)
library(ggplot2)
source(here::here("R/const.R"))

################################################################################

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

if (!file.exists("tables/ST21_gwas_summary.tsv")) {
  export_table(df.trait, "tables/ST21_gwas_summary.tsv", "ST21")
}

################################################################################
if (!file.exists("data/df.coloc.all.rds")) {
  df.coloc.all = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/coloc/susie/*.colocQC.tsv.gz", function(df, path) {
    dplyr::filter(
      df,
      low_purity1 == 0 &
        low_purity2 == 0 &
        cs_overlap > 0 &
        PP.H4.abf > 0.8 &
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
    dplyr::filter(dataset1 %in% c("FinnGen-R12--GWAS", "FinnGen-KANTA--GWAS", "FinnGen-R12"))
  saveRDS(df.coloc.all, "data/df.coloc.all.rds")
  
  df.coloc.pair =
    readRDS("data/df.coloc.all.rds") %>%
    dplyr::filter(trait1 %in% df.trait$phenocode) %>%
    dplyr::distinct(cell_type, trait1, trait2)
  rgsutil::write_gsfile(
    df.coloc.pair,
    "gs://expansion_areas/multiome/batch1_5/coloc/susie/FinnGen-R12-KANTA.caQTL_eQTL.coloc_trait_pairs.tsv.gz",
    overwrite = TRUE
  )
}

if (!file.exists("data/df.smr.all.rds")) {
  df.smr.all = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/smr/results/coloc/*.coloc.smr.gz", function(df, path) {
    dplyr::mutate(df, QTL = ifelse(stringr::str_detect(path, "atac"), "caQTL", "eQTL"))
  }) %>%
    dplyr::mutate(sig = p_SMR < (0.05 / n()) & p_HEIDI > 0.05)
  
  saveRDS(df.smr, "data/df.smr.all.rds")
}

df.smr = readRDS("data/df.smr.all.rds") %>%
  dplyr::filter(sig)

df.coloc = readRDS("data/df.coloc.all.rds") %>%
  dplyr::filter(trait1 %in% df.trait$phenocode) %>%
  dplyr::left_join(
    df.smr,
    by = c(
      "cell_type" = "cell_type",
      "trait1" = "phenotype",
      "trait2" = "probeID",
      "QTL" = "QTL"
    )
  )

df.coloc.gwas_only = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/coloc/GWAS_only/FinnGen-R12-KANTA.colocQC.tsv.gz"
) %>%
  dplyr::filter(
    low_purity1 == 0 &
      low_purity2 == 0 &
      cs_overlap > 0 &
      PP.H4.abf > 0.8 &
      pmin(probmass_1, probmass_2) > 0.5
  ) %>%
  dplyr::filter(trait1 %in% df.trait$phenocode &
                  trait2 %in% df.trait$phenocode)

if (FALSE) {
  df.all_pairs = rgsutil::read_gsfile(
    "gs://expansion_areas/multiome/batch1_5/coloc/susie/FinnGen-R12-KANTA.caQTL_eQTL.all_pairs.tsv.gz"
  ) %>%
    dplyr::mutate(
      df.all_pairs,
      QTL = dplyr::case_when(
        stringr::str_detect(dataset2, "_atac_") ~ "caQTL",
        stringr::str_detect(dataset2, "_gex_") ~ "eQTL",
        TRUE ~ NA_character_
      )
    )
  
  df.all_pairs.n_tested =
    dplyr::group_by(df.all_pairs, QTL, trait1, region1) %>%
    dplyr::summarize(coloc_tested = TRUE, n_tested = length(unique(trait2))) %>%
    dplyr::ungroup()
  rgsutil::write_gsfile(
    df.all_pairs.n_tested,
    "gs://expansion_areas/multiome/batch1_5/coloc/susie/FinnGen-R12-KANTA.caQTL_eQTL.all_pairs.n_tested.tsv.gz"
  )
  
  df.all_pairs.n_tested.qtl =
    dplyr::group_by(df.all_pairs, QTL, trait2, region2) %>%
    dplyr::summarize(coloc_tested = TRUE, n_tested = length(unique(trait1))) %>%
    dplyr::ungroup()
  rgsutil::write_gsfile(
    df.all_pairs.n_tested.qtl,
    "gs://expansion_areas/multiome/batch1_5/coloc/susie/FinnGen-R12-KANTA.caQTL_eQTL.all_pairs.n_tested.qtl.tsv.gz"
  )
}

df.all_pairs.n_tested = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/coloc/susie/FinnGen-R12-KANTA.caQTL_eQTL.all_pairs.n_tested.tsv.gz"
)
df.all_pairs.n_tested.qtl = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/coloc/susie/FinnGen-R12-KANTA.caQTL_eQTL.all_pairs.n_tested.qtl.tsv.gz"
)

if (!file.exists("data/df.gwas.cred.all.rds")) {
  df.gwas.cred.all =
    rgsutil::map_dfr_gsfiles(c(
      "gs://finngen-production-library-green/finngen_R12/finngen_R12_analysis_data/finemap/summary/*.SUSIE.cred.summary.tsv",
      "gs://finngen-production-library-green/lab_values/gwas_release_2025_3_10/finemap/summary/*.SUSIE.cred.summary.tsv"
    ), function(df, path) {
      dplyr::filter(df, good_cs)
    })
  saveRDS(df.gwas.cred.all, "data/df.gwas.cred.all.rds")
}

df.gwas.cred = readRDS("data/df.gwas.cred.all.rds") %>%
  dplyr::filter(trait %in% df.trait$phenocode)

df.gwas.cred.n_cs = dplyr::group_by(df.gwas.cred, trait) %>%
  dplyr::summarize(n_cs = n())

dplyr::distinct(df.gwas.cred, trait, region, cs) %>%
  nrow()
# 31983

dplyr::filter(df.coloc, QTL %in% c("caQTL", "eQTL")) %>%
  dplyr::distinct(trait1, region1, cs1) %>%
  nrow()
# 15947

dplyr::group_by(df.coloc, trait1, region1, cs1) %>%
  dplyr::summarize(
    caQTL = any(QTL == "caQTL"),
    eQTL = any(QTL == "eQTL"),
    pQTL = any(QTL == "pQTL")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::summarize(caQTL = mean(caQTL),
                   eQTL = mean(eQTL),
                   pQTL = mean(pQTL))
# caQTL  eQTL   pQTL
# <dbl> <dbl>  <dbl>
#   1 0.930 0.435 0.0456


df.annot = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/metadata/annotations/vep.most_severe.consequence.tsv.bgz"
)

df.gex.max_pip = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/saige_qtl/susie/full.max_pip.annot.most_severe.txt.bgz"
)
df.atac.max_pip = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/atac_results/susie/full.max_pip.annot.most_severe.txt.bgz"
)
df.pqtl.max_pip = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/olink/Olink_pQTL_2023_10_11.SUSIE.full.max_pip.tsv.bgz"
)
df.gwas.max_pip = rgsutil::read_gsfile("gs://expansion_areas/multiome/misc/finemap_r12/SUSIE.snp.max_pip.tsv.bgz")

df.any.coloc =
  dplyr::group_by(df.coloc, hit2) %>%
  dplyr::summarize(
    eQTL = any(QTL == "eQTL"),
    caQTL = any(QTL == "caQTL"),
    pQTL = any(QTL == "pQTL"),
    any_coloc = TRUE,
    eQTL_smr =  any(QTL == "eQTL" & !is.na(p_SMR)),
    caQTL_smr =  any(QTL == "caQTL" & !is.na(p_SMR)),
    any_smr = eQTL_smr | caQTL_smr,
    gwas_traits = stringr::str_c(sort(unique(trait1)), collapse = ",")
  )

df.gwas.snp = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/misc/finemap_r12/FinnGen_R12_KANTA.SUSIE.snp.filter.tsv.gz"
) %>%
  dplyr::left_join(dplyr::select(df.annot, variant, consequence))

df.tested = dplyr::left_join(
  df.all_pairs.n_tested.qtl,
  dplyr::group_by(df.coloc, trait2, region2) %>%
    dplyr::summarize(coloc = TRUE)
) %>%
  tidyr::replace_na(list(coloc = FALSE))

dplyr::group_by(df.tested, QTL) %>%
  dplyr::summarize(
    frac = mean(coloc),
    n_coloc = sum(coloc),
    n_total = length(unique(trait2))
  )
# QTL    frac n_coloc n_total
# <chr> <dbl>   <int>   <int>
# 1 caQTL 0.110   23025  208659
# 2 eQTL  0.188    3868   20598

df.indep = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/grm/FG_EA5_batch1_5_GRM_LD_0.2_withX.extract",
  header = FALSE
) %>% dplyr::rename(variant_id = V1) %>%
  dplyr::left_join(df.annot, by = c("variant_id" = "variant")) %>%
  dplyr::left_join(dplyr::transmute(
    df.gex.max_pip,
    variant_id = variant,
    gex_max_pip = max_pip
  )) %>%
  dplyr::left_join(dplyr::transmute(
    df.atac.max_pip,
    variant_id = variant,
    atac_max_pip = max_pip
  ))

df.null =
  dplyr::mutate(
    df.indep,
    gex_max_pip = tidyr::replace_na(gex_max_pip, 0),
    atac_max_pip = tidyr::replace_na(atac_max_pip, 0)
  ) %>%
  dplyr::filter(gex_max_pip < 0.01 & atac_max_pip < 0.01) %>%
  dplyr::left_join(
    dplyr::reframe(df.coloc.gwas_only, variant_id = c(hit1, hit2)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(any_coloc = TRUE)
  ) %>%
  dplyr::mutate(
    any_coloc = !is.na(any_coloc),
    qtl_mechanism_category = factor("No molQTL")
  )

df.variant = readRDS("data/df.variant.rds")

df.variant.coloc = dplyr::left_join(df.variant, df.any.coloc, by = c("variant_id" = "hit2")) %>%
  dplyr::mutate(any_coloc = !is.na(any_coloc),
                any_smr = !is.na(any_smr) & any_smr)

df.variant.pqtl =
  dplyr::filter(df.pqtl.max_pip, max_pip > 0.5) %>%
  dplyr::left_join(
    dplyr::filter(df.any.coloc, pQTL) %>% dplyr::transmute(variant = hit2, any_coloc = TRUE),
    by = "variant"
  ) %>%
  dplyr::mutate(
    any_coloc = !is.na(any_coloc),
    qtl_mechanism_category = factor("Any pQTL")
  )
df.variant.pqtl.only =
  dplyr::filter(df.variant.pqtl, !(variant %in% df.variant$variant_id)) %>%
  dplyr::mutate(qtl_mechanism_category = factor("Only pQTL"))




df.variant2 = dplyr::bind_rows(df.variant.coloc,
                               df.null,
                               df.variant.pqtl,
                               df.variant.pqtl.only)

df.cascade.coloc = dplyr::group_split(df.variant2, qtl_mechanism_category) %>%
  purrr::map_dfr(function(data) {
    n_total <- nrow(data)
    n_coloc <- sum(data$any_coloc)
    n_smr <- sum(data$any_smr)
    
    # Return tibble with all results
    tibble::tibble(
      qtl_mechanism_category = data$qtl_mechanism_category[1],
      n_coloc = n_coloc,
      n_smr = n_smr,
      n_total = n_total,
      locusviz::binom_ci(n_coloc, n_total, colname = "frac_coloc"),
      locusviz::binom_ci(n_smr, n_total, colname = "frac_smr"),
    )
  })

if (!file.exists("tables/ST22_cascade_coloc.tsv")) {
  dplyr::mutate(df.cascade.coloc,
                qtl_mechanism_category = forcats::fct_rev(qtl_mechanism_category)) %>%
    dplyr::arrange(qtl_mechanism_category) %>%
    export_table("tables/ST22_cascade_coloc.tsv", "ST22")
}

if (!file.exists("tables/ST23_cascade_coloc_variants.tsv")) {
  dplyr::filter(
    df.variant.coloc,
    any_coloc &
      qtl_mechanism_category %in% c(
        "Local Cascade",
        "Positional Cascade",
        "Distal Cascade",
        "caQTL + eQTL (No Link)"
      )
  ) %>%
    dplyr::select(-qtl_pattern, -(eQTL:any_coloc)) %>%
    dplyr::mutate(locusviz::parse_variant(variant_id, sep = "_")) %>%
    dplyr::arrange(dplyr::desc(qtl_mechanism_category), chromosome, position) %>%
    dplyr::select(-(chromosome:alt), -significant_cts) %>%
    dplyr::mutate(
      best_cell_types = stringr::str_replace_all(best_cell_types, "predicted\\.celltype\\.", ""),
      gene_affected_cell_types = stringr::str_replace_all(gene_affected_cell_types, "predicted\\.celltype\\.", ""),
      peak_affected_cell_types = stringr::str_replace_all(peak_affected_cell_types, "predicted\\.celltype\\.", ""),
      cascade_gene_symbols = purrr::map_chr(stringr::str_split(cascade_peak_genes, ","), ~
                                              {
                                                stringr::str_c(gene_symbols[unique(stringr::str_split_fixed(., "->", 2)[, 2])], collapse = ",")
                                              })
    ) %>%
    export_table("tables/ST23_cascade_coloc_variants.tsv", "ST23")
}

if (!file.exists("tables/ST24_smr.tsv")) {
  dplyr::filter(df.coloc, !is.na(p_SMR)) %>%
    dplyr::select(
      -dataset1,
      -dataset2,
      -low_purity1,
      -low_purity2,
      -(cs1_size:probmass_2),
      -colocRes,
      -(ProbeChr:p_eQTL),
      -sig
    ) %>%
    dplyr::select(QTL, cell_type, tidyselect::everything()) %>%
    export_table("tables/ST24_coloc_smr.tsv", save_googlesheet = FALSE)
}

p.cascade.coloc =
  dplyr::mutate(df.cascade.coloc,
                y = interaction(
                  scales::comma(n_total),
                  scales::comma(n_coloc),
                  qtl_mechanism_category
                )) %>%
  ggplot(aes(frac_coloc, y, color = qtl_mechanism_category)) +
  geom_hline(
    yintercept = c(1.5, 3.5),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point() +
  geom_errorbarh(aes(xmin = frac_coloc_lower, xmax = frac_coloc_upper), height = 0) +
  locusviz::get_default_theme(hide.ytitle = TRUE, legend.position = "none") +
  scale_x_continuous(labels = scales::label_percent(drop0trailing = TRUE)) +
  scale_color_manual(values = qtl_mechanism_colors) +
  guides(y = legendry::guide_axis_nested(
    key = legendry::key_range_auto(sep = "\\."),
    levels_text = list(
      element_text(color = "grey50"),
      element_text(color = "grey50"),
      element_text()
    )
  )) +
  labs(x = "% GWAS coloc") +
  coord_cartesian(xlim = c(0, 0.13))
p.cascade.coloc

cowplot::save_plot(
  "local/figures/cascade_coloc_slide.pdf",
  p.cascade.coloc,
  base_height = 2,
  base_width = 3.6
)
cowplot::save_plot(
  "local/figures/cascade_coloc_slide.png",
  p.cascade.coloc,
  base_height = 2,
  base_width = 3.6,
  dpi = 300
)

################################################################################

df.open4gene.sig = readRDS("data/open4gene.sig.rds")
df.open4gene.sig.linked_genes =
  dplyr::select(df.open4gene.sig, cell_type, peak_id, symbol) %>%
  dplyr::distinct()  %>%
  dplyr::group_by(cell_type, peak_id) %>%
  dplyr::summarize(symbol = stringr::str_c(symbol, collapse = ","),
                   .groups = "drop") %>%
  dplyr::rename(caQTL_linked_symbol = symbol)

assign_cs_category = function(trait, df.coloc, df.gwas.snp) {
  data = dplyr::filter(df.gwas.snp, trait == .env$trait)
  
  df.coloc.trait = dplyr::filter(df.coloc, trait1 == .env$trait) %>%
    dplyr::left_join(
      dplyr::select(df.features, phenotype_id, symbol) %>%
        dplyr::rename(eQTL_symbol = symbol),
      by = c("trait2" = "phenotype_id")
    ) %>%
    dplyr::left_join(df.open4gene.sig.linked_genes,
                     by = c("cell_type", "trait2" = "peak_id"))
  
  dplyr::group_by(data, trait, region, cs) %>%
    dplyr::summarize(
      n_cs_variants = n(),
      pLoF = any(consequence == "pLoF"),
      missense = any(consequence == "Missense"),
      nonsyn_coding = pLoF | missense,
      max_pip = max(cs_specific_prob),
      max_pip_variant = variant[which.max(cs_specific_prob)],
      max_pip_most_severe = most_severe[which.max(cs_specific_prob)],
      max_pip_gene_most_severe = gene_most_severe[which.max(cs_specific_prob)],
      nonsyn_coding_genes = stringr::str_c(sort(unique(
        gene_most_severe[consequence %in% c("pLoF", "Missense")]
      )), collapse = ","),
      .groups = "drop"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      dplyr::group_by(df.coloc.trait, region1, cs1) %>%
        dplyr::summarize(
          eQTL_coloc_genes = stringr::str_c(sort(unique(eQTL_symbol)), collapse = ","),
          caQTL_linked_genes = stringr::str_c(sort(unique(
            purrr::flatten_chr(stringr::str_split(caQTL_linked_symbol, ","))
          )), collapse = ","),
          caQTL_coloc = any(QTL == "caQTL"),
          caQTL_link = caQTL_coloc &
            !is.na(caQTL_linked_genes) & caQTL_linked_genes != "",
          eQTL_coloc = any(QTL == "eQTL"),
          pQTL_coloc = any(QTL == "pQTL"),
          coloc_any = caQTL_coloc | eQTL_coloc,
          coloc_both = caQTL_coloc & eQTL_coloc,
          eQTL_smr =  any(QTL == "eQTL" & !is.na(p_SMR)),
          caQTL_smr =  any(QTL == "caQTL" & !is.na(p_SMR)),
          smr_any = eQTL_smr | caQTL_smr,
          smr_both = eQTL_smr & caQTL_smr,
          max_coloc_prob = max(PP.H4.abf, na.rm = T),
          .groups = "drop"
        ),
      by = c("region" = "region1", "cs" = "cs1")
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.logical), ~ tidyr::replace_na(.x, FALSE))) %>%
    dplyr::mutate(
      cs_category = dplyr::case_when(
        pLoF ~ "pLoF | Missense",
        missense ~ "pLoF | Missense",
        coloc_both ~ "caQTL + eQTL",
        caQTL_link ~ "caQTL (with link)",
        caQTL_coloc ~ "caQTL (no link)",
        eQTL_coloc ~ "eQTL",
        TRUE ~ "Unknown"
      ),
      cs_category = factor(cs_category, levels = rev(
        c(
          "caQTL + eQTL",
          "caQTL (with link)",
          "caQTL (no link)",
          "eQTL",
          "pLoF | Missense",
          "Unknown"
        )
      )),
      cs_category_with_pqtl = dplyr::case_when(
        pLoF ~ "pLoF | Missense",
        missense ~ "pLoF | Missense",
        caQTL_coloc &
          eQTL_coloc & pQTL_coloc ~ "caQTL + eQTL + pQTL",
        caQTL_coloc & eQTL_coloc ~ "caQTL + eQTL",
        caQTL_coloc & pQTL_coloc ~ "caQTL + pQTL",
        eQTL_coloc & pQTL_coloc ~ "eQTL + pQTL",
        caQTL_link ~ "caQTL (with link)",
        caQTL_coloc ~ "caQTL (no link)",
        eQTL_coloc ~ "eQTL",
        pQTL_coloc ~ "pQTL",
        TRUE ~ "Unknown"
      ),
      cs_category_with_pqtl = factor(cs_category_with_pqtl, levels = rev(
        c(
          "caQTL + eQTL + pQTL",
          "caQTL + eQTL",
          "caQTL + pQTL",
          "eQTL + pQTL",
          "caQTL (with link)",
          "caQTL (no link)",
          "eQTL",
          "pQTL",
          "pLoF | Missense",
          "Unknown"
        )
      ))
    )
}

if (!file.exists("data/df.all.cs_category.rds")) {
  df.all.cs_category = purrr::map_dfr(unique(df.trait$phenocode), function(trait) {
    assign_cs_category(trait, df.coloc, df.gwas.snp)
  })
  saveRDS(df.all.cs_category, "data/df.all.cs_category.rds")
}

df.all.cs_category = readRDS("data/df.all.cs_category.rds")
if (!file.exists("tables/ST26_all_cs_category.tsv")) {
  export_table(df.all.cs_category,
               "tables/ST26_all_cs_category.tsv",
               save_googlesheet = FALSE)
}

df.aiht.cs_category = dplyr::filter(df.all.cs_category, trait == "E4_HYTHY_AI_STRICT")
df.skc.cs_category = dplyr::filter(df.all.cs_category, trait == "C3_SKIN_EXALLC")
df.t2d.cs_category = dplyr::filter(df.all.cs_category, trait == "T2D")
df.height.cs_category = dplyr::filter(df.all.cs_category, trait == "HEIGHT_IRN")

df.atopic.cs_category = dplyr::filter(df.all.cs_category, trait == "L12_ATOPIC")
df.asthma.cs_category = dplyr::filter(df.all.cs_category, trait == "J10_ASTHMA_EXMORE")
df.copd.cs_category = dplyr::filter(df.all.cs_category, trait == "J10_COPD")
df.asthma_copd.cs_category = dplyr::filter(df.all.cs_category, trait == "J10_ASTHMACOPDKELA")
df.ad.cs_category = dplyr::filter(df.all.cs_category, trait == "G6_AD_WIDE")

if (!file.exists("tables/ST25_AIHT_cs_category.tsv")) {
  export_table(df.aiht.cs_category,
               "tables/ST25_AIHT_cs_category.tsv",
               "ST25")
}

df.cs_category = dplyr::bind_rows(df.aiht.cs_category, df.skc.cs_category, df.t2d.cs_category)

trait_description = c(
  "E4_HYTHY_AI_STRICT" = "Autoimmune hypothyroidism",
  "C3_SKIN_EXALLC" = "Skin cancer",
  "T2D" = "Type 2 diabetes"
)

coloc.qtl.colors = c(
  "caQTL + eQTL" = BuenColors::jdb_palette("corona")[16],
  as.character(shades::opacity(unname(
    qtl_mechanism_colors["Full Cascade"]
  ), 0.5)),
  "caQTL (with link)" = BuenColors::jdb_palette("corona")[21],
  as.character(shades::opacity(unname(
    qtl_mechanism_colors["Only caQTL (With Link)"]
  ), 0.5)),
  "caQTL (no link)" = BuenColors::jdb_palette("corona")[19],
  as.character(shades::opacity(unname(
    qtl_mechanism_colors["Only caQTL (No Link)"]
  ), 0.5)),
  "eQTL" = BuenColors::jdb_palette("corona")[18],
  unname(qtl_mechanism_colors["Only eQTL"]),
  "pLoF | Missense" = BuenColors::jdb_palette("corona")[17],
  "Unknown" = "grey80"
)

bar_titles =
  dplyr::group_by(df.cs_category, trait) %>%
  dplyr::summarize(title = sprintf("%s (%d CSs)", trait_description[trait[1]], n())) %>%
  tibble::deframe()

df.coloc.bar =
  dplyr::mutate(df.cs_category, trait = bar_titles[trait]) %>%
  dplyr::group_by(trait, cs_category) %>%
  dplyr::summarize(count = n(), .groups = "drop") %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(frac = count / sum(count))
p.coloc.bar =
  ggplot(df.coloc.bar, aes(count, trait, fill = cs_category)) +
  geom_col(position = position_stack()) +
  geom_text(aes(label = sprintf(
    "%s\n%s", cs_category, scales::percent(frac, accuracy = 0.1)
  )),
  position = position_stack(vjust = 0.5),
  size = 2) +
  scale_fill_manual(values =     coloc.qtl.colors) +
  scale_x_continuous(expand = expansion()) +
  scale_y_discrete(expand = expansion()) +
  guides(fill = guide_legend(reverse = TRUE)) +
  locusviz::get_default_theme(legend.position = "none", hide.ylab = TRUE) +
  theme(
    plot.title = element_text(margin = margin(), hjust = 0),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.4, hjust = 0),
    panel.spacing = unit(0, "mm")
  ) +
  labs(x = "# 95% CS") +
  facet_wrap(trait ~ ., scale = "free", ncol = 1)
p.coloc.bar

p.coloc.bar.with_pqtl =
  dplyr::mutate(df.cs_category, trait = bar_titles[trait]) %>%
  dplyr::group_by(trait, cs_category_with_pqtl) %>%
  dplyr::summarize(count = n(), .groups = "drop") %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(frac = count / sum(count)) %>%
  ggplot(aes(count, trait, fill = cs_category_with_pqtl)) +
  geom_col(position = position_stack()) +
  geom_text(aes(label = sprintf(
    "%s\n%s",
    cs_category_with_pqtl,
    scales::percent(frac, accuracy = 0.1)
  )),
  position = position_stack(vjust = 0.5),
  size = 2) +
  scale_fill_manual(values = coloc.qtl.colors) +
  scale_x_continuous(expand = expansion()) +
  scale_y_discrete(expand = expansion()) +
  guides(fill = guide_legend(reverse = TRUE)) +
  locusviz::get_default_theme(legend.position = "none", hide.ylab = TRUE) +
  theme(
    plot.title = element_text(margin = margin(), hjust = 0),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.4, hjust = 0),
    panel.spacing = unit(0, "mm")
  ) +
  labs(x = "# 95% CS") +
  facet_wrap(trait ~ ., scale = "free", ncol = 1)
p.coloc.bar.with_pqtl
################################################################################
df.n_coloc_cs.eQTL =
  dplyr::filter(df.coloc, QTL == "eQTL") %>%
  dplyr::distinct(trait1, region1, cs1) %>%
  dplyr::group_by(trait1) %>%
  dplyr::summarize(n_coloc_cs = n(), .groups = "drop") %>%
  dplyr::left_join(df.gwas.cred.n_cs, by = c("trait1" = "trait")) %>%
  dplyr::filter(n_cs > 10)

df.n_coloc_cs.caQTL =
  dplyr::filter(df.coloc, QTL == "caQTL") %>%
  dplyr::distinct(trait1, region1, cs1) %>%
  dplyr::group_by(trait1) %>%
  dplyr::summarize(n_coloc_cs = n(), .groups = "drop") %>%
  dplyr::left_join(df.gwas.cred.n_cs, by = c("trait1" = "trait")) %>%
  dplyr::filter(n_cs > 10)

df.n_coloc_cs.pQTL =
  dplyr::filter(df.coloc, QTL == "pQTL") %>%
  dplyr::distinct(trait1, region1, cs1) %>%
  dplyr::group_by(trait1) %>%
  dplyr::summarize(n_coloc_cs = n(), .groups = "drop") %>%
  dplyr::left_join(df.gwas.cred.n_cs, by = c("trait1" = "trait")) %>%
  dplyr::filter(n_cs > 10)

df.n_coloc_cs.caQTL_eQTL =
  dplyr::distinct(df.coloc, trait1, region1, cs1) %>%
  dplyr::group_by(trait1) %>%
  dplyr::summarize(n_coloc_cs = n(), .groups = "drop") %>%
  dplyr::left_join(df.gwas.cred.n_cs, by = c("trait1" = "trait")) %>%
  dplyr::filter(n_cs > 10)

plot_coloc_slope = function(df) {
  dplyr::mutate(
    df,
    trait_category = dplyr::case_when(
      trait1 %in% immune_traits ~ "Immune",
      trait1 %in% cancer_traits ~ "Cancer",
      trait1 %in% cardiometabolic_traits ~ "Cardiometabolic",
      TRUE ~ NA_character_
    ),
    trait_category = factor(
      trait_category,
      levels = c("Immune", "Cancer", "Cardiometabolic")
    )
  ) %>%
    dplyr::filter(!is.na(trait_category)) %>%
    ggplot(aes(n_cs, n_coloc_cs, color = trait_category)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_smooth(method = "lm", show.legend = FALSE) +
    geom_point() +
    locusviz::get_default_theme() +
    theme(
      legend.position.inside = c(1, 0),
      legend.justification.inside = c(1, 0)
    ) +
    labs(x = "# 95% CS", y = "# colocalized 95% CS", color = "Category") +
    scale_color_manual(values = BuenColors::jdb_palette("corona"))
}

df.coloc.rate =
  dplyr::bind_rows(
    df.n_coloc_cs.eQTL %>% dplyr::mutate(QTL = "eQTL"),
    df.n_coloc_cs.caQTL %>% dplyr::mutate(QTL = "caQTL"),
    df.n_coloc_cs.pQTL %>% dplyr::mutate(QTL = "pQTL")
  ) %>%
  dplyr::mutate(
    trait_category = dplyr::case_when(
      trait1 %in% immune_traits ~ "Immune",
      trait1 %in% cancer_traits ~ "Cancer",
      trait1 %in% cardiometabolic_traits ~ "Cardiometabolic",
      trait1 %in% blood_counts ~ "Blood count",
      trait1 %in% lipids_glucose ~ "Lipid / Glucose",
      TRUE ~ NA_character_
    ),
    trait_category = factor(
      trait_category,
      levels = c(
        "Immune",
        "Blood count",
        "Cancer",
        "Cardiometabolic",
        "Lipid / Glucose"
      )
    )
  ) %>%
  dplyr::filter(!is.na(trait_category))

dplyr::bind_rows(
  df.n_coloc_cs.eQTL %>% dplyr::mutate(QTL = "eQTL"),
  df.n_coloc_cs.caQTL %>% dplyr::mutate(QTL = "caQTL"),
) %>%
  dplyr::group_split(QTL) %>%
  purrr::map_dfr(function(data) {
    n_total <- sum(data$n_cs)
    n_coloc <- sum(data$n_coloc_cs)
    frac_coloc <- n_coloc / n_total
    
    ci_result <- binom::binom.confint(
      x = n_coloc,
      n = n_total,
      conf.level = 0.95,
      methods = "wilson"
    )
    
    # Return tibble with all results
    tibble::tibble(
      QTL = data$QTL[1],
      n_coloc = n_coloc,
      n_total = n_total,
      frac_coloc = frac_coloc,
      lower = ci_result$lower,
      upper = ci_result$upper
    )
  })

p.coloc.rate =
  dplyr::bind_rows(
    df.n_coloc_cs.eQTL %>% dplyr::mutate(QTL = "eQTL"),
    df.n_coloc_cs.caQTL %>% dplyr::mutate(QTL = "caQTL")
  ) %>%
  ggplot(aes(n_cs, n_coloc_cs, color = QTL, shape = QTL)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_smooth(method = "lm",
              alpha = 0.1,
              show.legend = FALSE) %>%
  ggblend::partition(vars(QTL)) %>%
  ggblend::blend("multiply") +
  geom_point() +
  locusviz::get_default_theme() +
  theme(
    legend.position.inside = c(0, 1),
    legend.justification.inside = c(0, 1),
    legend.spacing = unit(0, "mm")
  ) +
  labs(x = "# 95% CS", y = "# colocalized 95% CS", color = "QTL") +
  scale_color_manual(values = qtl.colors)
p.coloc.rate

df.coloc.rate.category =
  dplyr::group_by(df.coloc.rate, QTL, trait_category) %>%
  dplyr::summarize(
    n_total = sum(n_cs),
    n_coloc = sum(n_coloc_cs),
    locusviz::binom_ci(n_coloc, n_total),
    .groups = "drop"
  ) %>%
  dplyr::rename(frac_coloc = frac) %>%
  dplyr::mutate(QTL = factor(QTL, levels = c("pQTL", "eQTL", "caQTL")))

if (!file.exists("tables/ST27_coloc_rate_category.tsv")) {
  export_table(df.coloc.rate.category,
               "tables/ST27_coloc_rate_category.tsv",
               "ST27")
  
}

pd = position_dodge(width = 0.75)
p.coloc.rate.category =
  dplyr::filter(df.coloc.rate.category, QTL %in% c("caQTL", "eQTL")) %>%
  ggplot(aes(
    frac_coloc,
    trait_category,
    color = QTL,
    shape = QTL
  )) +
  geom_hline(
    yintercept = c(2.5, 3.5),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(position = pd) +
  geom_errorbarh(aes(xmin = frac_lower, xmax = frac_upper),
                 height = 0,
                 position = pd) +
  locusviz::get_default_theme(hide.ytitle = TRUE, legend.position = "none") +
  scale_x_continuous(labels = scales::label_percent(drop0trailing = TRUE)) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = qtl.colors) +
  scale_shape_manual(values = qtl.shapes) +
  labs(x = "% GWAS coloc") +
  coord_cartesian(xlim = c(0, 0.7))
p.coloc.rate.category

################################################################################
# SMR

dplyr::filter(df.coloc, QTL %in% c("caQTL", "eQTL")) %>%
  dplyr::distinct(trait1, region1) %>%
  nrow()
# 11584

dplyr::distinct(df.coloc, trait2, QTL) %>%
  dplyr::group_by(QTL) %>%
  dplyr::summarize(count = n(), .groups = "drop")

dplyr::group_by(df.coloc, trait1, region1) %>%
  dplyr::summarize(
    eQTL_coloc = any(QTL == "eQTL"),
    caQTL_coloc = any(QTL == "caQTL"),
    eQTL_smr = any(!is.na(p_SMR) & QTL == "eQTL"),
    caQTL_smr = any(!is.na(p_SMR) & QTL == "caQTL"),
    .groups = "drop"
  ) %>%
  dplyr::summarize(
    n_eQTL_coloc = sum(eQTL_coloc),
    n_caQTL_coloc = sum(caQTL_coloc),
    n_eQTL_smr = sum(eQTL_smr),
    n_caQTL_smr = sum(caQTL_smr)
  )
# n_eQTL_coloc n_caQTL_coloc n_eQTL_smr n_caQTL_smr
#        <int>         <int>      <int>       <int>
#         5974         11104       1698        4048

p.cascade.coloc.smr =
  dplyr::mutate(
    df.cascade.coloc,
    y = interaction(
      scales::comma(n_total),
      scales::comma(n_coloc),
      scales::comma(n_smr),
      qtl_mechanism_category
    )
  ) %>%
  dplyr::rename_with( ~ sub("^(frac_(coloc|smr))$", "\\1_value", .x)) %>%
  tidyr::pivot_longer(
    cols = tidyselect::starts_with("frac_"),
    names_to = c("type", ".value"),
    names_pattern = "frac_(coloc|smr)_(value|lower|upper)"
  ) %>%
  dplyr::rename(frac = value) %>%
  dplyr::mutate(type = ifelse(type == "smr", "SMR", "Coloc")) %>%
  tidyr::drop_na() %>%
  ggplot(aes(frac, y, color = type)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  locusviz::get_default_theme(
    hide.ytitle = TRUE,
    legend.position = c(0, 1),
    legend.justification = c(0, 1)
  ) +
  scale_x_continuous(labels = scales::label_percent(drop0trailing = TRUE)) +
  scale_color_manual(values = BuenColors::jdb_palette("corona")[c(3, 4)]) +
  guides(y = legendry::guide_axis_nested(
    key = legendry::key_range_auto(sep = "\\."),
    levels_text = list(
      element_text(color = "grey50"),
      element_text(color = "grey50"),
      element_text(color = "grey50"),
      element_text()
    )
  )) +
  labs(x = "% variants", color = "Analysis") +
  coord_cartesian(xlim = c(0, 0.13))

p.cascade.coloc.smr.rate =
  dplyr::mutate(
    df.cascade.coloc,
    y = interaction(
      scales::comma(n_total),
      scales::comma(n_coloc),
      scales::comma(n_smr),
      qtl_mechanism_category
    )
  ) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(locusviz::binom_ci(n_smr, n_coloc)) %>%
  ggplot(aes(frac, y)) +
  geom_point(color = "grey20") +
  geom_errorbarh(aes(xmin = frac_lower, xmax = frac_upper), height = 0) +
  locusviz::get_default_theme(hide.ylab = TRUE, legend.position = "none") +
  scale_x_continuous(labels = scales::label_percent(drop0trailing = TRUE)) +
  labs(x = "# SMR / # GWAS coloc") +
  coord_cartesian(xlim = c(0, 0.71))

df.n_smr =
  dplyr::filter(df.coloc, QTL %in% c("eQTL", "caQTL")) %>%
  dplyr::group_by(QTL, trait1) %>%
  dplyr::summarize(
    n_coloc = length(unique(region1)),
    n_smr = length(unique(region1[!is.na(p_SMR)])),
    .groups = "drop"
  )

p.smr.rate =
  ggplot(df.n_smr, aes(n_coloc, n_smr, color = QTL, shape = QTL)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_smooth(method = "lm",
              alpha = 0.1,
              show.legend = FALSE) %>%
  ggblend::partition(vars(QTL)) %>%
  ggblend::blend("multiply") +
  geom_point() +
  locusviz::get_default_theme() +
  theme(
    legend.position.inside = c(0, 1),
    legend.justification.inside = c(0, 1),
    legend.spacing = unit(0, "mm")
  ) +
  labs(x = "# colocalized loci", y = "# SMR loci", color = "QTL") +
  scale_color_manual(values = qtl.colors)
p.smr =
  p.smr.rate + p.cascade.coloc.smr + p.cascade.coloc.smr.rate +
  patchwork::plot_annotation(tag_levels = "a")
p.smr

cowplot::save_plot(
  "figures/SFig11_smr_rate.pdf",
  p.smr,
  base_height = 60,
  base_width = 180,
  units = "mm"
)

################################################################################
if (!file.exists("data/expressed_genes.rds")) {
  expressed_genes = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/saige_qtl/step3/integrated_gex_batch1_5.fgid.qc.*.mean.inv.SAIGE.acat.txt.gz", function(df, remote_path) {
    dplyr::mutate(df, cell_type = parse_cell_type(remote_path))
  }) %>%
    dplyr::distinct(phenotype_id) %>%
    dplyr::pull(phenotype_id)
  saveRDS(expressed_genes, "data/expressed_genes.rds")
}
expressed_genes = readRDS("data/expressed_genes.rds")

df.loeuf.v4 = rgsutil::read_gsfile(
  "gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
) %>%
  dplyr::inner_join(df.features, by = c("gene_id" = "phenotype_id")) %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile) %>%
  dplyr::filter(gene_id %in% expressed_genes)

if (!file.exists("data/df.gex.in_cs.max_abs_beta.tsv.gz")) {
  df.gex.in_cs = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/saige_qtl/susie/full/integrated_gex_batch1_5.fgid.qc.*.mean.inv.SAIGE.chr*.SUSIE.in_cs.snp.bgz", function(df, path) {
    if (nrow(df) == 0) {
      return(NULL)
    }
    dplyr::group_by(df, region, cs) %>%
      dplyr::mutate(abs_beta = abs(beta)) %>%
      dplyr::filter(!low_purity) %>%
      dplyr::filter(max(abs_beta) == abs_beta) %>%
      dplyr::ungroup()
  })
  data.table::fwrite(
    df.gex.in_cs,
    "data/df.gex.in_cs.max_abs_beta.tsv.gz",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    na = "NA"
  )
}

if (!file.exists("data/df.atac.in_cs.max_abs_beta.tsv.gz")) {
  df.atac.in_cs = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/atac_results/susie/full/integrated_atac_batch1_5.fgid.*.sum.inv.chr*.SUSIE.in_cs.snp.bgz", function(df, path) {
    if (nrow(df) == 0) {
      return(NULL)
    }
    dplyr::group_by(df, region, cs) %>%
      dplyr::mutate(abs_beta = abs(beta)) %>%
      dplyr::filter(!low_purity) %>%
      dplyr::filter(max(abs_beta) == abs_beta) %>%
      dplyr::ungroup()
  })
  data.table::fwrite(
    df.atac.in_cs,
    "data/df.atac.in_cs.max_abs_beta.tsv.gz",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    na = "NA"
  )
}

df.gex = readRDS("data/df.gex.rds") %>%
  dplyr::select(phenotype_id, cell_type, ACAT_q)

df.atac = readRDS("data/df.atac.rds") %>%
  dplyr::rename(ACAT_q = qval) %>%
  dplyr::select(phenotype_id, cell_type, ACAT_q)

df.gex.in_cs.all = data.table::fread("data/df.gex.in_cs.max_abs_beta.tsv.gz", data.table = FALSE) %>%
  dplyr::mutate(cell_type = parse_cell_type(trait)) %>%
  dplyr::left_join(df.gex, by = c("region" = "phenotype_id", "cell_type"))

df.gex.in_cs =
  dplyr::group_by(df.gex.in_cs.all, region, cs) %>%
  # dplyr::group_by(df.gex.in_cs.all, region) %>%
  dplyr::slice_max(abs_beta, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(abs_beta_norm = sqrt(2 * maf * (1 - maf)) * abs_beta)

df.atac.in_cs.all = data.table::fread("data/df.atac.in_cs.max_abs_beta.tsv.gz", data.table = FALSE) %>%
  dplyr::mutate(cell_type = parse_cell_type(trait)) %>%
  dplyr::left_join(df.atac, by = c("region" = "phenotype_id", "cell_type")) %>%
  dplyr::inner_join(
    dplyr::distinct(df.open4gene.sig, peak_id, gene_id, cell_type),
    by = c("region" = "peak_id", "cell_type"),
    relationship = "many-to-many"
  )

df.atac.in_cs =
  dplyr::group_by(df.atac.in_cs.all, region) %>%
  dplyr::slice_max(abs_beta, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    peak_id = region,
    region = gene_id,
    abs_beta_norm = sqrt(2 * maf * (1 - maf)) * abs_beta
  )

df.atac.combined.in_cs = data.table::fread("data/df.atac.in_cs.max_abs_beta.tsv.gz", data.table = FALSE) %>%
  dplyr::mutate(cell_type = parse_cell_type(trait)) %>%
  dplyr::inner_join(
    dplyr::filter(df.open4gene.sig, mode %in% c("dual", "rheostat")) %>%
      dplyr::group_by(peak_id, gene_id, cell_type) %>%
      dplyr::slice_max(abs(hurdle_count_beta), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup(),
    by = c("region" = "peak_id", "cell_type"),
    relationship = "many-to-many"
  ) %>%
  dplyr::mutate(abs_beta = abs(beta * hurdle_count_beta)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::slice_max(abs_beta, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(peak_id = region, region = gene_id)

df.link.abs_beta =
  dplyr::filter(df.open4gene.sig, mode %in% c("dual", "rheostat")) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarize(n_peaks = length(unique(peak_id)), abs_beta = max(abs(hurdle_count_beta))) %>%
  dplyr::rename(region = gene_id)

df.link.abs_beta.mode =
  dplyr::group_by(df.open4gene.sig, mode, gene_id) %>%
  dplyr::summarize(
    n_peaks = length(unique(peak_id)),
    abs_beta = max(abs(hurdle_count_beta)),
    distance = distance[which.max(abs(hurdle_count_beta))],
    .groups = "drop"
  ) %>%
  dplyr::rename(region = gene_id)

################################################################################
plot_loeuf_beta = function(df.abs_beta,
                           df.loeuf,
                           ylab = NULL,
                           ylim = NULL,
                           ybreaks = ggplot2::waiver(),
                           color.lab = NULL) {
  df = dplyr::left_join(df.loeuf, df.abs_beta, by = c("gene_id" = "region"))
  print(
    dplyr::filter(df, lof.oe_ci.upper_bin_decile %in% c(0, 9)) %>%
      dplyr::group_by(lof.oe_ci.upper_bin_decile) %>%
      dplyr::summarize(
        n_annot = length(unique(gene_id[!is.na(abs_beta)])),
        n_total = length(unique(gene_id)),
        locusviz::binom_ci(n_annot, n_total),
        locusviz::median_ci(abs_beta)
      )
  )
  ct = with(df, cor.test(lof.oe_ci.upper, abs_beta, method = "spearman"))
  print(ct)
  print(ct$p.value)
  
  dplyr::group_by(df, lof.oe_ci.upper_bin_decile) %>%
    dplyr::summarize(
      n_annot = length(unique(gene_id[!is.na(abs_beta)])),
      n_total = length(unique(gene_id)),
      locusviz::binom_ci(n_annot, n_total),
      locusviz::median_ci(abs_beta)
    ) %>%
    dplyr::mutate(x = (lof.oe_ci.upper_bin_decile + 0.5) / 10) %>%
    ggplot(aes(x, median, color = frac)) +
    geom_errorbar(aes(ymin = median_lower, ymax = median_upper), width = 0) +
    geom_point() +
    coord_cartesian(xlim = c(0, 1), ylim = ylim) +
    scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
    scale_y_continuous(expand = expansion(), breaks = ybreaks) +
    locusviz::get_default_theme(legend.position = c(0.5, 0),
                                legend.justification = c(0.5, 0)) +
    theme(
      legend.direction = "horizontal",
      legend.title = element_text(
        size = 6.4,
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      legend.title.position = "top",
      legend.key.width = unit(5, "mm"),
      legend.box.margin = margin(0, 0, 0, 8)
    ) +
    labs(x = "LOEUF decile", y = ylab, color = color.lab)
}

eqtl_color = scale_color_gradientn(
  colors = BuenColors::jdb_palette("brewer_red"),
  labels = scales::label_percent(),
  limit = c(0.6, 0.8)
)
caqtl_color = scale_color_gradientn(
  colors = BuenColors::jdb_palette("brewer_blue"),
  labels = scales::label_percent(),
  limit = c(0.3, 0.85)
)

# max per gene
p.constraint.gex =
  plot_loeuf_beta(
    df.gex.in_cs,
    df.loeuf.v4,
    ylab = expression(paste("Median |", italic(beta)[eQTL], "|")),
    ylim = c(0, 0.24),
    color.lab = "% genes with eQTL"
  ) +
  eqtl_color

length(unique(df.atac.in_cs$peak_id))
# 67790
length(unique(df.atac.in_cs$gene_id))
# 9044

# max per peak
p.constraint.atac =
  plot_loeuf_beta(
    df.atac.in_cs,
    df.loeuf.v4,
    ylab = expression(paste("Median |", italic(beta)[caQTL], "|")),
    ylim = c(0, 0.61),
    color.lab = "% genes with linked caPeaks"
  ) +
  caqtl_color

p.constraint.link =
  plot_loeuf_beta(df.link.abs_beta,
                  df.loeuf.v4,
                  ylab = expression(paste("Median |", italic(beta)[link], "|")),
                  ylim = c(0, 0.59)) +
  caqtl_color +
  theme(legend.position = "none")

p.constraint.link.n_peaks =
  dplyr::mutate(df.link.abs_beta, abs_beta = n_peaks) %>%
  plot_loeuf_beta(df.loeuf.v4, ylab = "Median # linked peaks", ylim = c(0, 60)) +
  caqtl_color +
  theme(legend.position = "none")

# max per gene
p.constraint.atac.combined =
  plot_loeuf_beta(
    df.atac.combined.in_cs,
    df.loeuf.v4,
    ylab = expression(paste("Median |", italic(beta)[combined], "|")),
    ylim = c(0, 0.28)
  ) +
  caqtl_color + theme(legend.position = "none")

df.eqtl.coloc =
  dplyr::filter(df.coloc, QTL == "eQTL") %>%
  dplyr::group_by(trait2) %>%
  dplyr::summarize(coloc = TRUE, .groups = "drop")

df.caqtl.coloc =
  dplyr::filter(df.coloc, QTL == "caQTL") %>%
  dplyr::inner_join(
    dplyr::distinct(df.open4gene.sig, peak_id, gene_id, cell_type),
    by = c("trait2" = "peak_id", "cell_type"),
    relationship = "many-to-many"
  ) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarize(coloc = TRUE, .groups = "drop") %>%
  dplyr::rename(trait2 = gene_id)

df.constraint.coloc =
  dplyr::bind_rows(
    dplyr::left_join(df.loeuf.v4, df.eqtl.coloc, by = c("gene_id" = "trait2")) %>%
      dplyr::mutate(QTL = "eQTL"),
    dplyr::left_join(df.loeuf.v4, df.caqtl.coloc, by = c("gene_id" = "trait2")) %>%
      dplyr::mutate(QTL = "caQTL")
  ) %>%
  dplyr::group_by(QTL, lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(
    n_coloc = length(unique(gene_id[coloc])),
    n_total = length(unique(gene_id)),
    locusviz::binom_ci(n_coloc, n_total),
    .groups = "drop"
  ) %>%
  dplyr::mutate(x = (lof.oe_ci.upper_bin_decile + 0.5) / 10)

p.constraint.coloc =
  ggplot(df.constraint.coloc, aes(x, frac, color = QTL)) +
  geom_errorbar(aes(ymin = frac_lower, ymax = frac_upper), width = 0) +
  geom_point() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.8)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_color_manual(values = qtl.colors) +
  locusviz::get_default_theme() +
  labs(x = "LOEUF decile", y = "% genes with GWAS coloc")


################################################################################
# Extended
p.constraint.gene.frac =
  dplyr::bind_rows(
    dplyr::left_join(df.loeuf.v4, df.gex.in_cs, by = c("gene_id" = "region")) %>%
      dplyr::mutate(QTL = "with eQTL"),
    dplyr::left_join(df.loeuf.v4, df.atac.in_cs, by = c("gene_id" = "region")) %>%
      dplyr::mutate(QTL = "with linked caPeak")
  ) %>%
  dplyr::group_by(QTL, lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(
    n_annot = length(unique(gene_id[!is.na(abs_beta)])),
    n_total = length(unique(gene_id)),
    locusviz::binom_ci(n_annot, n_total),
    .groups = "drop"
  ) %>%
  dplyr::mutate(x = (lof.oe_ci.upper_bin_decile + 0.5) / 10) %>%
  ggplot(aes(x, frac, color = QTL)) +
  geom_errorbar(aes(ymin = frac_lower, ymax = frac_upper), width = 0) +
  geom_point() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_color_manual(values = c(
    "with linked caPeak" = unname(qtl.colors["caQTL"]),
    "with eQTL" = unname(qtl.colors["eQTL"])
  )) +
  locusviz::get_default_theme() +
  theme(legend.title = element_blank()) +
  labs(x = "LOEUF decile", y = "% genes")

# link mechanism
df.constraint.link.mode =
  dplyr::left_join(df.loeuf.v4, df.link.abs_beta.mode, by = c("gene_id" = "region")) %>%
  dplyr::filter(mode %in% c("dual", "rheostat")) %>%
  dplyr::group_by(mode, lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(locusviz::median_ci(abs_beta), .groups = "drop")

dplyr::left_join(df.loeuf.v4, df.link.abs_beta.mode, by = c("gene_id" = "region")) %>%
  dplyr::filter(mode %in% c("dual", "rheostat")) %>%
  dplyr::group_by(mode) %>%
  dplyr::summarize(locusviz::spearman_ci(lof.oe_ci.upper, abs_beta), .groups = "drop")

p.constraint.link.mode =
  dplyr::mutate(
    df.constraint.link.mode,
    x = (lof.oe_ci.upper_bin_decile + 0.5) / 10,
    mode = factor(mode, levels = link_mode_levels)
  ) %>%
  ggplot(aes(x, median)) +
  # geom_errorbar(aes(ymin = median_lower, ymax = median_upper), width = 0) +
  geom_ribbon(aes(ymin = median_lower, ymax = median_upper, fill = mode), alpha = 0.1) +
  geom_line(aes(color = mode)) +
  geom_point(aes(color = mode)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.5)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  scale_color_manual(values = link_mode_colors, labels = stringr::str_to_title) +
  scale_fill_manual(values = link_mode_colors, guide = "none") +
  locusviz::get_default_theme(legend.position = c(0, 1),
                              legend.justification = c(0, 1)) +
  theme(legend.title = element_blank()) +
  labs(x = "LOEUF decile", y = expression(paste("Median |", italic(beta)[link], "|")))
p.constraint.link.mode

df.constraint.link.n_peaks.mode =
  dplyr::left_join(df.loeuf.v4, df.link.abs_beta.mode, by = c("gene_id" = "region")) %>%
  dplyr::filter(mode %in% c("dual", "rheostat", "switch")) %>%
  dplyr::group_by(mode, lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(locusviz::median_ci(n_peaks), .groups = "drop")

df.constraint.link.n_peaks.mode %>% dplyr::filter(lof.oe_ci.upper_bin_decile %in% c(0, 9))

p.constraint.link.n_peaks.mode =
  dplyr::mutate(
    df.constraint.link.n_peaks.mode,
    x = (lof.oe_ci.upper_bin_decile + 0.5) / 10,
    mode = factor(mode, levels = link_mode_levels)
  ) %>%
  ggplot(aes(x, median)) +
  # geom_errorbar(aes(ymin = median_lower, ymax = median_upper), width = 0) +
  geom_ribbon(aes(ymin = median_lower, ymax = median_upper, fill = mode), alpha = 0.1) +
  geom_line(aes(color = mode)) +
  geom_point(aes(color = mode)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 35)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  scale_color_manual(values = link_mode_colors, labels = stringr::str_to_title) +
  scale_fill_manual(values = link_mode_colors, guide = "none") +
  locusviz::get_default_theme() +
  theme(legend.title = element_blank()) +
  labs(x = "LOEUF decile", y = "Median # linked peaks")
p.constraint.link.n_peaks.mode


df.constraint.link.distance.mode =
  dplyr::left_join(df.loeuf.v4, df.link.abs_beta.mode, by = c("gene_id" = "region")) %>%
  dplyr::filter(mode %in% c("dual", "rheostat", "switch")) %>%
  dplyr::group_by(mode, lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(locusviz::median_ci(abs(distance)), .groups = "drop")


p.constraint.link.distance.mode =
  dplyr::mutate(
    df.constraint.link.distance.mode,
    x = (lof.oe_ci.upper_bin_decile + 0.5) / 10,
    mode = factor(mode, levels = link_mode_levels)
  ) %>%
  ggplot(aes(x, median)) +
  # geom_errorbar(aes(ymin = median_lower, ymax = median_upper), width = 0) +
  geom_ribbon(aes(ymin = median_lower, ymax = median_upper, fill = mode), alpha = 0.1) +
  geom_line(aes(color = mode)) +
  geom_point(aes(color = mode)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 55e4)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(
    expand = expansion(),
    labels = scales::label_number(scale = 1e-6, drop0trailing = TRUE)
  ) +
  scale_color_manual(values = link_mode_colors, labels = stringr::str_to_title) +
  scale_fill_manual(values = link_mode_colors, guide = "none") +
  locusviz::get_default_theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  labs(x = "LOEUF decile", y = "Median abs. distance to TSS (Mb)")
p.constraint.link.distance.mode

p.ext.res =
  list(
    p.constraint.link.n_peaks.mode + labs(tag = "h"),
    p.constraint.link.mode + labs(tag = "i"),
    p.constraint.link.distance.mode + labs(tag = "j")
  ) %>%
  purrr::reduce(`+`)

cowplot::save_plot(
  "local/figures/ExtendedDataFig9_constraint_loeuf_reponse_hij.pdf",
  p.ext.res,
  base_height = 60,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "local/figures/ExtendedDataFig9_constraint_loeuf_reponse_hij.png",
  p.ext.res,
  base_height = 60,
  base_width = 180,
  units = "mm",
  dpi = 300
)


p.constraint.gex.norm =
  dplyr::mutate(df.gex.in_cs, abs_beta = abs_beta_norm) %>%
  plot_loeuf_beta(
    df.loeuf.v4,
    ylab = expression(paste("Median normalized |", italic(beta)[eQTL], "|")),
    ylim = c(0, 0.1),
    ybreaks = seq(0, 0.1, by = 0.02),
    color.lab = "% genes with eQTL"
  ) +
  eqtl_color

p.constraint.atac.norm =
  dplyr::mutate(df.atac.in_cs, abs_beta = abs_beta_norm) %>%
  plot_loeuf_beta(
    df.loeuf.v4,
    ylab = expression(paste("Median normalized |", italic(beta)[caQTL], "|")),
    ylim = c(0, 0.22),
    color.lab = "% genes with linked caPeaks"
  ) +
  caqtl_color

################################################################################

df.gex.raw.mean = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/pseudobulk/integrated_gex_batch1_5.fgid.qc.predicted.celltype.l1.PBMC.mean.raw.bed.gz"
) %>%
  dplyr::transmute(
    gene_id = gene_id,
    mean_expression = rowMeans(dplyr::select(., tidyselect::starts_with("FG")), na.rm = TRUE),
    var_expression = matrixStats::rowSds(as.matrix(
      dplyr::select(., tidyselect::starts_with("FG"))
    ), na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    log_mean = log(mean_expression),
    log_var = log(var_expression),
    resi_var = residuals(loess(
      log_var ~ log_mean, data = ., span = 0.1
    ))
  )

# take top 20% highly expressed genes and match based on 10 x 10 deciles
set.seed(42)
df.gex.raw.mean.subset <-
  dplyr::inner_join(df.gex.raw.mean, df.loeuf.v4, by = "gene_id") %>%
  dplyr::filter(mean_expression >= quantile(mean_expression, 0.8)) %>%
  dplyr::mutate(expr_quantile = Hmisc::cut2(mean_expression, g = 10)) %>%
  dplyr::group_by(lof.oe_ci.upper_bin_decile, expr_quantile) %>%
  dplyr::mutate(n_available = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(expr_quantile) %>%
  dplyr::mutate(min_available = min(n_available)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(lof.oe_ci.upper_bin_decile, expr_quantile) %>%
  dplyr::sample_n(size = min(dplyr::first(min_available), n()), replace = FALSE) %>%
  dplyr::ungroup()

set.seed(42)
df.gex.raw.mean_var_double.subset <-
  dplyr::inner_join(df.gex.raw.mean, df.loeuf.v4, by = "gene_id") %>%
  dplyr::filter(mean_expression >= quantile(mean_expression, 0.8)) %>%
  dplyr::mutate(
    expr_quantile = Hmisc::cut2(mean_expression, g = 10),
    var_quantile = Hmisc::cut2(resi_var, g = 10),
  ) %>%
  dplyr::group_by(lof.oe_ci.upper_bin_decile, expr_quantile, var_quantile) %>%
  dplyr::mutate(n_available = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(expr_quantile, var_quantile) %>%
  dplyr::mutate(min_available = min(n_available)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(lof.oe_ci.upper_bin_decile, expr_quantile, var_quantile) %>%
  dplyr::sample_n(size = min(dplyr::first(min_available), n()), replace = FALSE) %>%
  dplyr::ungroup()


set.seed(42)
df.gex.raw.mean.subset.for_loeuf <-
  dplyr::inner_join(df.gex.raw.mean, df.loeuf.v4, by = "gene_id") %>%
  dplyr::mutate(
    expr_quantile = Hmisc::cut2(mean_expression, g = 10),
    eQTL = gene_id %in% df.gex.in_cs$region
  ) %>%
  dplyr::group_by(eQTL, expr_quantile) %>%
  dplyr::mutate(n_available = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(expr_quantile) %>%
  dplyr::mutate(min_available = min(n_available)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(eQTL, expr_quantile) %>%
  dplyr::sample_n(size = min(dplyr::first(min_available), n()), replace = FALSE) %>%
  dplyr::ungroup()

dplyr::group_by(df.gex.raw.mean.subset.for_loeuf, eQTL) %>%
  dplyr::summarize(n = n(), locusviz::median_ci(lof.oe_ci.upper))
# eQTL      n median median_lower median_upper
# <lgl> <int>  <dbl>        <dbl>        <dbl>
#   1 FALSE  3256  0.953        0.934        0.967
# 2 TRUE   3256  0.991        0.973        1.00

with(
  df.gex.raw.mean.subset,
  cor(mean_expression, lof.oe_ci.upper_bin_decile, method = "spearman")
)
# 0.001312212


dplyr::group_by(df.gex.raw.mean.subset, lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(
    n_genes = n(),
    median_expr = median(mean_expression),
    mean_expr = mean(mean_expression)
  )

df.constraint.gex.mean =
  dplyr::bind_rows(
    dplyr::inner_join(df.loeuf.v4, df.gex.raw.mean, by = "gene_id") %>%
      dplyr::group_by(lof.oe_ci.upper_bin_decile) %>%
      dplyr::summarize(exp_category = "All genes", locusviz::median_ci(mean_expression)),
    dplyr::inner_join(df.loeuf.v4, df.gex.raw.mean, by = "gene_id") %>%
      dplyr::filter(gene_id %in% df.gex.raw.mean.subset$gene_id) %>%
      dplyr::group_by(lof.oe_ci.upper_bin_decile) %>%
      dplyr::summarize(exp_category = "Matched (top 20%)", locusviz::median_ci(mean_expression))
  )

df.gex.raw.mean.loeuf = dplyr::inner_join(df.loeuf.v4, df.gex.raw.mean, by = "gene_id")
result =  with(
  df.gex.raw.mean.loeuf,
  cor.test(lof.oe_ci.upper, mean_expression, method = "spearman")
)
rho <- result$estimate
S <- result$statistic
n <- nrow(df.gex.raw.mean.loeuf)

z_stat <- Rmpfr::mpfr(rho, precBits = 200) * sqrt(Rmpfr::mpfr(n, precBits = 200) - 1)
p_value <- 2 * Rmpfr::pnorm(z_stat, lower.tail = TRUE)
print(p_value, digits = 50)

p.constraint.gex.mean =
  dplyr::mutate(df.constraint.gex.mean,
                x = (lof.oe_ci.upper_bin_decile + 0.5) / 10) %>%
  ggplot(aes(x, median, color = exp_category)) +
  geom_errorbar(aes(ymin = median_lower, ymax = median_upper), width = 0) +
  geom_point() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 3200)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  locusviz::get_default_theme(legend.position = c(1, 0.3),
                              legend.justification = c(1, 0.3)) +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c(
    "All genes" = unname(qtl.colors["eQTL"]),
    "Matched (top 20%)" = BuenColors::jdb_palette("brewer_purple")[7]
  )) +
  labs(x = "LOEUF decile", y = "Median mean expr. (CP10K)")
p.constraint.gex.mean

p.constraint.gex.var =
  dplyr::mutate(df.gex.raw.mean, abs_beta = resi_var, region = gene_id) %>%
  plot_loeuf_beta(df.loeuf.v4, ylim = c(-0.17, 0.02), ylab = "Median expr. variance residuals") +
  scale_color_gradientn(colors = unname(qtl.colors["eQTL"])) +
  theme(legend.position = "none")

df.gex.raw.mean.matched.subset =
  dplyr::inner_join(
    df.loeuf.v4,
    dplyr::bind_rows(
      dplyr::filter(df.gex.in_cs, region %in% df.gex.raw.mean.subset$gene_id) %>%
        dplyr::mutate(match = "Expr. matched (top 20%)"),
      dplyr::filter(
        df.gex.in_cs,
        region %in% df.gex.raw.mean_var_double.subset$gene_id
      ) %>%
        dplyr::mutate(match = "Expr. matched (top 20%)\n& variance matched")
    ),
    by = c("gene_id" = "region")
  ) %>%
  dplyr::group_by(match, lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(
    n_annot = length(unique(gene_id[!is.na(abs_beta)])),
    n_total = length(unique(gene_id)),
    locusviz::binom_ci(n_annot, n_total),
    locusviz::median_ci(abs_beta),
    .groups = "drop"
  )

p.constraint.gex.exp_matched =
  dplyr::mutate(df.gex.raw.mean.matched.subset,
                x = (lof.oe_ci.upper_bin_decile + 0.5) / 10) %>%
  ggplot(aes(x, median)) +
  geom_ribbon(aes(ymin = median_lower, ymax = median_upper, fill = match), alpha = 0.1) +
  geom_line(aes(color = match)) +
  geom_point(aes(color = match)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.2)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  locusviz::get_default_theme(legend.position = c(0, 1),
                              legend.justification = c(0, 1)) +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c(
    BuenColors::jdb_palette("brewer_purple")[7],
    BuenColors::jdb_palette("brewer_fire")[7]
  )) +
  scale_fill_manual(
    values = c(
      BuenColors::jdb_palette("brewer_purple")[7],
      BuenColors::jdb_palette("brewer_fire")[7]
    ),
    guide = "none"
  ) +
  labs(x = "LOEUF decile",
       y = expression(paste("Median |", italic(beta)[eQTL], "|")),
       color = "Subset")
p.constraint.gex.exp_matched

################################################################################

get_cell_type_level <- function(cell_type) {
  dplyr::case_when(
    cell_type == "predicted.celltype.l1.PBMC" ~ "L0",
    filter_l1_cell_types(cell_type) &
      cell_type != "predicted.celltype.l1.PBMC" ~ "L1",!filter_l1_cell_types(cell_type) ~ "L2"
  )
}

df.eqtl.coloc.per_level =
  dplyr::filter(df.coloc, QTL == "eQTL") %>%
  dplyr::mutate(cell_type_level = get_cell_type_level(cell_type)) %>%
  dplyr::group_by(trait2, cell_type_level) %>%
  dplyr::summarize(coloc = TRUE, .groups = "drop")

df.caqtl.coloc.per_level =
  dplyr::filter(df.coloc, QTL == "caQTL") %>%
  dplyr::mutate(cell_type_level = get_cell_type_level(cell_type)) %>%
  dplyr::inner_join(
    dplyr::distinct(df.open4gene.sig, peak_id, gene_id, cell_type),
    by = c("trait2" = "peak_id", "cell_type"),
    relationship = "many-to-many"
  ) %>%
  dplyr::group_by(gene_id, cell_type_level) %>%
  dplyr::summarize(coloc = TRUE, .groups = "drop") %>%  # <-- add this
  dplyr::rename(trait2 = gene_id)

df.constraint.coloc.per_level =
  dplyr::bind_rows(
    dplyr::left_join(
      tidyr::crossing(df.loeuf.v4, cell_type_level = c("L0", "L1", "L2")),
      df.eqtl.coloc.per_level,
      by = c("gene_id" = "trait2", "cell_type_level")
    ) %>% dplyr::mutate(QTL = "eQTL"),
    dplyr::left_join(
      tidyr::crossing(df.loeuf.v4, cell_type_level = c("L0", "L1", "L2")),
      df.caqtl.coloc.per_level,
      by = c("gene_id" = "trait2", "cell_type_level")
    ) %>% dplyr::mutate(QTL = "caQTL")
  ) %>%
  dplyr::group_by(QTL, cell_type_level, lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(
    n_coloc = length(unique(gene_id[coloc])),
    n_total = length(unique(gene_id)),
    locusviz::binom_ci(n_coloc, n_total),
    .groups = "drop"
  ) %>%
  dplyr::mutate(x = (lof.oe_ci.upper_bin_decile + 0.5) / 10)


eQTL.per_level.shades = BuenColors::jdb_palette("brewer_red")[c(8, 6, 4)]
caQTL.per_level.shades = BuenColors::jdb_palette("brewer_blue")[c(8, 6, 4)]

p.constraint.coloc.eQTL.per_level =
  dplyr::filter(df.constraint.coloc.per_level, QTL == "eQTL") %>%
  ggplot(aes(x, frac)) +
  geom_ribbon(aes(ymin = frac_lower, ymax = frac_upper, fill = cell_type_level),
              alpha = 0.1) +
  geom_line(aes(color = cell_type_level)) +
  geom_point(aes(color = cell_type_level)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.28)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_color_manual(values = eQTL.per_level.shades) +
  scale_fill_manual(values = eQTL.per_level.shades, guide = "none") +
  locusviz::get_default_theme() +
  labs(
    x = "LOEUF decile",
    y = "% genes with GWAS coloc",
    color = "Cell type level",
    title = "eQTL"
  )

p.constraint.coloc.caQTL.per_level =
  dplyr::filter(df.constraint.coloc.per_level, QTL == "caQTL") %>%
  ggplot(aes(x, frac)) +
  geom_ribbon(aes(ymin = frac_lower, ymax = frac_upper, fill = cell_type_level),
              alpha = 0.1) +
  geom_line(aes(color = cell_type_level)) +
  geom_point(aes(color = cell_type_level)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.75)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_color_manual(values = caQTL.per_level.shades) +
  scale_fill_manual(values = caQTL.per_level.shades, guide = "none") +
  locusviz::get_default_theme() +
  labs(
    x = "LOEUF decile",
    y = "% genes with GWAS coloc",
    color = "Cell type level",
    title = "caQTL"
  )

p.constraint.coloc.per_level =
  p.constraint.coloc.eQTL.per_level + p.constraint.coloc.caQTL.per_level +
  patchwork::plot_annotation(tag_levels = "a")
p.constraint.coloc.per_level

cowplot::save_plot(
  "figures/SFig18_coloc_per_level.pdf",
  p.constraint.coloc.per_level,
  base_height = 60,
  base_width = 120,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig18_coloc_per_level.png",
  p.constraint.coloc.per_level,
  base_height = 60,
  base_width = 120,
  units = "mm",
  dpi = 300
)

################################################################################

# RINT vs log2(CP10K)
df.gex.rint = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/gex_results/cis/integrated_gex_batch1_5.fgid.qc.predicted.celltype.l1.PBMC.mean.inv.cis_qtl_egenes.tsv.gz"
) %>%
  dplyr::mutate(
    region = phenotype_id,
    abs_beta = abs(slope),
    normalization = "RINT"
  ) %>%
  dplyr::filter(qval < 0.05)

df.gex.log2cp10k = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/log2cp10k_test/gex_results/cis/integrated_gex_batch1_5.fgid.qc.predicted.celltype.l1.PBMC.mean.log2cp10k.cis_qtl_egenes.tsv.gz"
) %>%
  dplyr::mutate(
    region = phenotype_id,
    abs_beta = abs(slope),
    normalization = "log2(CP10K)"
  ) %>%
  dplyr::filter(qval < 0.05)

df.contraint.gex.loeuf.method =
  dplyr::inner_join(
    df.loeuf.v4,
    dplyr::bind_rows(df.gex.rint, df.gex.log2cp10k),
    by = c("gene_id" = "region")
  )

dplyr::group_by(df.contraint.gex.loeuf.method, normalization) %>%
  dplyr::summarize(locusviz::spearman_ci(abs_beta, lof.oe_ci.upper))

p.contraint.gex.loeuf.method =
  dplyr::group_by(df.contraint.gex.loeuf.method,
                  normalization,
                  lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(locusviz::median_ci(abs_beta), .groups = "drop") %>%
  dplyr::mutate(
    x = (lof.oe_ci.upper_bin_decile + 0.5) / 10,
    normalization = factor(normalization, levels = c("RINT", "log2(CP10K)"))
  ) %>%
  ggplot(aes(x, median, color = normalization)) +
  geom_errorbar(aes(ymin = median_lower, ymax = median_upper), width = 0) +
  geom_point() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.35)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  scale_color_manual(
    values = c(
      "RINT" = BuenColors::jdb_palette("corona")[11],
      "log2(CP10K)" = BuenColors::jdb_palette("corona")[10]
    )
  ) +
  locusviz::get_default_theme(legend.position = c(1, 0),
                              legend.justification = c(1, 0)) +
  labs(
    x = "LOEUF decile",
    y = expression(paste("Median |", italic(beta)[eQTL], "|")),
    color = "Normalization",
    title = "PBMC (tensorQTL)"
  )
p.contraint.gex.loeuf.method

################################################################################
# gnocchi
df.gnocchi = rgsutil::read_gsfile(
  "gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/genomic_constraint/constraint_z_genome_1kb.qc.download.txt.gz"
) %>%
  dplyr::mutate(q = dplyr::percent_rank(z))

df.consensus = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/atac/peaks/integrated_atac_batch1_5.iterative_merged.bed.gz",
  header = FALSE
) %>%
  dplyr::transmute(
    chrom = V1,
    start = V2,
    end = V3,
    peak_id = paste(chrom, start + 1, end, sep = "-")
  )

gr.gnocchi = GenomicRanges::makeGRangesFromDataFrame(df.gnocchi, keep.extra.columns = TRUE)
gr.consensus = GenomicRanges::makeGRangesFromDataFrame(df.consensus, keep.extra.columns = TRUE)

gr.consensus.overlaps <- GenomicRanges::findOverlaps(gr.consensus, gr.gnocchi)

df.consensus.overlap <- tibble::tibble(
  consensus_idx = S4Vectors::queryHits(gr.consensus.overlaps),
  gnocchi_idx = S4Vectors::subjectHits(gr.consensus.overlaps),
  z = df.gnocchi$z[S4Vectors::subjectHits(gr.consensus.overlaps)],
  q = df.gnocchi$q[S4Vectors::subjectHits(gr.consensus.overlaps)]
) %>%
  dplyr::group_by(consensus_idx) %>%
  # dplyr::summarize(z = mean(z, na.rm = TRUE), .groups = 'drop')
  dplyr::summarize(
    z = max(z, na.rm = TRUE),
    q = max(q, na.rm = TRUE),
    .groups = 'drop'
  )

df.consensus.gnocchi <-
  dplyr::mutate(df.consensus, consensus_idx = dplyr::row_number()) %>%
  dplyr::left_join(df.consensus.overlap) %>%
  dplyr::select(-consensus_idx) %>%
  tidyr::drop_na()

df.gnocchi.combined =
  dplyr::inner_join(
    df.loeuf.v4,
    dplyr::bind_rows(
      dplyr::filter(df.open4gene.sig, peak_id %in% df.consensus.gnocchi$peak_id) %>%
        dplyr::inner_join(df.consensus.gnocchi, by = "peak_id") %>%
        dplyr::group_by(gene_id) %>%
        dplyr::slice_max(z, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::rename(region = gene_id) %>%
        dplyr::mutate(abs_beta = z, peak_group = "Most constrained"),
      dplyr::filter(
        df.open4gene.sig,
        peak_id %in% df.consensus.gnocchi$peak_id &
          mode %in% c("dual", "rheostat")
      ) %>%
        dplyr::mutate(abs_beta = abs(hurdle_count_beta)) %>%
        dplyr::inner_join(df.consensus.gnocchi, by = "peak_id") %>%
        dplyr::group_by(gene_id) %>%
        dplyr::slice_max(abs_beta, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::rename(region = gene_id) %>%
        dplyr::mutate(abs_beta = z, peak_group = "Strongest effect"),
      dplyr::filter(df.open4gene.sig, peak_id %in% df.consensus.gnocchi$peak_id) %>%
        dplyr::inner_join(df.consensus.gnocchi, by = "peak_id") %>%
        dplyr::group_by(gene_id) %>%
        dplyr::sample_n(1) %>%
        dplyr::ungroup() %>%
        dplyr::rename(region = gene_id) %>%
        dplyr::mutate(abs_beta = z, peak_group = "Random")
    ),
    by = c("gene_id" = "region")
  ) %>%
  dplyr::group_by(peak_group, lof.oe_ci.upper_bin_decile) %>%
  dplyr::summarize(
    n_annot = length(unique(gene_id[!is.na(abs_beta)])),
    n_total = length(unique(gene_id)),
    locusviz::binom_ci(n_annot, n_total),
    locusviz::median_ci(abs_beta)
  ) %>%
  dplyr::mutate(
    x = (lof.oe_ci.upper_bin_decile + 0.5) / 10,
    peak_group = factor(
      peak_group,
      levels = c("Most constrained", "Strongest effect", "Random")
    )
  )

p.gnocchi.combined =
  ggplot(df.gnocchi.combined, aes(x, median, color = peak_group)) +
  geom_errorbar(aes(ymin = median_lower, ymax = median_upper), width = 0) +
  geom_point() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 6.3)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  scale_color_manual(
    values = c(
      "Strongest effect" = BuenColors::jdb_palette("corona")[3],
      "Most constrained" = BuenColors::jdb_palette("corona")[4],
      "Random" = "grey50"
    )
  ) +
  locusviz::get_default_theme() +
  theme(legend.title = element_blank()) +
  labs(x = "LOEUF decile", y = "Gnocchi", color = "Peak per gene")
p.gnocchi.combined

df.z =
  dplyr::bind_rows(
    dplyr::mutate(df.gnocchi, category = "All regions"),
    dplyr::mutate(df.consensus.gnocchi, category = "Consensus peaks"),
    dplyr::filter(df.consensus.gnocchi, peak_id %in% df.open4gene.sig$peak_id) %>%
      dplyr::mutate(category = "Linked peaks")
  )

p.gnocchi.z =
  ggplot(df.z, aes(z, fill = category)) +
  geom_histogram(bins = 100) +
  geom_vline(
    aes(xintercept = median),
    data = dplyr::group_by(df.z, category) %>% dplyr::summarize(median = median(z)),
    linetype = "dashed"
  ) +
  locusviz::get_default_theme() +
  theme(legend.title = element_blank()) +
  scale_x_continuous(expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  scale_fill_manual(
    values = c(
      "All regions" = "grey80",
      "Consensus peaks" = unname(qtl.colors["caQTL"]),
      "Linked peaks" = BuenColors::jdb_palette("corona")[3]
    )
  ) +
  labs(x = "Gnocchi", y = "Count")

if (!file.exists("data/df.gnocchi.loeuf.rds")) {
  df.gnocchi.loeuf =
    dplyr::filter(
      df.open4gene.sig,
      peak_id %in% df.consensus.gnocchi$peak_id &
        mode %in% c("dual", "rheostat")
    ) %>%
    dplyr::mutate(abs_beta = abs(hurdle_count_beta)) %>%
    # dplyr::group_by(peak_id, gene_id) %>%
    # dplyr::slice_max(abs_beta, n = 1, with_ties = FALSE) %>%
    # dplyr::ungroup() %>%
    data.table::as.data.table() %>%
    (\(dt) dt[dt[, .I[which.max(abs_beta)], by = .(peak_id, gene_id)]$V1])() %>%
    tibble::as_tibble() %>%
    dplyr::inner_join(df.consensus.gnocchi, by = "peak_id") %>%
    dplyr::inner_join(df.loeuf.v4, by = "gene_id") %>%
    dplyr::mutate(z_bin = cut(z, c(-4, 4, 10)),
                  x = (lof.oe_ci.upper_bin_decile %/% 2 + 0.5) / 5) %>%
    tidyr::drop_na(z_bin) %>%
    dplyr::select(x, z_bin, abs_beta) %>%
    dplyr::group_by(x, z_bin) %>%
    dplyr::summarize(n = n(),
                     locusviz::median_ci(abs_beta),
                     .groups = "drop")
  saveRDS(df.gnocchi.loeuf, "data/df.gnocchi.loeuf.rds")
}
df.gnocchi.loeuf = readRDS("data/df.gnocchi.loeuf.rds")

p.gnocchi.loeuf =
  dplyr::filter(df.gnocchi.loeuf, z_bin %in% c("(-4,4]", "(4,10]")) %>%
  ggplot(aes(x, median, color = z_bin)) +
  geom_errorbar(aes(ymin = median_lower, ymax = median_upper), width = 0) +
  geom_point() +
  labs(x = "LOEUF quintile",
       y = expression(paste("Median |", italic(beta)[link], "|")),
       color = "Gnocchi") +
  locusviz::get_default_theme(legend.position = c(1, 0),
                              legend.justification = c(1, 0)) +
  scale_color_manual(values = BuenColors::jdb_palette("brewer_violet")[c(5, 8)]) +
  scale_x_continuous(
    labels = scales::label_percent(),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion()
  ) +
  scale_y_continuous(expand = expansion()) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.151))
p.gnocchi.loeuf

################################################################################
p.constraint.ext =
  list(
    p.constraint.gene.frac,
    p.constraint.gex.mean,
    p.constraint.gex.var,
    p.constraint.gex.exp_matched,
    p.contraint.gex.loeuf.method,
    p.constraint.gex.norm,
    p.constraint.atac.norm,
    p.constraint.link.n_peaks.mode,
    p.constraint.link.mode,
    p.constraint.link.distance.mode,
    p.gnocchi.combined,
    p.gnocchi.loeuf
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 4) +
  patchwork::plot_annotation(tag_levels = "a")

p.constraint.ext

cowplot::save_plot(
  "figures/ExtendedDataFig9_constraint_loeuf.pdf",
  p.constraint.ext,
  base_height = 170,
  base_width = 180,
  units = "mm"
)
cowplot::save_plot(
  "figures/ExtendedDataFig9_constraint_loeuf.png",
  p.constraint.ext,
  base_height = 170,
  base_width = 180,
  units = "mm",
  dpi = 300
)

################################################################################
p1 = p.constraint.gene.frac + p.constraint.gex + p.constraint.atac +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(tag_levels = "a")
p2 = (p.constraint.link.n_peaks + labs(tag = "d")) + (p.constraint.link + labs(tag = "e")) +
  patchwork::plot_layout(nrow = 1)
p3 = ((p.constraint.atac + labs(tag = "c")) + p.constraint.link + labs(tag = "e")) + (p.constraint.atac.combined + labs(tag = "f")) + (p.constraint.gex + labs(tag = "b")) +
  patchwork::plot_layout(nrow = 1)

cowplot::save_plot(
  "~/Downloads/p.constaint.slides.p1.pdf",
  p1,
  base_height = 2,
  base_width = 7.08661 / 4 * 3
)
cowplot::save_plot(
  "~/Downloads/p.constaint.slides.p2.pdf",
  p2,
  base_height = 2,
  base_width = 7.08661 / 4 * 2
)
cowplot::save_plot(
  "~/Downloads/p.constaint.slides.p3.pdf",
  p3,
  base_height = 2,
  base_width = 7.08661 / 4 * 4
)

################################################################################
df.nfee = rgsutil:::fread_wrapper("data/R12_annotated_variants_v1.small.filt.AF.gz") %>%
  dplyr::transmute(
    variant = `#variant`,
    AF = AF,
    GENOME_enrichment_nfee = GENOME_enrichment_nfee,
    AFE = dplyr::case_when(GENOME_enrichment_nfee > 5 ~ "AFE > 5", TRUE ~ "AFE ≤ 5")
  ) %>%
  tidyr::drop_na(GENOME_enrichment_nfee)

df.variant.nfee = dplyr::left_join(df.variant, df.nfee, by = c("variant_id" = "variant"))

if (!file.exists("tables/ST32_FIN_molQTL.tsv")) {
  dplyr::filter(df.variant.nfee, GENOME_enrichment_nfee > 5) %>%
    dplyr::select(
      variant_id,
      AF,
      GENOME_enrichment_nfee,
      qtl_mechanism_category,
      best_cell_types:cascade_peak_genes
    ) %>%
    export_table("tables/ST32_FIN_molQTL.tsv", "ST32", save_googlesheet = FALSE)
}

dplyr::filter(df.variant.nfee, GENOME_enrichment_nfee > 5) %>%
  nrow()
# 8124

dplyr::filter(df.variant.nfee, GENOME_enrichment_nfee > 5) %>%
  dplyr::group_by(qtl_mechanism_category) %>%
  dplyr::summarize(count = n())
# qtl_mechanism_category count
# <fct>                  <int>
#   1 Only eQTL                493
# 2 Only caQTL (No Link)    2892
# 3 Only caQTL (With Link)  4475
# 4 caQTL + eQTL (No Link)   101
# 5 Distal Cascade            54
# 6 Positional Cascade        49
# 7 Local Cascade             60

df.eqtl_catalogue = read.table("data/finemapped.AFE5.tsv",
                               T,
                               sep = "\t",
                               comment.char = "") %>%
  dplyr::filter(X.resource == "eQTL_Catalogue_R7" & pip > 0.5)

dplyr::distinct(df.eqtl_catalogue, variant_id) %>%
  nrow()
# 78

df.gex.max_pip.nfee = dplyr::left_join(df.gex.max_pip, df.nfee) %>%
  dplyr::mutate(QTL = "eQTL")
df.atac.max_pip.nfee = dplyr::left_join(df.atac.max_pip, df.nfee) %>%
  dplyr::mutate(QTL = "caQTL")

data =
  dplyr::bind_rows(df.gex.max_pip.nfee, df.atac.max_pip.nfee) %>%
  tidyr::drop_na(max_pip, GENOME_enrichment_nfee) %>%
  dplyr::filter(GENOME_enrichment_nfee > 0) %>%
  dplyr::mutate(
    MAF = 0.5 - abs(0.5 - AF),
    max_pip_bin = cut(max_pip, c(-Inf, 0.01, 0.5, 1)),
    GENOME_enrichment_nfee = ifelse(
      is.infinite(GENOME_enrichment_nfee),
      GENOME_enrichment_nfee[is.finite(GENOME_enrichment_nfee)],
      GENOME_enrichment_nfee
    )
  )

p.MAF =
  dplyr::filter(data, max_pip_bin == "(0.5,1]") %>%
  ggplot(aes(MAF, QTL)) +
  ggridges::geom_density_ridges(aes(fill = AFE), alpha = 0.8, linewidth = 0) +
  scale_x_log10(
    labels = scales::label_number(drop0trailing = TRUE),
    breaks = c(0.001, 0.01, 0.1, 0.5)
  ) +
  locusviz::get_default_theme(
    hide.ytitle = TRUE,
    legend.position = c(0.72, 1.1),
    legend.justification = c(0.72, 1.1)
  ) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(
    values = c(BuenColors::jdb_palette("corona")[3], "grey50"),
    labels = c("AFE > 5", expression(AFE <= 5))
  ) +
  scale_y_discrete(expand = expansion())

df.enrichment.AFE =
  dplyr::bind_rows(
    locusviz::compute_functional_enrichment(
      df.gex.max_pip.nfee,
      annot_levels = c("AFE > 5", "AFE ≤ 5"),
      pip_bin_breaks = c(-Inf, 0.01, 0.5, 1),
      consequence_col = "AFE"
    ) %>%
      dplyr::mutate(QTL = "eQTL"),
    locusviz::compute_functional_enrichment(
      df.atac.max_pip.nfee,
      annot_levels = c("AFE > 5", "AFE ≤ 5"),
      pip_bin_breaks = c(-Inf, 0.01, 0.5, 1),
      consequence_col = "AFE"
    ) %>%
      dplyr::mutate(QTL = "caQTL")
  ) %>%
  dplyr::mutate(QTL = factor(QTL, levels = c("caQTL", "eQTL")))

if (!file.exists("tables/ST33_AFE_enrichment.tsv")) {
  dplyr::rename(df.enrichment.AFE, AFE = consequence) %>%
    dplyr::select(QTL, tidyselect::everything()) %>%
    export_table("tables/ST33_AFE_enrichment.tsv", "ST33")
}

pd = position_dodge(width = 0.75)
p.enrichment.AFE =
  ggplot(df.enrichment.AFE, aes(enrichment, consequence, color = QTL)) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0,
    position = pd,
    show.legend = FALSE
  ) +
  geom_point(aes(shape = QTL), position = pd) +
  scale_x_continuous(expand = expansion(),
                     labels = scales::label_number(drop0trailing = TRUE)) +
  scale_y_discrete(limits = rev, labels = c(expression(AFE <= 5), "AFE > 5")) +
  scale_color_manual(values = qtl.colors) +
  scale_shape_manual(values = qtl.shapes) +
  locusviz::get_default_theme(
    hide.ytitle = TRUE,
    legend.position = c(1, 0),
    legend.justification = c(1, 0)
  ) +
  coord_cartesian(xlim = c(0.8, 3.2)) +
  labs(x = "Enrichment")

p.AFE =
  p.MAF + p.enrichment.AFE + patchwork::plot_layout(nrow = 1) + patchwork::plot_annotation(tag_levels = "a")
cowplot::save_plot(
  "figures/SFig22_AFE_MAF.pdf",
  p.AFE,
  base_height = 60,
  base_width = 120,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig22_AFE_MAF.png",
  p.AFE,
  base_height = 60,
  base_width = 120,
  units = "mm",
  dpi = 300
)


df.gex.in_cs.AFE =
  dplyr::mutate(df.gex.in_cs,
                variant = stringr::str_c("chr", stringr::str_replace_all(v, ":", "_"))) %>%
  dplyr::inner_join(df.nfee, by = "variant")

df.atac.in_cs.AFE =
  dplyr::mutate(df.atac.in_cs,
                variant = stringr::str_c("chr", stringr::str_replace_all(v, ":", "_"))) %>%
  dplyr::inner_join(df.nfee, by = "variant")

df.gex.in_cs.AFE.loeuf =
  dplyr::inner_join(df.gex.in_cs.AFE, df.loeuf.v4, by = c("region" = "gene_id")) %>%
  dplyr::mutate(
    AFE_bin = dplyr::case_when(
      is.na(GENOME_enrichment_nfee) ~ NA_character_,
      GENOME_enrichment_nfee <= 5 ~ "<= 5",
      TRUE ~ "> 5"
    ),
    AFE_bin = factor(AFE_bin, levels = c("> 5", "<= 5")),
    x = (lof.oe_ci.upper_bin_decile %/% 2 + 0.5) / 5
  ) %>%
  tidyr::drop_na(AFE_bin) %>%
  dplyr::select(AFE_bin, x, abs_beta) %>%
  dplyr::group_by(AFE_bin, x) %>%
  dplyr::summarize(locusviz::median_ci(abs_beta), .groups = "drop")

df.atac.in_cs.AFE.loeuf =
  dplyr::inner_join(df.atac.in_cs.AFE, df.loeuf.v4, by = c("region" = "gene_id")) %>%
  dplyr::mutate(
    AFE_bin = dplyr::case_when(
      is.na(GENOME_enrichment_nfee) ~ NA_character_,
      GENOME_enrichment_nfee <= 5 ~ "<= 5",
      TRUE ~ "> 5"
    ),
    AFE_bin = factor(AFE_bin, levels = c("> 5", "<= 5")),
    x = (lof.oe_ci.upper_bin_decile %/% 2 + 0.5) / 5
  ) %>%
  tidyr::drop_na(AFE_bin) %>%
  dplyr::select(AFE_bin, x, abs_beta) %>%
  dplyr::group_by(AFE_bin, x) %>%
  dplyr::summarize(locusviz::median_ci(abs_beta), .groups = "drop")

df.gex.in_cs.AFE.rho =
  dplyr::inner_join(df.gex.in_cs.AFE, df.loeuf.v4, by = c("region" = "gene_id")) %>%
  dplyr::mutate(
    AFE_bin = dplyr::case_when(GENOME_enrichment_nfee <= 5 ~ "<= 5", TRUE ~ "> 5"),
    AFE_bin = factor(AFE_bin, levels = c("> 5", "<= 5")),
    x = (lof.oe_ci.upper_bin_decile %/% 2 + 0.5) / 5
  ) %>%
  dplyr::group_by(AFE_bin) %>%
  dplyr::summarize(locusviz::spearman_ci(abs_beta, lof.oe_ci.upper), .groups = "drop")

p.constraint.gex.AFE =
  ggplot(df.gex.in_cs.AFE.loeuf, aes(x, median, color = AFE_bin)) +
  geom_errorbar(aes(ymin = median_lower, ymax = median_upper), width = 0) +
  geom_point() +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  scale_color_manual(values = c(unname(qtl.colors["eQTL"]), "grey50"),
                     labels = c("AFE > 5", expression(AFE <= 5))) +
  locusviz::get_default_theme(legend.position = c(0, 1),
                              legend.justification = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.34)) +
  labs(x = "LOEUF quintile",
       y = expression(paste("Median |", italic(beta)[eQTL], "|")),
       color = "Finnish AFE")

p.constraint.atac.AFE =
  ggplot(df.atac.in_cs.AFE.loeuf, aes(x, median, color = AFE_bin)) +
  geom_errorbar(aes(ymin = median_lower, ymax = median_upper), width = 0) +
  geom_point() +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  scale_color_manual(values = c(unname(qtl.colors["caQTL"]), "grey50"),
                     labels = c("AFE > 5", expression(AFE <= 5))) +
  locusviz::get_default_theme(legend.position = c(0, 1),
                              legend.justification = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.2)) +
  labs(x = "LOEUF quintile",
       y = expression(paste("Median |", italic(beta)[caQTL], "|")),
       color = "Finnish AFE")

layout =
  c(
    patchwork::area(1, 1, 1, 3),
    patchwork::area(1, 4, 1, 9),
    # coloc bar
    patchwork::area(2, 1, 2, 12),
    # d
    patchwork::area(1, 10, 1, 12),
    # constraints
    patchwork::area(3, 1, 3, 3),
    patchwork::area(4, 1, 4, 3),
    patchwork::area(3, 4, 3, 6),
    patchwork::area(3, 7, 3, 9),
    patchwork::area(4, 4, 4, 6),
    patchwork::area(4, 7, 4, 9),
    # FIN enrichment
    patchwork::area(3, 10, 3, 12),
    patchwork::area(4, 10, 4, 12)
  )
plot(layout)

plt.coloc =
  list(
    A = p.coloc.rate,
    B = patchwork::free(p.cascade.coloc, side = "l"),
    C = patchwork::free(p.coloc.bar + theme(plot.tag = element_blank()), side = "l"),
    D = patchwork::free(p.coloc.rate.category, side = "l"),
    E = p.constraint.gex,
    F = p.constraint.coloc,
    G = p.constraint.atac,
    H = p.constraint.link,
    I = p.constraint.link.n_peaks,
    J = p.constraint.atac.combined,
    K = p.constraint.gex.AFE,
    L = p.constraint.atac.AFE
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(design = layout, heights = c(1.4, 1.1, 1, 1)) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/Fig5_coloc_rate.pdf",
  plt.coloc,
  base_height = 170,
  base_width = 180,
  units = "mm"
)
