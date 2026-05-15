library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R/const.R"))

df.loeuf.v4 = rgsutil::read_gsfile(
  "gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
) %>%
  dplyr::inner_join(df.features, by = c("gene_id" = "phenotype_id")) %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile)

df.downsample.conditions =
  rgsutil::read_gsfile("gs://expansion_areas/multiome/batch1_5/downsample/conditions/conditions.tsv")

if (!file.exists("data/df.gex.downsample.rds")) {
  df.gex.downsample = rgsutil::map_dfr_gsfiles("gs://expansion_areas/multiome/batch1_5/downsample/gex_results/cis/*.cis_qtl_acat.tsv.gz", function(df, path) {
    dplyr::mutate(
      df,
      label = stringr::str_replace(path, "^.*(d[0-9]+_c[0-9]+_rep[0-9]+).*$", "\\1"),
      condition_id = stringr::str_remove(label, "_rep[0-9]+"),
      rep = as.numeric(stringr::str_extract(label, "_rep([0-9]+)", group = 1)),
      normalization = stringr::str_extract(path, "mean\\.(inv|log2cp10k)\\.", group = 1)
    )
  })
  saveRDS(df.gex.downsample, "data/df.gex.downsample.rds")
}

df.gex.downsample = readRDS("data/df.gex.downsample.rds")
df.gex.downsample.loeuf =
  dplyr::left_join(df.gex.downsample,
                   df.loeuf.v4,
                   by = c("phenotype_id" = "gene_id")) %>%
  tidyr::drop_na(lof.oe_ci.upper_bin_decile)

meta_analyze_rho = ~ {
  d <- dplyr::filter(.x, !is.na(rho), !is.na(rho_lower), !is.na(rho_upper))
  
  # No valid reps
  if (nrow(d) == 0) {
    return(tibble::tibble(
      rho = NA_real_,
      rho_lower = NA_real_,
      rho_upper = NA_real_,
      k = 0L
    ))
  }
  
  # Only one valid rep — fall back to that rep's own CI
  if (nrow(d) == 1) {
    return(tibble::tibble(
      rho = d$rho,
      rho_lower = d$rho_lower,
      rho_upper = d$rho_upper,
      k = 1L
    ))
  }
  
  # Fisher z with floor on zero variance (at rho = ±1, CI width is 0)
  z <- atanh(pmin(pmax(d$rho, -0.9999), 0.9999))
  z_var <- ((atanh(pmin(
    pmax(d$rho_upper, -0.9999), 0.9999
  )) -
    atanh(pmin(
      pmax(d$rho_lower, -0.9999), 0.9999
    ))) / (2 * 1.96))^2
  z_var <- pmax(z_var, 1e-4)  # floor to avoid degenerate studies
  
  fit <- tryCatch(
    metafor::rma(yi = z, vi = z_var, method = "REML"),
    error = function(e)
      NULL
  )
  
  if (is.null(fit)) {
    # Fall back to simple Fisher z + t-CI across reps
    n_rep <- nrow(d)
    z_mean <- mean(z)
    z_se <- sd(z) / sqrt(n_rep)
    t_crit <- qt(0.975, df = n_rep - 1)
    tibble::tibble(
      rho = tanh(z_mean),
      rho_lower = tanh(z_mean - t_crit * z_se),
      rho_upper = tanh(z_mean + t_crit * z_se),
      k = n_rep
    )
  } else {
    tibble::tibble(
      rho = tanh(fit$b[1, 1]),
      rho_lower = tanh(fit$ci.lb),
      rho_upper = tanh(fit$ci.ub),
      k = fit$k
    )
  }
}

df.gex.downsample.frac.rho =
  dplyr::group_by(
    df.gex.downsample.loeuf,
    label,
    condition_id,
    rep,
    normalization,
    lof.oe_ci.upper_bin_decile
  ) %>%
  dplyr::summarize(locusviz::binom_ci(sum(qval < 0.05), n()), .groups = "drop") %>%
  dplyr::group_by(label, condition_id, rep, normalization) %>%
  dplyr::summarize(locusviz::spearman_ci(frac, lof.oe_ci.upper_bin_decile),
                   .groups = "drop") %>%
  dplyr::left_join(df.downsample.conditions, by = "condition_id")

df.gex.downsample.frac.rho.meta <-
  df.gex.downsample.frac.rho %>%
  dplyr::group_by(condition_id, normalization) %>%
  dplyr::group_modify(meta_analyze_rho) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(df.downsample.conditions, by = "condition_id") %>%
  dplyr::mutate(Normalization = ifelse(normalization == "inv", "RINT", "log2(CP10K)"))

df.gex.downsample.rho =
  dplyr::filter(df.gex.downsample.loeuf, qval < 0.05) %>%
  dplyr::group_by(label, condition_id, rep, normalization) %>%
  dplyr::summarize(locusviz::spearman_ci(abs(slope), lof.oe_ci.upper_bin_decile),
                   .groups = "drop") %>%
  dplyr::mutate(w = 1 / (rho_upper - rho_lower)^2) %>%
  dplyr::left_join(df.downsample.conditions, by = "condition_id")

df.gex.downsample.rho.meta <-
  df.gex.downsample.rho %>%
  dplyr::group_by(condition_id, normalization) %>%
  dplyr::group_modify(meta_analyze_rho) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(df.downsample.conditions, by = "condition_id") %>%
  dplyr::mutate(Normalization = ifelse(normalization == "inv", "RINT", "log2(CP10K)"))

plot_downsample_rho = function(df, xvar, ylim, legend.position = c(1, 0.9)) {
  ggplot(df, aes(.data[[xvar]], rho)) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "grey50") +
    geom_ribbon(aes(
      ymin = rho_lower,
      ymax = rho_upper,
      fill = Normalization
    ),
    alpha = 0.1) +
    geom_point(aes(color = Normalization)) +
    geom_line(aes(color = Normalization)) +
    scale_y_continuous(expand = expansion(c(0, 0.1))) +
    coord_cartesian(xlim = c(0, 1000), ylim = ylim) +
    locusviz::get_default_theme(legend.position = legend.position, legend.justification = legend.position) +
    theme(plot.title = element_text(hjust = 0.02)) +
    scale_color_manual(values = BuenColors::jdb_palette("corona")) +
    scale_fill_manual(values = BuenColors::jdb_palette("corona"))
}

p.downsample =
  list(
    dplyr::filter(df.gex.downsample.frac.rho.meta, cells_per_donor == 1000) %>%
      plot_downsample_rho("n_donors", c(-1.2, 0.5)) +
      labs(
        x = "# donors",
        y = expression(paste("% eGene-LOEUF Spearman's ", italic(rho))),
        title = "# cells per donor = 1000"
      ),
    dplyr::filter(df.gex.downsample.rho.meta, cells_per_donor == 1000) %>%
      plot_downsample_rho("n_donors", c(0, 0.35), legend.position = "none") +
      labs(x = "# donors", y = expression(
        paste("|", italic(beta)[eQTL], "|-LOEUF Spearman's ", italic(rho))
      )),
    dplyr::filter(df.gex.downsample.frac.rho.meta, n_donors == 1000) %>%
      plot_downsample_rho("cells_per_donor", c(-1.2, 0.5), legend.position = "none") +
      labs(
        x = "# cells per donor",
        y = expression(paste("% eGene-LOEUF Spearman's ", italic(rho))),
        title = "# donors = 1000"
      ),
    dplyr::filter(df.gex.downsample.rho.meta, n_donors == 1000) %>%
      plot_downsample_rho("cells_per_donor", c(0, 0.35), legend.position = "none") +
      labs(x = "# cells per donor", y = expression(
        paste("|", italic(beta)[eQTL], "|-LOEUF Spearman's ", italic(rho))
      ))
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig19_buffering_downsample.pdf",
  p.downsample,
  base_height = 120,
  base_width = 120,
  units = "mm"
)
cowplot::save_plot(
  "figures/SFig19_buffering_downsample.png",
  p.downsample,
  base_height = 120,
  base_width = 120,
  units = "mm",
  dpi = 300
)
