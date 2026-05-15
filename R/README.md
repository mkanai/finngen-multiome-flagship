# Figure generation scripts

Each script writes its outputs to [`../figures/`](../figures/). `const.R` provides shared constants and helpers (sourced by the figure scripts).

For figures assembled in Illustrator from multiple R-generated panels, only the combined PDF is linked.

## Main figures

- [`Fig1_study_overview.R`](Fig1_study_overview.R):
  - [Fig1_study_overview.pdf](../figures/Fig1_study_overview.pdf)
  - [SFig1_qc.pdf](../figures/SFig1_qc.pdf)
  - [SFig2_open4gene_l2.pdf](../figures/SFig2_open4gene_l2.pdf)
  - [SFig4_open4gene_mode.pdf](../figures/SFig4_open4gene_mode.pdf)
  - [SFig5_open4gene_encode_all.pdf](../figures/SFig5_open4gene_encode_all.pdf)
- [`Fig2_QTL_overview.R`](Fig2_QTL_overview.R):
  - [Fig2.overview.pdf](../figures/Fig2.overview.pdf)
  - [ExtendedDataFig3_qtl.pdf](../figures/ExtendedDataFig3_qtl.pdf)
  - [ExtendedDataFig4_qtl2.pdf](../figures/ExtendedDataFig4_qtl2.pdf)
  - [SFig7_trans_qtl.pdf](../figures/SFig7_trans_qtl.pdf)
- [`Fig3_celltype_specificity.R`](Fig3_celltype_specificity.R): [Fig3_celltype_specificity.pdf](../figures/Fig3_celltype_specificity.pdf)
- [`Fig4_cascade.R`](Fig4_cascade.R):
  - [Fig4_cascade.pdf](../figures/Fig4_cascade.pdf)
  - [ExtendedDataFig5_cascade.pdf](../figures/ExtendedDataFig5_cascade.pdf)
- [`Fig5_gwas_coloc.R`](Fig5_gwas_coloc.R):
  - [Fig5_gwas_coloc.pdf](../figures/Fig5_gwas_coloc.pdf)
  - [ExtendedDataFig9_constraint_loeuf.pdf](../figures/ExtendedDataFig9_constraint_loeuf.pdf)
  - [SFig11_smr.pdf](../figures/SFig11_smr.pdf)
  - [SFig18_coloc_per_level.pdf](../figures/SFig18_coloc_per_level.pdf)
  - [SFig22_AFE_MAF.pdf](../figures/SFig22_AFE_MAF.pdf)
- [`Fig6_vignettes.R`](Fig6_vignettes.R): [Fig6_vignettes.pdf](../figures/Fig6_vignettes.pdf)

## Extended Data figures

- [`ExDataFig2_co_accessibility.R`](ExDataFig2_co_accessibility.R): [ExtendedDataFig2_coacc.pdf](../figures/ExtendedDataFig2_coacc.pdf)
- [`ExDataFig6_MPRA_replication.R`](ExDataFig6_MPRA_replication.R):
  - [ExtendedDataFig6_mpra.pdf](../figures/ExtendedDataFig6_mpra.pdf)
  - [ExtendedDataFig10_mpra_alasoo_loeuf.pdf](../figures/ExtendedDataFig10_mpra_alasoo_loeuf.pdf)
- [`ExDataFig7_TICAM1_RNASET2_RHOH.R`](ExDataFig7_TICAM1_RNASET2_RHOH.R):
  - [ExtendedDataFig7_TICAM1_RNASET2_RHOH.pdf](../figures/ExtendedDataFig7_TICAM1_RNASET2_RHOH.pdf)
  - [SFig12_TICAM1_forest.pdf](../figures/SFig12_TICAM1_forest.pdf)
- [`ExDataFig8_PILRB_APOBEC3A.R`](ExDataFig8_PILRB_APOBEC3A.R):
  - [ExtendedDataFig8_PILRB_APOBEC3A.pdf](../figures/ExtendedDataFig8_PILRB_APOBEC3A.pdf)
  - [SFig14_PILRB_forest.pdf](../figures/SFig14_PILRB_forest.pdf)
  - [SFig15_APOBEC3A_3B.pdf](../figures/SFig15_APOBEC3A_3B.pdf)

## Supplementary figures

- [`SFig3_open4gene_schematic.R`](SFig3_open4gene_schematic.R): [SFig3_open4gene_schematic.pdf](../figures/SFig3_open4gene_schematic.pdf)
- [`SFig6_replication.R`](SFig6_replication.R):
  - [SFig6_peak_gene_replication.pdf](../figures/SFig6_peak_gene_replication.pdf)
  - [SFig8_TenK10K_OneK1K_replication.pdf](../figures/SFig8_TenK10K_OneK1K_replication.pdf)
  - [SFig9_lead_replication.pdf](../figures/SFig9_lead_replication.pdf)
- [`SFig10_stimulus_responsive_QTL.R`](SFig10_stimulus_responsive_QTL.R): [SFig10_stimulus_responsive_QTL.pdf](../figures/SFig10_stimulus_responsive_QTL.pdf)
- [`SFig13_coloc_replication.R`](SFig13_coloc_replication.R): [SFig13_coloc_replication.pdf](../figures/SFig13_coloc_replication.pdf)
- [`SFig16_ldsc_mesc.R`](SFig16_ldsc_mesc.R):
  - [SFig16_mesc_upset.pdf](../figures/SFig16_mesc_upset.pdf)
  - [SFig17_ldsc.pdf](../figures/SFig17_ldsc.pdf)
- [`SFig19_buffering_downsample.R`](SFig19_buffering_downsample.R): [SFig19_buffering_downsample.pdf](../figures/SFig19_buffering_downsample.pdf)
- [`SFig20_buffering_replication.R`](SFig20_buffering_replication.R): [SFig20_buffering_replication.pdf](../figures/SFig20_buffering_replication.pdf)
- [`SFig21_GRAMD1B_FAS.R`](SFig21_GRAMD1B_FAS.R): [SFig21_GRAMD1B_FAS.pdf](../figures/SFig21_GRAMD1B_FAS.pdf)
- [`SFig23_batch_LISI.R`](SFig23_batch_LISI.R):
  - [SFig23_batch_variance.pdf](../figures/SFig23_batch_variance.pdf)
  - [SFig24_batch_lisi.pdf](../figures/SFig24_batch_lisi.pdf)
- [`SFig25_mashr_sensitivity.R`](SFig25_mashr_sensitivity.R): [SFig25_mashr_sensitivity.pdf](../figures/SFig25_mashr_sensitivity.pdf)

## Supplementary tables

- [`ST34_n_peers.R`](ST34_n_peers.R): [ST34_n_peers.tsv](../tables/ST34_n_peers.tsv)
