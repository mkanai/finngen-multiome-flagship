# Population-scale immune multiome atlas reveals regulatory disease mechanisms

Code to reproduce the figures and tables for:

Kanai, M. et al. [Population-scale multiome immune cell atlas reveals complex disease drivers](https://doi.org/10.1101/2025.11.25.25340489). medRxiv (2025)

## Repository layout

- `R/` — figure-generation scripts. See [`R/README.md`](R/README.md) for the script-to-figure index.
- `figures/` — generated figures (.pdf).

## Related repositories

- **CASCADE framework**: <https://github.com/mkanai/cascade>
- **fasthurdle** (hurdle model for peak-gene links): <https://github.com/mkanai/fasthurdle>
- **ldcov** (covariate-adjusted LD): <https://github.com/mkanai/ldcov>
- **Single-cell processing pipeline**: <https://github.com/FINNGEN/singlecell-pipeline>
- **FinnGen GWAS pipeline**: <https://github.com/FINNGEN/regenie-pipelines>
- **FinnGen fine-mapping pipeline**: <https://github.com/FINNGEN/finemapping-pipeline>
- **FinnGen colocalization pipeline**: <https://github.com/FINNGEN/coloc.susie.direct>

## Data availability

FinnGen as a research project is granted use of national healthcare data and biospecimens according to national and European regulations (GDPR), which preclude the research project from distributing individual level data. The FinnGen single-nucleus caQTL and eQTL summary statistics are currently subject to embargo according to the FinnGen consortium agreement. They are available upon request and being prepared for public release by Q4 2026 or at the time of publication. Individual-level data, including single-cell count matrices, are available to approved researchers through FinnGen’s data access procedures (<https://www.finngen.fi/en/access/results>). Academic users wishing to work with the FinnGen project resource directly can follow the procedures described at: <https://www.finngen.fi/en/how-we-collaborate>. The interactive CASCADE browser will be publicly available upon public release.

The following resources are currently publicly available:

- **FinnGen GWAS summary statistics (R12)**: <https://r12.finngen.fi/>
- **KANTA lab value GWAS summary statistics**: <https://labvalues.finngen.fi/>
- **FinnGen (R12) + MVP + UKBB meta-analysis results**: <https://mvp-ukbb.finngen.fi/>

## License

MIT License

## Contact

Masahiro Kanai (<mkanai@broadinstitute.org>)
