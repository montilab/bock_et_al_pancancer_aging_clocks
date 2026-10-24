# Tissue-of-origin aging clocks reveal composite aging states across cancers

Analysis code for Bock *et al.*

> [Citation to be added upon publication]

Tissue-specific transcriptomic aging clocks are trained on GTEx normal tissue using Hallmark ssGSEA pathway scores, then applied to matched TCGA tumors. Each tumor's age-matched deviation is decomposed into per-pathway contributions with [ClockSHAP](https://github.com/montilab/ClockSHAP), distributed separately as an R package.

## Contents

**`part1_pancancer_aging/`** — age associations across six assay layers and 33 TCGA projects

- `acquisition/` — download and SummarizedExperiment construction, one pair of scripts per layer
- `metadata/` — clinical harmonization and per-cancer covariate selection
- `templates/` — per-cancer age-association analysis, with its renderer, per layer
- `Rmd/` — aggregation, Hallmark enrichment, tissue-of-origin clustering
- `figures/`, `html/`

**`part2_gtex_tcga_comparison/`** — age trends compared between GTEx normal tissue and TCGA tumors

- `Rmd/` — processing, pathway age meta-analysis, cell-type deconvolution
- `figures/`, `html/`

**`part3_clocks/`** — clocks, deviation, and ClockSHAP

- `Rmd/` — pipeline templates parameterized by cancer–tissue pair, in execution order: `01` processing, `02` clock training, `03` ClockSHAP, `04` survival and deviation. `06_scanb_processing_BRCA.Rmd` scores SCAN-B and writes the inputs for the breast figure scripts (its `save_*` params default to `false`)
- `render/` — runs each template across the eight pairs; holds the pairings, hyperparameters, and seeds
- `figures/` — main figure scripts, numbered by manuscript figure
- `supplement/` — supplementary figure and table scripts
- `html/` — rendered reports, including the full pipeline for BRCA

Rendered `.html` reports are tracked. They carry the code, diagnostics, and numeric output of each step, and are the readable record for anyone who cannot run the pipeline.

## Manuscript outputs

| Output | Script |
|---|---|
| Figure 1 | `part1_pancancer_aging/figures/fig1_pancancer_aging_landscape.Rmd` |
| Figure 2 | `part2_gtex_tcga_comparison/figures/fig2_GTEx-TCGA_AgeTrend-Comparison.Rmd` |
| Figure 3 | `part3_clocks/figures/fig3_clockshap_intro_GTEx_LUNG.Rmd` |
| Figure 4, Supplementary Figure 3 | `part3_clocks/figures/fig4_cross_cancer_deviation.Rmd` |
| Figure 5 | `part3_clocks/figures/fig5_deviation_BRCA_SCANB.Rmd` |
| Figure 6 | `part3_clocks/figures/fig6_survival_OS.Rmd` |
| Supplementary Figure 1 | `part3_clocks/supplement/supp_predicted_vs_chronological_age.Rmd` |
| Supplementary Figure 2, Supplementary Table 4 | `part3_clocks/supplement/supp_purity_deviation.Rmd` |
| Supplementary Figure 4 | `part3_clocks/supplement/supp_METABRIC_validation.Rmd` |
| Supplementary Figure 5 | `part3_clocks/supplement/supp_BRCA_clock_stability.Rmd` |
| Supplementary Figure 6 | `part3_clocks/supplement/supp_scBRCA_pseudobulk_deviation.Rmd` |
| Supplementary Table 1 | `part3_clocks/supplement/supp_GTEx_clock_nestedCV_performance.Rmd` |
| Supplementary Table 2 | `part3_clocks/supplement/supp_deviation_summary_by_cohort.Rmd` |
| Supplementary Table 3 | `part3_clocks/supplement/supp_OS_survival_PH_diagnostics.Rmd` |
| Supplementary Table 5 | `part3_clocks/supplement/supp_GTEx_clock_coefficients.Rmd` |

Scripts may produce panels that do not appear in the manuscript; these are marked in the script. Figure panels are written as individual PDFs and assembled in Adobe Illustrator.

## Notes

Part 1 runs download and metadata harmonization first, then per-layer processing, covariate selection, and analysis. Part 3 templates and renderers are numbered in execution order; the BRCA reports in `html/` were rendered before the final renumbering, so step numbers mentioned inside them may differ from the filenames.

## Data

No data are included. Sources: GTEx and TCGA expression via recount3, GTEx donor ages via dbGaP under controlled access, TCGA molecular data from the NCI Genomic Data Commons, copy number from Broad Firehose GISTIC, clinical and survival data from the TCGA Pan-Cancer Clinical Data Resource, tumor purity from the PanCanAtlas ABSOLUTE calls, SCAN-B from Mendeley Data (Staaf et al. 2022), METABRIC from cBioPortal, the single-cell breast cancer atlas from the published studies named in its script, Hallmark gene sets from MSigDB.

This pipeline cannot be run from a clone. It requires controlled-access donor ages, several hundred gigabytes of expression data, and a compute cluster. Paths are set through `params` in each script.

`renv.lock` pins package versions for Part 3. Parts 1 and 2 were run under R 4.3.1 and 4.4.3.

## License

GPL-3.0. See [LICENSE](LICENSE).
