# Tissue-of-origin aging clocks reveal composite aging states across cancers

Analysis code for Bock *et al.*

> [Citation to be added upon publication]

Tissue-specific transcriptomic aging clocks are trained on GTEx normal tissue using Hallmark ssGSEA pathway scores, then applied to matched TCGA tumors. Each tumor's age-matched deviation is decomposed into per-pathway contributions with [ClockSHAP](https://github.com/montilab/ClockSHAP), distributed separately as an R package.

## Contents

**`part1_pancancer_aging/`** — age associations across six assay layers and 33 TCGA projects (Figure 1)

- `acquisition/` — download and SummarizedExperiment construction, one pair of scripts per layer
- `metadata/` — clinical harmonization and per-cancer covariate selection
- `templates/` — per-cancer age-association analysis, with its renderer, per layer
- `Rmd/` — aggregation, Hallmark enrichment, tissue-of-origin clustering
- `figures/`, `html/`

**`part2_gtex_tcga_comparison/`** — age trends compared between GTEx normal tissue and TCGA tumors (Figure 2)

- `Rmd/` — processing, pathway age meta-analysis, cell-type deconvolution
- `figures/`, `html/`

**`part3_clocks/`** — clocks, deviation, and ClockSHAP (Figures 3 to 6, all supplements)

- `Rmd/` — templates parameterized by cancer and tissue pair
- `render/` — runs each template across pairs; holds the pairings, hyperparameters, and seeds
- `figures/` — main figure and supplement scripts
- `html/` — rendered reports, including a full worked example for BRCA

Rendered `.html` reports are tracked. They carry the code, diagnostics, and numeric output of each step, and are the readable record for anyone who cannot run the pipeline.

## Notes

Part 1 runs download and metadata harmonization first, then per-layer processing, covariate selection, and analysis. Part 3 templates run in the order `01`, `02`, `04`, `03`, `07`, since step `03` reads objects written by step `04`; numbering reflects development order and is kept so filenames match the rendered reports.

Figure panels are written as individual PDFs and assembled in Adobe Illustrator.

## Data

No data are included. Sources: GTEx and TCGA expression via recount3, GTEx donor ages via dbGaP under controlled access, TCGA molecular data from the NCI Genomic Data Commons, copy number from Broad Firehose GISTIC, clinical and survival data from the TCGA Pan-Cancer Clinical Data Resource, SCAN-B from SCAN-B from Mendeley Data (Staaf et al. 2022), METABRIC from cBioPortal, Hallmark gene sets from MSigDB.

This pipeline cannot be run from a clone. It requires controlled-access donor ages, several hundred gigabytes of expression data, and a compute cluster. Paths are set through `params` in each script.

`renv.lock` pins package versions for Part 3. Parts 1 and 2 were run under R 4.3.1 and 4.4.3.

## License

GPL-3.0. See [LICENSE](LICENSE).
