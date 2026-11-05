# StockAssessments

This repo processes age-structured stock assessment outputs for South Atlantic
reef fish into inputs for the SAR Ecopath with Ecosim (EwE) model. Inputs
include biomass, F, M and Z by EwE stanza, landings and dead discards by
fleet, relative abundance indices with selectivity summaries, and
multi-stanza parameters.

## Status

The code is being moved here from the OneDrive `SAR EwE UF/BAM2EwE` folder
and refactored into a reproducible pipeline.

- [x] **Phase 1, migrate and capture a baseline:** rdats copied, repo set up,
      legacy script run end-to-end, baseline outputs saved
- [ ] **Phase 2, config tables:** `config/species.csv` and `config/fleets.csv`
      replace hard-coded row positions and species-name special cases
- [ ] **Phase 3, functions and bug fixes:** split the code into `R/` functions;
      bring back the EwE fleet grouping and stanza start lengths
- [ ] **Phase 4, verify:** regression tests against `tests/baseline/` and
      landings/discards reconciliation tests
- [ ] **Phase 5, document:** `renv` lockfile and a methods write-up

## Layout

```
data/bam_rdats/      BAM .rdat output files (see data/bam_rdats/README.md)
legacy/              original script, lightly edited to produce the baseline
tests/baseline/      baseline outputs used for regression testing
```

## Species

Black sea bass, gag, greater amberjack, gray triggerfish, golden tilefish,
red grouper, red porgy, red snapper, snowy grouper and vermilion snapper.
Blueline tilefish and scamp rdats are in `data/` but not processed yet.
Mutton snapper will be added from the new Stock Synthesis assessment.

## Running the legacy baseline

From the repo root (needs the `xlsx` package, which needs Java, plus `plyr`):

```bash
Rscript legacy/BAM2EwE_age0_corr.R
```
