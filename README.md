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
- [x] **Phase 2, config tables:** `config/` tables replace hard-coded row
      positions and species-name special cases; `R/config.R` checks them
      against the rdats
- [ ] **Phase 3, functions and bug fixes:** split the code into `R/` functions;
      bring back the EwE fleet grouping and stanza start lengths
- [ ] **Phase 4, verify:** regression tests against `tests/baseline/` and
      landings/discards reconciliation tests
- [ ] **Phase 5, document:** `renv` lockfile and a methods write-up

## Layout

```
config/              species, fleet and index tables, plus global settings
data/bam_rdats/      BAM .rdat output files (see data/bam_rdats/README.md)
R/                   pipeline functions (config.R: load and validate config)
legacy/              original script, lightly edited to produce the baseline
tests/baseline/      baseline outputs used for regression testing
tests/testthat/      automated tests
```

## Configuration (`config/`)

All species-specific choices live here, not in the code. Species are matched
by `species_id` and file name, never by the species name stored inside the
rdat.

**`settings.yml`:** data folder, Ecopath base year (1995), model area
(336,192 km²), and which index CV to report (`weighted` = `cv.U.*`, the CV
used in the assessment fit).

**`species.csv`:** one row per rdat.

| Column | Meaning |
|---|---|
| `species_id` | Short code; matches the file prefix and the index names |
| `file` | rdat file in `data/bam_rdats/` |
| `include` | `TRUE` to process the species |
| `ewe_stanza_ages` | Start age of each EwE stanza, e.g. `"0,1,4"`; `0` = single stanza |
| `maturity_col` | `a.series` column with female maturity |
| `gutted_to_whole` | `none`, `by_age` (whole/gutted ratio at each age), or `mean_ratio_gt1` (mean of ratios > 1; gag) |
| `altM_target` | M the Lorenzen curve is scaled to: `M.constant` or `M.oldest_age` |
| `lw_length_unit`, `lw_weight_unit` | Units of the `wgt.a`/`wgt.b` length-weight parameters |

**`fleets.csv`:** one row per BAM fleet, with its SEDAR number, `catch_type`
(`landings`, `discards`, or `both` when the same fleet code has both), sector,
and the EwE fleet it rolls up to. In the next Ecopath model each species is
its own fleet, so `ewe_fleet` is the `species_id`.

**`indices.csv`:** every abundance index and the `sel.age` selectivity used
to summarize it.

Check the config against the rdats:

```r
source("R/config.R")
validate_config(load_config())   # zero rows = everything matches
```

## Tests

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

## Open questions

- **Missing selectivities:** GAJ headboat and video indices, and the VSN CVID
  index, have no matching selectivity in the rdat (legacy output was NA). Is
  `sel.m.GR` right for GAJ HB, and `sel.m.CVT` for VSN CVID?
- **Fleet and index descriptions** in `config/` were inferred from BAM codes.
  Check the ones flagged "confirm" (e.g. VSN `cHTR`).

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
