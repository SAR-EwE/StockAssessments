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
- [x] **Phase 3, functions and bug fixes:** `R/` functions driven by
      `config/`, run by `run_bam2ewe.R`; legacy bugs fixed; stanza start
      lengths included
- [ ] **Phase 4, verify:** regression tests against `tests/baseline/` and
      landings/discards reconciliation tests
- [ ] **Phase 5, document:** `renv` lockfile and a methods write-up

## Layout

```
run_bam2ewe.R        runs the whole pipeline
config/              species, fleet and index tables, plus global settings
data/bam_rdats/      BAM .rdat output files (see data/bam_rdats/README.md)
R/                   pipeline functions (see below)
output/              generated EwE inputs (CSV tables + per-species xlsx)
legacy/              original script, lightly edited to produce the baseline
tests/baseline/      baseline outputs used for regression testing
tests/testthat/      automated tests
```

## Running the pipeline

From the repo root (needs the `yaml` and `openxlsx` packages):

```bash
Rscript run_bam2ewe.R
```

The script checks the config against the rdats, processes every species
with `include = TRUE`, confirms the landings, discards and biomass add back
up to the assessment totals, and writes to `output/`:

| File | Contents |
|---|---|
| `species.csv` | Species processed: file, SEDAR, years, stanzas, units |
| `stanza_timeseries.csv` | Year x stanza: biomass (mt), F, M, Z, landings and discards (mt) |
| `catch_by_fleet.csv` | Landings and discards (mt, whole weight) by BAM fleet, EwE fleet and stanza |
| `ecopath_inputs.csv` | Base-year values by stanza and in total, with t/km² and biomass accumulation |
| `indices.csv`, `indices_wide.csv` | Index observations (long, and year x index) |
| `index_summary.csv` | Mean CV and selectivity summary (A50, A95, Amax) for each index |
| `multistanza_params.csv` | vbK, Wmat, Winf, Wmat/Winf, proportion mature by stanza |
| `stanza_lengths.csv` | Length (mm) at the start of each stanza |
| `alternative_M.csv` | Lorenzen-based M by stanza (unscaled and scaled to assessment M) |
| `checks.csv` | Year-by-year reconciliation with the assessment totals |
| `xlsx/<species_id>.xlsx` | All of the above for one species, one sheet per table |

## Methods

- **Stanzas:** each assessment age goes to the last EwE stanza starting at or
  before it. Assessments that start at age 1 have no age-0 stanza.
- **Biomass:** mid-year numbers at age (`N.age.mdyr`) x weight at age
  (`a.series$wgt.mt`), summed within stanzas.
- **Catch:** `CLD.est.mats` landings (`Lw.*`) and dead discards (`Dw.*`) by
  fleet, converted from gutted to whole weight where needed and from 1000 lb
  to mt.
- **F** = (landings + dead discards) / mid-year biomass, by stanza.
- **M:** `a.series$M`, averaged within stanzas and weighted by mid-year
  numbers. **Z** = F + M.
- **Alternative M:** Lorenzen M = 3W(g)^-0.288 at monthly ages from von
  Bertalanffy growth (t0 = -0.5 when the assessment t0 < -1). The scaled
  version has mean M over the fully selected ages (mean F >= half of F at the
  oldest age) equal to the assessment M.
- **Index CVs:** mean of `cv.U.<index>` over years with data.

## R functions (`R/`)

| File | Purpose |
|---|---|
| `config.R` | Load and validate `config/` |
| `read_bam.R` | Read an rdat, model years, unit constants |
| `stanzas.R` | Age-to-stanza mapping, biomass and M by stanza |
| `catch.R` | Landings and discards by fleet and stanza |
| `indices.R` | Index series, CVs and selectivity summaries |
| `life_history.R` | Multi-stanza parameters, stanza start lengths |
| `natural_mortality.R` | Lorenzen alternative M |
| `ecopath_inputs.R` | Base-year Ecopath inputs |
| `checks.R` | Reconciliation with assessment totals |
| `process_species.R` | All steps for one species |
| `output.R` | Combine species and write CSV/xlsx |

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

## Differences from the legacy outputs

Biomass, F, M, Z, landings, discards, multi-stanza parameters, stanza lengths
and the index series match the legacy baseline to machine precision. These
changes are deliberate bug fixes:

- **Index CVs:** each index now gets its own weighted CV (`cv.U.<index>`).
  The legacy file repeated rows and attached CVs to the wrong indices for 8
  species.
- **Ecopath total row:** `BA_1yr` and `BA_3yr` for the total now use all
  stanzas (the legacy code used only stanza 1). Column `BA.1985` is renamed
  `BA_1yr`, and stanzas are labeled `age1`, `age4`, ... instead of 1, 2, ....
  Single-stanza gray triggerfish now has an `age0` row as well as `total`.
- **Alternative M:**
  - t0 follows the rule above; the legacy code always used -0.5.
  - The length-weight units come from config, fixing weights that were off
    by up to 70%.
  - Fully selected ages are taken from the F-at-age column names. Legacy used
    column positions, which shifted black sea bass, the one assessment with an
    age-0 column, by a year.
  - Scaling uses the mean over fully selected ages.

  Together these change the scaled stanza M by 2–31%, mostly in the age-0
  stanza.
- **Stanza lengths:** gray triggerfish has one age-0 row (legacy wrote a
  second age-0 row at 258 mm).

## Running the legacy baseline

From the repo root (needs the `xlsx` package, which needs Java, plus `plyr`):

```bash
Rscript legacy/BAM2EwE_age0_corr.R
```
