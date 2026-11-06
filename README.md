# StockAssessments

This repo converts age-structured stock assessment outputs for South
Atlantic reef fish into inputs for the SAR Ecopath with Ecosim (EwE) model.

It reads Beaufort Assessment Model (BAM) output files (`.rdat`), groups
assessment ages into EwE stanzas, and writes:

- biomass, F, M and Z by stanza and year;
- landings and dead discards by fleet and stanza (one EwE fleet per species);
- base-year Ecopath inputs (biomass, F, M, Z, catch, biomass accumulation);
- relative abundance indices with CVs and selectivity summaries;
- multi-stanza parameters (vbK, Wmat/Winf, proportion mature, stanza start
  lengths);
- an alternative, Lorenzen-based M by stanza.

Species currently processed: black sea bass, gag, greater amberjack, gray
triggerfish, golden tilefish, red grouper, red porgy, red snapper, snowy
grouper and vermilion snapper. Blueline tilefish and scamp rdats are in the
repo but switched off. Mutton snapper will be added when its new Stock
Synthesis assessment is available (see [Other assessment models](#other-assessment-models)).

---

## Quick start

1. **Install R (4.4 or later) and the packages:**

   ```r
   install.packages(c("yaml", "openxlsx", "testthat"))
   ```

   To use the exact package versions the outputs were made with, run
   `renv::restore()` from the repo root. `renv.lock` records them.

2. **Run the pipeline** from the repo root (or open `StockAssessments.Rproj`
   in RStudio and source the script):

   ```bash
   Rscript run_bam2ewe.R
   ```

3. **Find the results in `output/`** (see [Outputs](#outputs)).

4. **Optionally, run the tests:**

   ```bash
   Rscript -e 'testthat::test_dir("tests/testthat")'
   ```

## What `run_bam2ewe.R` does

1. **Loads the config** from `config/` and **checks it against every rdat**
   with `validate_config()`. It stops with a list of problems if any fleet,
   index, column or unit doesn't match.
2. **Processes each species with `include = TRUE`** (`process_species()`):
   reads the rdat, assigns ages to stanzas, and calculates everything listed
   above.
3. **Checks that the results add back up** to the assessment's own totals:
   landings, discards and biomass for every year. It warns if any differ by
   more than 0.01%.
4. **Writes** the combined CSV tables and one Excel workbook per species to
   `output/`.

A full run takes about a minute.

---

## Updating to a new assessment

When a new BAM assessment, update or FES revision comes out for a species:

1. **Add the rdat.** Copy it into `data/bam_rdats/`. Use the file naming
   convention `<ID>_SEDAR<number>[_suffix].rdat`, e.g. `RSN_SEDAR99.rdat`.
   Delete the superseded file; git history keeps it.

2. **Draft the config rows** from the new file. In R, from the repo root:

   ```r
   for (f in list.files("R", full.names = TRUE)) source(f)
   sugg <- suggest_config("data/bam_rdats/RSN_SEDAR99.rdat", species_id = "RSN")
   sugg   # review the species, fleets and indices tables
   ```

   `suggest_config()` reads everything it can from the rdat: fleet codes and
   catch types, SEDAR number, indices and matching selectivities, maturity
   column, gutted vs. whole weight, and length-weight units. For an existing
   species it carries over the stanza ages, include flag, alternative-M
   target, and fleet and index descriptions. It prints a message about
   anything it couldn't fill in.

3. **Write the rows** into `config/`. This replaces only that species' rows:

   ```r
   write_suggested_config(sugg)
   ```

   Then open the three CSVs and check that species' rows. In particular:

   - **`ewe_stanza_ages`:** confirm the stanza breaks still suit the EwE
     model.
   - **`gutted_to_whole`:** if landings are in gutted weight, choose `by_age`
     or `mean_ratio_gt1`.
   - **New fleets:** add `description` text, and check `sector`, which is
     guessed from the code (codes starting with `c` or `COM` are
     commercial).
   - **Indices with a blank `selectivity`:** fill these in from the
     assessment report, or leave blank to report no selectivity summary.

4. **Validate:**

   ```r
   validate_config(load_config())   # zero rows = everything matches
   ```

   Fix any errors it lists. Warnings about indices without a selectivity are
   allowed.

5. **Run the pipeline and the tests.** The regression tests compare a
   species with `tests/baseline/` only while it uses the rdat the baseline
   was made from. An updated species is skipped automatically. The other
   checks (reconciliation, F = (L+D)/B, and so on) still run for it.

6. **Update `data/bam_rdats/README.md`** (the file table and citation), then
   commit the new rdat, the `config/` changes and the regenerated `output/`
   together.

### Adding a new species

Follow the same steps with a new `species_id`. `suggest_config()` sets
`include = FALSE` and leaves `ewe_stanza_ages` blank. Fill those in and set
`include` to `TRUE`.

### Other assessment models

The pipeline reads **BAM** output only. An assessment from another model,
such as Stock Synthesis (the new mutton snapper assessment) or ASAP, needs
its own reader that returns the same structure as `read_bam()`: numbers,
weight and M at age, catch at age by fleet, index series and selectivities,
growth and length-weight parameters. Everything downstream of the reader can
then be reused.

---

## Repository layout

```
run_bam2ewe.R        runs the whole pipeline
config/              species, fleet and index tables, plus global settings
data/bam_rdats/      BAM .rdat output files (see data/bam_rdats/README.md)
R/                   pipeline functions
output/              generated EwE inputs (CSV tables + per-species xlsx)
tests/testthat/      automated tests
tests/baseline/      outputs of the original script, used for regression tests
legacy/              original script, lightly edited to produce the baseline
renv.lock            package versions used
```

### R functions (`R/`)

| File | Purpose |
|---|---|
| `config.R` | Load `config/` and validate it against the rdats |
| `new_assessment.R` | `suggest_config()` / `write_suggested_config()` for new rdats |
| `read_bam.R` | Read an rdat; model years; unit constants |
| `stanzas.R` | Age-to-stanza mapping, biomass and M by stanza |
| `catch.R` | Landings and discards by fleet and stanza |
| `indices.R` | Index series, CVs and selectivity summaries |
| `life_history.R` | Multi-stanza parameters, stanza start lengths |
| `natural_mortality.R` | Lorenzen alternative M |
| `ecopath_inputs.R` | Base-year Ecopath inputs |
| `checks.R` | Reconciliation with the assessment totals |
| `process_species.R` | All steps for one species |
| `output.R` | Combine species and write CSV/xlsx |

---

## Configuration (`config/`)

All species-specific choices live here, not in the code. Species are matched
by `species_id` and file name, never by the species name stored inside the
rdat (e.g. the scamp rdat says "SA Red Grouper").

**`settings.yml`**

| Setting | Meaning |
|---|---|
| `data_dir` | Folder with the rdats |
| `base_year` | Ecopath base year (1995) |
| `model_area_km2` | Model area used for t/km² (336,192) |
| `index_cv` | `weighted` = `cv.U.*` (the CV used in the assessment fit), or `unweighted` |

**`species.csv`:** one row per rdat.

| Column | Meaning |
|---|---|
| `species_id` | Short code; matches the file prefix, index names and EwE fleet |
| `file` | rdat file in `data/bam_rdats/` |
| `common_name` | Name used in outputs |
| `include` | `TRUE` to process the species |
| `ewe_stanza_ages` | Start age of each EwE stanza, e.g. `"0,1,4"`; `0` alone = single stanza |
| `maturity_col` | `a.series` column with female maturity |
| `gutted_to_whole` | `none`; `by_age` (whole/gutted ratio at each age); or `mean_ratio_gt1` (mean of ratios > 1; gag) |
| `altM_target` | M the Lorenzen curve is scaled to: `M.constant` or `M.oldest_age` |
| `lw_length_unit`, `lw_weight_unit` | Units of the `wgt.a`/`wgt.b` length-weight parameters (checked against `a.series$weight`) |
| `notes` | Free text |

**`fleets.csv`:** one row per BAM fleet.

| Column | Meaning |
|---|---|
| `sedar` | SEDAR number; must match the rdat file name |
| `bam_fleet` | Fleet code in `CLD.est.mats` (`Lw.<code>`, `Dw.<code>`) |
| `catch_type` | `landings`, `discards`, or `both` |
| `sector` | `commercial` or `recreational` |
| `ewe_fleet` | EwE fleet; one fleet per species, so the `species_id` |
| `description` | Free text |

**`indices.csv`:** one row per abundance index.

| Column | Meaning |
|---|---|
| `bam_index` | Index code in `t.series` (`U.<code>.ob`) |
| `selectivity` | `sel.age` element summarized for the index (blank = none) |

---

## Outputs

`output/` (all weights in metric tons, whole weight):

| File | Contents |
|---|---|
| `species.csv` | Species processed: file, SEDAR, years, stanzas, assessment units |
| `stanza_timeseries.csv` | Year x stanza: biomass, F, M, Z, landings, discards |
| `catch_by_fleet.csv` | Landings and discards by BAM fleet, EwE fleet and stanza |
| `ecopath_inputs.csv` | Base-year values by stanza and in total, with t/km² and biomass accumulation (`BA_1yr`, `BA_3yr`) |
| `indices.csv`, `indices_wide.csv` | Index observations (long, and year x index) |
| `index_summary.csv` | Mean CV and selectivity summary (A50, A95, Amax) for each index |
| `multistanza_params.csv` | vbK, Wmat, Winf, Wmat/Winf, proportion mature by stanza |
| `stanza_lengths.csv` | Length (mm) at the start of each stanza |
| `alternative_M.csv` | Lorenzen-based M by stanza, unscaled and scaled to the assessment M |
| `checks.csv` | Year-by-year reconciliation with the assessment totals |
| `xlsx/<species_id>.xlsx` | All of the above for one species, one sheet per table |

## Methods

- **Stanzas:** each assessment age goes to the last EwE stanza that starts
  at or before it. An assessment that starts at age 1 has no age-0 stanza in
  the time series.
- **Biomass:** mid-year numbers at age (`N.age.mdyr`) x weight at age
  (`a.series$wgt.mt`), summed within stanzas.
- **Catch:** `CLD.est.mats` landings (`Lw.*`) and dead discards (`Dw.*`) by
  fleet and age, converted from gutted to whole weight where needed and from
  1000 lb to mt (x 0.453592), then summed within stanzas.
- **F** = (landings + dead discards) / mid-year biomass, by stanza.
- **M:** `a.series$M`, averaged within stanzas and weighted by mid-year
  numbers at age. **Z** = F + M.
- **Ecopath inputs:** values in `base_year`, plus `BA_1yr` = B(y+1)/B(y) − 1
  and `BA_3yr` = ΣB(y..y+2) / ΣB(y−3..y−1) − 1. The total row sums all
  stanzas.
- **Index summaries:**
  - CV is the mean of `cv.U.<index>` over years with data.
  - Selectivity (averaged over years if it varies by year) is summarized by
    the value closest to 1 and its age (Amax), and, on the ascending limb,
    the values closest to 0.5 and 0.95 and their ages.
- **Multi-stanza parameters:**
  - Wmat is the weight (lb) at the first age with ≥ 50% female maturity,
    and Winf the largest weight of mature ages.
  - Proportion mature is the mean maturity of the ages in each stanza.
  - Each stanza's start length is the smallest mean length at age in that
    stanza (0 for age 0).
- **Alternative M (Lorenzen 1996):**
  - M = 3W^−0.288 (W in g) at monthly ages. Weight comes from von
    Bertalanffy growth, with t0 = −0.5 when the assessment t0 < −1, and the
    length-weight relationship.
  - The scaled version is rescaled so mean M over the fully selected ages
    (mean F ≥ half of F at the oldest age) equals the assessment M.
  - Each stanza reports the survivorship-weighted mean, the M at its start
    age, and the M at its mid-age.

## Differences from the original OneDrive script

Biomass, F, M, Z, landings, discards, multi-stanza parameters, stanza lengths
and the index series match the original `BAM2EwE age0 corr.R` to machine
precision (the tests check this). These differences are deliberate bug
fixes:

- **Index CVs:** each index now gets its own weighted CV. The original
  repeated rows and attached CVs to the wrong indices for 8 of 10 species.
- **Ecopath total row:** `BA_1yr` and `BA_3yr` now use all stanzas; the
  original used only stanza 1. Gray triggerfish also gets an `age0` row.
- **Alternative M:** these fixes change the scaled stanza M by 2–31%,
  mostly in the age-0 stanza.
  - The t0 rule above is applied as intended; the original always used
    −0.5.
  - The length-weight units are correct; some were off by up to 70%.
  - Fully selected ages are taken by age, not by column position.
  - Scaling uses the mean over the fully selected ages.
- **Stanza lengths:** the duplicate gray triggerfish age-0 row is gone.
- **Fleets:** catch is assigned to one EwE fleet per species. The original
  comm/rec split misclassified several commercial fleets.

`tests/baseline/README.md` has details. To regenerate the original outputs
(needs the `xlsx` package, which needs Java, plus `plyr`):

```bash
Rscript legacy/BAM2EwE_age0_corr.R
```

## Open questions

- **Missing selectivity:** the GAJ video index (VID) has no matching
  selectivity in its rdat, so its selectivity summary is blank. (GAJ HB
  uses `sel.m.GR`, and VSN CVID uses `sel.m.CVT`.)
- **Fleet and index descriptions** in `config/` were inferred from BAM codes.
  Check the ones flagged "confirm" (e.g. VSN `cHTR`).
- **Citations:** `data/bam_rdats/README.md` still needs the full SEDAR
  citations.
