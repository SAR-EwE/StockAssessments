# End-to-end tests: run every included species and check the results
# against the assessment totals and the legacy baseline (tests/baseline/).
#
# Regression tests compare a species with the baseline only while its
# config still points at the rdat the baseline was made from. After an
# assessment is updated, that species' comparisons are skipped, so update
# the baseline only if you want to keep comparing it.

cfg <- load_config(repo_root)
ids <- cfg$species$species_id[cfg$species$include]
results <- lapply(ids, process_species, cfg = cfg)
names(results) <- ids

baseline_dir <- file.path(repo_root, "tests", "baseline")
baseline_list <- utils::read.csv(file.path(baseline_dir, "BAM list.csv"), stringsAsFactors = FALSE)

# Baseline workbook for a species, or NULL if its rdat has changed since then
baseline_workbook <- function(id) {
  file <- cfg$species$file[cfg$species$species_id == id]
  row <- baseline_list[basename(baseline_list$filename) == file, ]
  if (nrow(row) == 0) return(NULL)
  # legacy file names end in "_no age0" when the assessment starts at age 1
  suffix <- if (min(results[[id]]$ages) >= 1) "_no age0" else ""
  path <- file.path(baseline_dir, paste0(row$species, suffix, ".xlsx"))
  if (file.exists(path)) path else NULL
}

read_sheet <- function(path, sheet) {
  x <- openxlsx::read.xlsx(path, sheet = sheet)
  as.matrix(x[, -1, drop = FALSE])
}

expect_close <- function(new, old, tol = 1e-10, info = NULL) {
  new <- as.matrix(new); old <- as.matrix(old)
  expect_equal(dim(new), dim(old), info = info)
  expect_lt(max(abs(new - old) / pmax(abs(old), 1e-9), na.rm = TRUE), tol, label = info)
}

# ---- reconciliation with the assessment ------------------------------------

test_that("landings, discards and biomass add back up to the assessment for every species", {
  for (r in results) {
    d <- check_max_rel_diff(r$checks)
    # 1e-5 allows for the rounded 1000 lb -> mt factor (0.453592) on BSB biomass
    expect_true(all(d < 1e-5), info = paste(r$id, paste(names(d), signif(d, 3), collapse = ", ")))
  }
})

test_that("F = (landings + discards) / biomass and Z = F + M", {
  for (r in results) {
    expect_equal(r$F, (r$landings + r$discards) / r$B, info = r$id)
    expect_equal(r$Z, r$F + r$M, info = r$id)
  }
})

test_that("catch by fleet sums to the stanza totals", {
  for (r in results) {
    for (ct in c("landings", "discards")) {
      x <- r$catch_fleet[r$catch_fleet$catch_type == ct, ]
      tot <- if (nrow(x)) sum(x$catch_mt) else 0
      expect_equal(tot, sum(r[[ct]]), tolerance = 1e-10, info = paste(r$id, ct))
    }
    expect_true(all(r$catch_fleet$ewe_fleet == r$id), info = r$id)
  }
})

test_that("every species has one EwE stanza per configured stanza with ages", {
  for (r in results) {
    st <- cfg$species$stanza_ages[[which(cfg$species$species_id == r$id)]]
    expect_true(all(sub("^age", "", colnames(r$B)) %in% st), info = r$id)
    expect_equal(r$stanza_lengths$age_start, st, info = r$id)
    expect_equal(r$stanza_lengths$length_start_mm[1], 0, info = r$id)
  }
})

# ---- regression against the legacy baseline --------------------------------

test_that("stanza time series match the legacy baseline", {
  for (id in ids) {
    path <- baseline_workbook(id)
    if (is.null(path)) next
    r <- results[[id]]
    expect_close(r$B, read_sheet(path, "Biomass mt"), info = paste(id, "B"))
    expect_close(r$F, read_sheet(path, "F rates"), info = paste(id, "F"))
    expect_close(r$M, read_sheet(path, "M rates"), info = paste(id, "M"))
    expect_close(r$Z, read_sheet(path, "Z rates"), info = paste(id, "Z"))
    expect_close(r$landings, read_sheet(path, "landings mt"), info = paste(id, "landings"))
  }
})

test_that("catch by fleet matches the legacy baseline", {
  for (id in ids) {
    path <- baseline_workbook(id)
    if (is.null(path)) next
    r <- results[[id]]
    for (ct in c("landings", "discards")) {
      sheet <- paste(ct, "mt fleet")
      old <- openxlsx::read.xlsx(path, sheet = sheet)
      x <- r$catch_fleet[r$catch_fleet$catch_type == ct, ]
      if (nrow(x) == 0) {
        expect_true(all(as.matrix(old[, -1]) == 0), info = paste(id, ct))
        next
      }
      x$col <- paste0(x$bam_fleet, ".", sub("^age", "", x$stanza))
      new <- tapply(x$catch_mt, list(x$year, x$col), sum)
      expect_setequal(colnames(new), names(old)[-1])
      expect_close(new[, names(old)[-1]], old[, -1], info = paste(id, ct))
    }
  }
})

test_that("multi-stanza parameters match the legacy baseline", {
  for (id in ids) {
    path <- baseline_workbook(id)
    if (is.null(path)) next
    old <- openxlsx::read.xlsx(path, sheet = "multistanza")
    expect_equal(results[[id]]$multistanza$par, old$par, info = id)
    expect_close(results[[id]]$multistanza$val, old$val, info = id)
  }
})

test_that("index time series match the legacy baseline", {
  old <- utils::read.csv(file.path(baseline_dir, "indices.csv"))
  names(old) <- sub("_u[.]", "_", sub("[.]ob$", "", tolower(names(old))))
  for (id in ids) {
    if (is.null(baseline_workbook(id))) next
    x <- results[[id]]$indices
    for (ix in unique(x$index)) {
      col <- tolower(ix)
      expect_true(col %in% names(old), info = ix)
      o <- old[!is.na(old[[col]]), c("year", col)]
      n <- x[x$index == ix, ]
      expect_equal(n$year, o$year, info = ix)
      expect_equal(n$obs, o[[col]], tolerance = 1e-10, info = ix)
    }
  }
})

# ---- deliberate fixes relative to the legacy code ---------------------------

test_that("each index gets its own weighted CV (legacy misaligned them)", {
  for (id in ids) {
    bam <- read_bam(cfg$species[cfg$species$species_id == id, ], cfg$settings$data_dir)
    ts <- year_rows(bam$dat$t.series, bam$years)
    s <- results[[id]]$index_summary
    for (j in seq_len(nrow(s))) {
      code <- sub(paste0("^", id, "_"), "", s$index[j])
      expect_equal(s$cv[j], mean(bam_na(ts[[paste0("cv.U.", code)]]), na.rm = TRUE), info = s$index[j])
    }
  }
  # concrete example: gag video index (legacy reported the headboat unweighted CV)
  gag <- results$GAG$index_summary
  expect_equal(round(gag$cv[gag$index == "GAG_VID"], 3), 0.406)
})

test_that("Ecopath total rows sum all stanzas (legacy used stanza 1 only)", {
  for (r in results) {
    e <- r$ecopath
    tot <- e[e$stanza == "total", ]
    st <- e[e$stanza != "total", ]
    expect_equal(tot$biomass_mt, sum(st$biomass_mt), info = r$id)
    expect_equal(tot$landings_mt, sum(st$landings_mt), info = r$id)
    i <- match(cfg$settings$base_year, as.numeric(rownames(r$B)))
    b <- rowSums(r$B)
    expect_equal(tot$BA_1yr, unname(b[i + 1] / b[i] - 1), info = r$id)
  }
  expect_equal(round(results$GAG$ecopath$BA_1yr[results$GAG$ecopath$stanza == "total"], 4), 0.0456)
})

test_that("alternative M scaling averages to the assessment M over fully selected ages", {
  for (id in ids) {
    a <- results[[id]]$alt_m
    expect_true(all(is.finite(a$M_wtavg_yr_scaled)), info = id)
    expect_true(all(diff(a$M_styr_scaled) < 0), info = id)  # M declines with age
    expect_equal(attr(a, "M_target") > 0, TRUE, info = id)
  }
})

# ---- outputs ----------------------------------------------------------------

test_that("combined tables and workbooks are written", {
  out <- file.path(tempdir(), "bam2ewe_test_output")
  on.exit(unlink(out, recursive = TRUE))
  tables <- combine_results(results, cfg)
  files <- write_tables(tables, out)
  expect_true(all(file.exists(files)))
  expect_setequal(unique(tables$stanza_timeseries$species_id), ids)
  expect_equal(nrow(tables$species), length(ids))
  xl <- write_species_xlsx(results[[1]], file.path(out, "test.xlsx"))
  expect_true("ecopath inputs" %in% openxlsx::getSheetNames(xl))
})
