cfg <- load_config(repo_root)
problems <- validate_config(cfg)

test_that("config matches the rdat files with no errors", {
  errors <- problems[problems$severity == "error", ]
  expect_equal(nrow(errors), 0, info = paste(errors$species_id, errors$message, collapse = "\n"))
})

test_that("the only warnings are indices without a selectivity", {
  warn <- problems[problems$severity == "warning", ]
  expect_true(all(warn$check == "index_selectivity"), info = paste(warn$message, collapse = "\n"))
})

test_that("every fleet has an EwE fleet", {
  expect_false(any(is.na(cfg$fleets$ewe_fleet) | cfg$fleets$ewe_fleet == ""))
})

test_that("fleets.csv has one row per BAM fleet and one SEDAR number per species", {
  expect_false(anyDuplicated(paste(cfg$fleets$species_id, cfg$fleets$bam_fleet)) > 0)
  n_sedar <- tapply(cfg$fleets$sedar, cfg$fleets$species_id, function(x) length(unique(x)))
  expect_true(all(n_sedar == 1))
})

test_that("read_fleet_config rejects separate landings and discards rows for one fleet", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  fl <- utils::read.csv(file.path(repo_root, "config/fleets.csv"), stringsAsFactors = FALSE)
  dup <- fl[fl$species_id == "RSN" & fl$bam_fleet == "cH", ]
  fl$catch_type[fl$species_id == "RSN" & fl$bam_fleet == "cH"] <- "landings"
  dup$catch_type <- "discards"
  utils::write.csv(rbind(fl, dup), tmp, row.names = FALSE)
  expect_error(read_fleet_config(tmp), "duplicated species/fleet")
})

test_that("parse_stanza_ages handles single, multiple and empty values", {
  expect_equal(parse_stanza_ages(c("0,1,4", "0", NA, "")),
               list(c(0L, 1L, 4L), 0L, integer(0), integer(0)))
})

# ---- the validator catches broken config ------------------------------------

errors_for <- function(cfg_mod, id) {
  cfg_mod$species <- cfg_mod$species[cfg_mod$species$species_id == id, ]
  p <- validate_config(cfg_mod)
  p[p$severity == "error", ]
}

test_that("validator flags a wrong maturity column", {
  bad <- cfg
  bad$species$maturity_col[bad$species$species_id == "RPGY"] <- "mat.female"
  expect_true("maturity_col" %in% errors_for(bad, "RPGY")$check)
})

test_that("validator flags gutted landings with no conversion", {
  bad <- cfg
  bad$species$gutted_to_whole[bad$species$species_id == "GAG"] <- "none"
  expect_true("gutted_to_whole" %in% errors_for(bad, "GAG")$check)
})

test_that("validator flags wrong length-weight units", {
  bad <- cfg
  bad$species$lw_weight_unit[bad$species$species_id == "BSB"] <- "kg"
  expect_true("lw_units" %in% errors_for(bad, "BSB")$check)
})

test_that("validator flags a fleet missing from fleets.csv", {
  bad <- cfg
  bad$fleets <- bad$fleets[!(bad$fleets$species_id == "VSN" & bad$fleets$bam_fleet == "cHTR"), ]
  e <- errors_for(bad, "VSN")
  expect_true(any(e$check == "fleets" & grepl("cHTR", e$message)))
})

test_that("validator flags a 'both' fleet that has no discards in the rdat", {
  bad <- cfg
  bad$fleets$catch_type[bad$fleets$species_id == "SNWY" & bad$fleets$bam_fleet == "cH"] <- "both"
  e <- errors_for(bad, "SNWY")
  expect_true(any(e$check == "fleets" & grepl("discards.*cH", e$message)))
})

test_that("validator flags a SEDAR number that does not match the file", {
  bad <- cfg
  bad$fleets$sedar[bad$fleets$species_id == "RGRP"] <- 5
  expect_true("sedar" %in% errors_for(bad, "RGRP")$check)
})

test_that("validator flags a bad stanza definition", {
  bad <- cfg
  bad$species$stanza_ages[[which(bad$species$species_id == "GAJ")]] <- c(1L, 3L)
  expect_true("stanzas" %in% errors_for(bad, "GAJ")$check)
})

test_that("validator flags a selectivity name not in the rdat", {
  bad <- cfg
  bad$indices$selectivity[bad$indices$species_id == "RSN" & bad$indices$bam_index == "VID"] <- "sel.m.nope"
  expect_true("index_selectivity" %in% errors_for(bad, "RSN")$check)
})
