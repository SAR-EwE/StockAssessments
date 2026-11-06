# suggest_config() should reproduce the hand-checked config for the current
# rdats, and write_suggested_config() should only touch one species' rows.

cfg <- load_config(repo_root)

test_that("suggest_config reproduces the current config for every rdat", {
  for (i in seq_len(nrow(cfg$species))) {
    sp <- cfg$species[i, ]
    s <- suppressMessages(suggest_config(file.path(cfg$settings$data_dir, sp$file), sp$species_id, root = repo_root))
    cols <- c("maturity_col", "gutted_to_whole", "altM_target", "lw_length_unit",
              "lw_weight_unit", "ewe_stanza_ages", "include", "common_name")
    for (cc in cols) expect_identical(s$species[[cc]], sp[[cc]], info = paste(sp$species_id, cc))

    fl <- cfg$fleets[cfg$fleets$species_id == sp$species_id, ]
    expect_setequal(s$fleets$bam_fleet, fl$bam_fleet)
    m <- match(fl$bam_fleet, s$fleets$bam_fleet)
    expect_equal(s$fleets$catch_type[m], fl$catch_type, info = sp$species_id)
    expect_equal(s$fleets$sector[m], fl$sector, info = sp$species_id)
    expect_equal(s$fleets$sedar[m], fl$sedar, info = sp$species_id)

    ix <- cfg$indices[cfg$indices$species_id == sp$species_id, ]
    m <- match(ix$bam_index, s$indices$bam_index)
    expect_equal(s$indices$selectivity[m], ix$selectivity, info = sp$species_id)
  }
})

test_that("suggest_config fills what it can for a species not in the config", {
  s <- suppressMessages(suggest_config(file.path(cfg$settings$data_dir, "BSB_SEDAR76.rdat"), "NEW", root = repo_root))
  expect_false(s$species$include)
  expect_true(is.na(s$species$ewe_stanza_ages))
  expect_equal(c(s$species$lw_length_unit, s$species$lw_weight_unit), c("mm", "g"))
  expect_equal(unique(s$fleets$sedar), 76L)
  expect_true(all(s$fleets$ewe_fleet == "NEW"))
  # CVID uses the MARMAP chevron trap selectivity, which can't be inferred
  expect_true(is.na(s$indices$selectivity[s$indices$bam_index == "CVID"]))
  expect_message(suggest_config(file.path(cfg$settings$data_dir, "BSB_SEDAR76.rdat"), "NEW", root = repo_root),
                 "ewe_stanza_ages")
})

test_that("write_suggested_config replaces only that species' rows", {
  root <- file.path(tempdir(), "bam2ewe_cfg_test")
  on.exit(unlink(root, recursive = TRUE))
  dir.create(file.path(root, "config"), recursive = TRUE)
  file.copy(list.files(file.path(repo_root, "config"), full.names = TRUE), file.path(root, "config"))

  s <- suppressMessages(suggest_config(file.path(cfg$settings$data_dir, "RSN_SEDAR73_update.rdat"), "RSN", root = root))
  s$fleets$description[s$fleets$bam_fleet == "GR"] <- "changed"
  s$fleets$ewe_fleet[s$fleets$bam_fleet == "cH"] <- "RSN_comm"
  write_suggested_config(s, root = root)

  fl <- read_fleet_config(file.path(root, "config", "fleets.csv"))
  expect_equal(fl$description[fl$species_id == "RSN" & fl$bam_fleet == "GR"], "changed")

  # a custom EwE fleet assignment survives re-running suggest_config()
  s2 <- suppressMessages(suggest_config(file.path(cfg$settings$data_dir, "RSN_SEDAR73_update.rdat"), "RSN", root = root))
  expect_equal(s2$fleets$ewe_fleet[s2$fleets$bam_fleet == "cH"], "RSN_comm")
  expect_equal(s2$fleets$ewe_fleet[s2$fleets$bam_fleet == "GR"], "RSN")
  others <- function(d) d[d$species_id != "RSN", ]
  key <- function(d) paste(d$species_id, d[[2]], d[[3]])
  old_fl <- others(cfg$fleets); new_fl <- others(fl)
  expect_equal(new_fl[order(key(new_fl)), ], old_fl[order(key(old_fl)), ], ignore_attr = TRUE)
  expect_equal(nrow(read_species_config(file.path(root, "config", "species.csv"))), nrow(cfg$species))
})
