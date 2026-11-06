# Draft config rows for a new or updated BAM assessment.
#
#   sugg <- suggest_config("data/bam_rdats/RSN_SEDAR99.rdat", species_id = "RSN")
#   sugg                                   # review the three tables
#   write_suggested_config(sugg)           # replace RSN's rows in config/*.csv
#
# Settings that can be read from the rdat are filled in automatically:
# fleets and catch types, indices and matching selectivities, maturity column,
# gutted-to-whole method, and length-weight units. If the species is already
# in the config (an updated assessment), its stanza ages, include flag, name,
# alternative-M target, and fleet/index descriptions are carried over.
# Anything that can't be inferred is left blank. Check it, then run
# validate_config().

#' Find length-weight units that reproduce a.series$weight (lb), or NA.
guess_lw_units <- function(dat, tolerance = 0.02) {
  a <- dat$a.series
  best <- c(NA_character_, NA_character_)
  best_err <- Inf
  for (lu in names(MM_PER_UNIT)) for (wu in names(LB_PER_UNIT)) {
    w <- dat$parms$wgt.a * (a$length / MM_PER_UNIT[[lu]])^dat$parms$wgt.b * LB_PER_UNIT[[wu]]
    err <- max(abs(w / a$weight - 1))
    if (is.finite(err) && err < best_err) {
      best_err <- err
      best <- c(lu, wu)
    }
  }
  if (best_err > tolerance) return(c(NA_character_, NA_character_))
  best
}

#' Sector from a BAM fleet code: codes starting with "c" or "com" are commercial.
guess_sector <- function(code) {
  ifelse(grepl("^(c|com)", code, ignore.case = TRUE), "commercial", "recreational")
}

#' Draft species, fleet and index config rows for one rdat.
#'
#' `file` is a path to the rdat (normally in data/bam_rdats/). Returns a list
#' with data frames `species`, `fleets` and `indices`, in the column layout
#' of the config CSVs.
suggest_config <- function(file, species_id, common_name = NA, sedar = NA, root = ".") {
  dat <- dget(file)
  if (is.null(dat$a.series) || is.null(dat$CLD.est.mats)) {
    stop(basename(file), " is not a BAM rdat (no a.series / CLD.est.mats)", call. = FALSE)
  }
  cfg_dir <- file.path(root, "config")
  old_sp <- read_species_config(file.path(cfg_dir, "species.csv"))
  old_sp <- old_sp[old_sp$species_id == species_id, ]
  old_fl <- read_fleet_config(file.path(cfg_dir, "fleets.csv"))
  old_fl <- old_fl[old_fl$species_id == species_id, ]
  old_ix <- read_index_config(file.path(cfg_dir, "indices.csv"))
  old_ix <- old_ix[old_ix$species_id == species_id, ]
  existing <- nrow(old_sp) == 1

  if (is.na(sedar)) {
    m <- regmatches(basename(file), regexpr("SEDAR[0-9]+", basename(file)))
    sedar <- if (length(m)) as.integer(sub("SEDAR", "", m)) else NA_integer_
  }
  if (is.na(common_name)) common_name <- if (existing) old_sp$common_name else dat$info$species

  a <- dat$a.series
  mat_col <- intersect(c("mat.female", "mat.fem.endyr"), names(a))[1]
  gutted <- grepl("gutted", dat$info$units.landings, ignore.case = TRUE)
  lw <- guess_lw_units(dat)

  species <- data.frame(
    species_id = species_id,
    file = basename(file),
    common_name = common_name,
    include = if (existing) old_sp$include else FALSE,
    ewe_stanza_ages = if (existing) old_sp$ewe_stanza_ages else NA_character_,
    maturity_col = mat_col,
    gutted_to_whole = if (!gutted) "none" else if (existing && old_sp$gutted_to_whole != "none") old_sp$gutted_to_whole else "by_age",
    altM_target = if (existing) old_sp$altM_target else if (!is.null(dat$parms$M.constant)) "M.constant" else "M.oldest_age",
    lw_length_unit = lw[1],
    lw_weight_unit = lw[2],
    notes = if (existing && !is.na(old_sp$notes)) old_sp$notes else
      paste0("Drafted by suggest_config() on ", Sys.Date(), if (!existing) "; set ewe_stanza_ages and include" else ""),
    stringsAsFactors = FALSE
  )

  land <- rdat_fleet_codes(dat, "Lw")
  disc <- rdat_fleet_codes(dat, "Dw")
  codes <- union(land, disc)
  fleets <- data.frame(
    species_id = species_id,
    sedar = sedar,
    bam_fleet = codes,
    catch_type = ifelse(codes %in% land & codes %in% disc, "both", ifelse(codes %in% land, "landings", "discards")),
    sector = guess_sector(codes),
    ewe_fleet = species_id,
    description = old_fl$description[match(codes, old_fl$bam_fleet)],
    stringsAsFactors = FALSE
  )

  idx <- rdat_index_codes(dat)
  sel_names <- names(dat$sel.age)
  guess_sel <- function(i) {
    old <- old_ix$selectivity[match(i, old_ix$bam_index)]
    if (!is.na(old) && old %in% sel_names) return(old)
    cand <- c(paste0("sel.m.", i), paste0("sel.v.", i))
    hit <- sel_names[tolower(sel_names) %in% tolower(cand)]
    if (length(hit)) hit[1] else NA_character_
  }
  indices <- data.frame(
    species_id = species_id,
    bam_index = idx,
    selectivity = vapply(idx, guess_sel, character(1), USE.NAMES = FALSE),
    description = old_ix$description[match(idx, old_ix$bam_index)],
    notes = old_ix$notes[match(idx, old_ix$bam_index)],
    stringsAsFactors = FALSE
  )

  missing <- c(
    if (is.na(species$ewe_stanza_ages)) "ewe_stanza_ages",
    if (is.na(mat_col)) "maturity_col",
    if (anyNA(lw)) "length-weight units",
    if (is.na(sedar)) "sedar",
    if (anyNA(indices$selectivity)) paste("selectivity for index", paste(idx[is.na(indices$selectivity)], collapse = ", "))
  )
  if (length(missing)) message("Fill in by hand: ", paste(missing, collapse = "; "))
  if (gutted) message("Landings are gutted weight: check gutted_to_whole (by_age or mean_ratio_gt1)")
  list(species = species, fleets = fleets, indices = indices)
}

#' Replace one species' rows in config/species.csv, fleets.csv and indices.csv
#' with the rows from suggest_config(). Other species are untouched.
write_suggested_config <- function(sugg, root = ".") {
  id <- sugg$species$species_id
  cfg_dir <- file.path(root, "config")
  replace_rows <- function(name, new) {
    path <- file.path(cfg_dir, name)
    old <- read_csv_config(path)
    out <- rbind_fill(list(old[old$species_id != id, , drop = FALSE], new[, intersect(names(old), names(new))]))
    out <- out[order(out$species_id), names(old)]
    utils::write.csv(out, path, row.names = FALSE, na = "")
  }
  replace_rows("species.csv", sugg$species)
  replace_rows("fleets.csv", sugg$fleets)
  replace_rows("indices.csv", sugg$indices)
  invisible(TRUE)
}
