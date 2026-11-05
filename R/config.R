# Read and validate the pipeline configuration in config/.
#
#   cfg <- load_config()                 # from the repo root
#   problems <- validate_config(cfg)     # checks the config against the rdats
#
# Species are identified by species_id / file name, never by info$species
# inside the rdat (e.g. SCMP_SEDAR68.rdat is labeled "SA Red Grouper").

# Allowed values for coded columns
CONFIG_CHOICES <- list(
  gutted_to_whole = c("none", "by_age", "mean_ratio_gt1"),
  altM_target     = c("M.constant", "M.oldest_age"),
  lw_length_unit  = c("mm", "cm"),
  lw_weight_unit  = c("g", "kg", "mt"),
  catch_type      = c("landings", "discards", "both"),
  sector          = c("commercial", "recreational"),
  index_cv        = c("weighted", "unweighted")
)

# Pounds per unit of weight, for checking length-weight parameters
LB_PER_UNIT <- c(g = 1 / 453.59237, kg = 2.2046226, mt = 2204.6226)
# Millimetres per unit of length
MM_PER_UNIT <- c(mm = 1, cm = 10)

read_csv_config <- function(path) {
  if (!file.exists(path)) stop("Config file not found: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "",
                  strip.white = TRUE, check.names = FALSE)
}

check_choices <- function(df, col, file) {
  bad <- setdiff(unique(df[[col]]), CONFIG_CHOICES[[col]])
  if (length(bad) > 0) {
    stop(sprintf("%s: invalid %s value(s): %s (allowed: %s)", file, col,
                 paste(bad, collapse = ", "),
                 paste(CONFIG_CHOICES[[col]], collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

#' Parse "0,1,4" into c(0, 1, 4); empty/NA -> integer(0)
parse_stanza_ages <- function(x) {
  lapply(x, function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return(integer(0))
    as.integer(trimws(strsplit(s, ",", fixed = TRUE)[[1]]))
  })
}

read_species_config <- function(path = "config/species.csv") {
  sp <- read_csv_config(path)
  required <- c("species_id", "file", "common_name", "include", "ewe_stanza_ages",
                "maturity_col", "gutted_to_whole", "altM_target",
                "lw_length_unit", "lw_weight_unit")
  missing <- setdiff(required, names(sp))
  if (length(missing) > 0) stop(path, ": missing column(s) ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(sp$species_id)) stop(path, ": duplicated species_id", call. = FALSE)
  if (anyDuplicated(sp$file)) stop(path, ": duplicated file", call. = FALSE)
  sp$include <- as.logical(sp$include)
  if (anyNA(sp$include)) stop(path, ": include must be TRUE or FALSE", call. = FALSE)
  for (col in c("gutted_to_whole", "altM_target", "lw_length_unit", "lw_weight_unit")) {
    check_choices(sp, col, path)
  }
  sp$stanza_ages <- parse_stanza_ages(sp$ewe_stanza_ages)
  sp
}

read_fleet_config <- function(path = "config/fleets.csv") {
  fl <- read_csv_config(path)
  for (col in c("catch_type", "sector")) check_choices(fl, col, path)
  # One row per BAM fleet; a fleet with both landings and discards uses "both"
  key <- paste(fl$species_id, fl$bam_fleet)
  if (anyDuplicated(key)) {
    stop(path, ": duplicated species/fleet row(s) (use catch_type 'both' for a fleet with landings and discards): ",
         paste(unique(key[duplicated(key)]), collapse = ", "), call. = FALSE)
  }
  if (anyNA(fl$ewe_fleet)) stop(path, ": every fleet needs an ewe_fleet", call. = FALSE)
  if (anyNA(fl$sedar)) stop(path, ": every fleet needs a sedar number", call. = FALSE)
  fl
}

#' Fleets in the config that have a given catch type ("landings" or "discards")
fleets_with_catch <- function(fleets, catch) {
  fleets[fleets$catch_type %in% c(catch, "both"), ]
}

read_index_config <- function(path = "config/indices.csv") {
  ix <- read_csv_config(path)
  key <- paste(ix$species_id, ix$bam_index)
  if (anyDuplicated(key)) stop(path, ": duplicated species/index row(s)", call. = FALSE)
  ix
}

read_settings <- function(path = "config/settings.yml") {
  if (!file.exists(path)) stop("Settings file not found: ", path, call. = FALSE)
  s <- yaml::read_yaml(path)
  check_choices(data.frame(index_cv = s$index_cv), "index_cv", path)
  s
}

#' Load all configuration. `root` is the repo root.
load_config <- function(root = ".") {
  cfg_dir <- file.path(root, "config")
  settings <- read_settings(file.path(cfg_dir, "settings.yml"))
  settings$data_dir <- file.path(root, settings$data_dir)
  cfg <- list(
    settings = settings,
    species  = read_species_config(file.path(cfg_dir, "species.csv")),
    fleets   = read_fleet_config(file.path(cfg_dir, "fleets.csv")),
    indices  = read_index_config(file.path(cfg_dir, "indices.csv"))
  )
  unknown <- setdiff(c(cfg$fleets$species_id, cfg$indices$species_id), cfg$species$species_id)
  if (length(unknown) > 0) {
    stop("fleets.csv/indices.csv reference unknown species_id: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  cfg
}

# ---- validation against the rdat files --------------------------------------

#' BAM fleet codes in CLD.est.mats for a catch type ("Lw" landings, "Dw" discards)
rdat_fleet_codes <- function(dat, prefix = c("Lw", "Dw")) {
  prefix <- match.arg(prefix)
  nm <- names(dat$CLD.est.mats)
  codes <- sub(paste0("^", prefix, "\\."), "", nm[startsWith(nm, paste0(prefix, "."))])
  setdiff(codes, "total")
}

#' Index codes in t.series, e.g. "U.cH.ob" -> "cH"
rdat_index_codes <- function(dat) {
  nm <- grep("^U\\..*\\.ob$", names(dat$t.series), value = TRUE)
  sub("\\.ob$", "", sub("^U\\.", "", nm))
}

#' Check the configuration against the rdat files.
#'
#' Returns a data frame of problems with columns species_id, severity
#' ("error" or "warning"), check and message. Zero rows means everything
#' matches. Species with include = FALSE are checked too.
validate_config <- function(cfg, lw_tolerance = 0.02) {
  out <- list()
  add <- function(id, severity, check, msg) {
    out[[length(out) + 1]] <<- data.frame(species_id = id, severity = severity,
                                          check = check, message = msg,
                                          stringsAsFactors = FALSE)
  }
  cv_prefix <- if (cfg$settings$index_cv == "weighted") "cv.U." else "cv.unwgted.U."

  for (i in seq_len(nrow(cfg$species))) {
    sp <- cfg$species[i, ]
    id <- sp$species_id
    path <- file.path(cfg$settings$data_dir, sp$file)
    if (!file.exists(path)) {
      add(id, "error", "file", paste("rdat not found:", path))
      next
    }
    dat <- dget(path)
    a <- dat$a.series
    if (is.null(a) || is.null(dat$CLD.est.mats)) {
      add(id, "error", "format", "not a BAM rdat (no a.series / CLD.est.mats)")
      next
    }

    # stanza ages
    ages <- sp$stanza_ages[[1]]
    if (sp$include) {
      if (length(ages) == 0) {
        add(id, "error", "stanzas", "included species has no ewe_stanza_ages")
      } else if (ages[1] != 0 || is.unsorted(ages, strictly = TRUE) || max(ages) > max(a$age)) {
        add(id, "error", "stanzas", paste0("ewe_stanza_ages must start at 0, increase, and be <= max age (",
                                          max(a$age), "): ", sp$ewe_stanza_ages))
      }
    }

    # maturity column
    if (!sp$maturity_col %in% names(a)) {
      add(id, "error", "maturity_col", paste0("'", sp$maturity_col, "' not in a.series"))
    }

    # gutted -> whole conversion
    gutted <- grepl("gutted", dat$info$units.landings, ignore.case = TRUE)
    if (gutted && sp$gutted_to_whole == "none") {
      add(id, "error", "gutted_to_whole", paste0("landings are '", dat$info$units.landings,
                                                 "' but gutted_to_whole is 'none'"))
    }
    if (!gutted && sp$gutted_to_whole != "none") {
      add(id, "error", "gutted_to_whole", paste0("landings are '", dat$info$units.landings,
                                                 "' but gutted_to_whole is '", sp$gutted_to_whole, "'"))
    }
    if (sp$gutted_to_whole != "none" && is.null(a$weight.gut)) {
      add(id, "error", "gutted_to_whole", "a.series has no weight.gut column")
    }

    # alternative-M target
    if (sp$altM_target == "M.constant" && is.null(dat$parms$M.constant)) {
      add(id, "error", "altM_target", "parms$M.constant is missing")
    }

    # length-weight units: a * L^b (converted to lb) should match a.series$weight
    len <- a$length / MM_PER_UNIT[[sp$lw_length_unit]]
    w_lb <- dat$parms$wgt.a * len^dat$parms$wgt.b * LB_PER_UNIT[[sp$lw_weight_unit]]
    rel <- max(abs(w_lb / a$weight - 1))
    if (!is.finite(rel) || rel > lw_tolerance) {
      add(id, "error", "lw_units", sprintf("a*L^b in %s->%s differs from a.series$weight by up to %.1f%%",
                                           sp$lw_length_unit, sp$lw_weight_unit, 100 * rel))
    }

    # fleets
    sp_fleets <- cfg$fleets[cfg$fleets$species_id == id, ]
    sedars <- unique(sp_fleets$sedar)
    matches_file <- vapply(sedars, function(s) grepl(paste0("SEDAR", s, "(_|\\.)"), sp$file), logical(1))
    bad_sedar <- sedars[!matches_file]
    if (length(bad_sedar)) {
      add(id, "error", "sedar", paste0("fleets.csv sedar ", paste(bad_sedar, collapse = ", "),
                                       " does not match file ", sp$file))
    }
    for (ct in c("landings", "discards")) {
      prefix <- if (ct == "landings") "Lw" else "Dw"
      in_rdat <- rdat_fleet_codes(dat, prefix)
      in_cfg <- fleets_with_catch(sp_fleets, ct)$bam_fleet
      miss <- setdiff(in_rdat, in_cfg)
      extra <- setdiff(in_cfg, in_rdat)
      if (length(miss)) add(id, "error", "fleets", paste0(ct, " fleet(s) in rdat but not fleets.csv: ", paste(miss, collapse = ", ")))
      if (length(extra)) add(id, "error", "fleets", paste0(ct, " fleet(s) in fleets.csv but not rdat: ", paste(extra, collapse = ", ")))
    }

    # indices
    in_rdat <- rdat_index_codes(dat)
    ix <- cfg$indices[cfg$indices$species_id == id, ]
    miss <- setdiff(in_rdat, ix$bam_index)
    extra <- setdiff(ix$bam_index, in_rdat)
    if (length(miss)) add(id, "error", "indices", paste("index(es) in rdat but not indices.csv:", paste(miss, collapse = ", ")))
    if (length(extra)) add(id, "error", "indices", paste("index(es) in indices.csv but not rdat:", paste(extra, collapse = ", ")))
    for (j in seq_len(nrow(ix))) {
      sel <- ix$selectivity[j]
      if (is.na(sel)) {
        add(id, "warning", "index_selectivity", paste0("no selectivity set for index ", ix$bam_index[j]))
      } else if (!sel %in% names(dat$sel.age)) {
        add(id, "error", "index_selectivity", paste0("'", sel, "' (index ", ix$bam_index[j], ") not in sel.age"))
      }
      cv_col <- paste0(cv_prefix, ix$bam_index[j])
      if (ix$bam_index[j] %in% in_rdat && !cv_col %in% names(dat$t.series)) {
        add(id, "error", "index_cv", paste0("CV column '", cv_col, "' not in t.series"))
      }
    }
    rm(dat)
  }

  if (length(out) == 0) {
    return(data.frame(species_id = character(), severity = character(),
                      check = character(), message = character()))
  }
  do.call(rbind, out)
}
