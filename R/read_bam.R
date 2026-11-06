# Read a BAM .rdat file and pull out the pieces the pipeline uses.

KLB_TO_MT <- 0.453592          # 1000 lb -> metric tons
KG_PER_UNIT <- c(g = 0.001, kg = 1, mt = 1000)

#' Read one species' rdat and attach its config row.
#'
#' Returns a list with the raw rdat (`dat`), the config row (`cfg`), the
#' model years and ages, and the stanza assignment for each age.
read_bam <- function(sp_cfg, data_dir) {
  stopifnot(nrow(sp_cfg) == 1)
  path <- file.path(data_dir, sp_cfg$file)
  dat <- dget(path)
  if (is.null(dat$a.series) || is.null(dat$CLD.est.mats)) {
    stop(sp_cfg$file, " is not a BAM rdat (no a.series / CLD.est.mats)", call. = FALSE)
  }
  if (!grepl("^1000 lb", dat$info$units.landings)) {
    stop(sp_cfg$file, ": expected landings in '1000 lb', got '", dat$info$units.landings, "'", call. = FALSE)
  }
  years <- dat$parms$styr:dat$parms$endyr
  ages <- dat$a.series$age
  if (nrow(dat$N.age.mdyr) != length(years)) {
    stop(sp_cfg$file, ": N.age.mdyr rows do not match model years", call. = FALSE)
  }
  list(
    id = sp_cfg$species_id,
    cfg = sp_cfg,
    dat = dat,
    years = years,
    ages = ages,
    stanza = stanza_of_age(ages, sp_cfg$stanza_ages[[1]])
  )
}

#' Rows of an age matrix or t.series that belong to the model years
#' (some BAM matrices carry an extra terminal-year row).
year_rows <- function(x, years) {
  rn <- if (is.data.frame(x) && "year" %in% names(x)) x$year else as.numeric(rownames(x))
  idx <- match(years, rn)
  if (anyNA(idx)) stop("Model years missing from matrix/table", call. = FALSE)
  x[idx, , drop = FALSE]
}

#' Replace BAM missing-value codes (negative, e.g. -99999) with NA
bam_na <- function(x) {
  x[!is.na(x) & x < 0] <- NA
  x
}
