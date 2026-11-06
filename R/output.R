# Combine per-species results into tidy tables and write them out.

#' Year x stanza matrices to long format.
stanza_long <- function(res) {
  vars <- c(biomass_mt = "B", F = "F", M = "M", Z = "Z", landings_mt = "landings", discards_mt = "discards")
  yrs <- as.numeric(rownames(res$B))
  stz <- colnames(res$B)
  out <- data.frame(species_id = res$id,
                    year = rep(yrs, times = length(stz)),
                    stanza = rep(stz, each = length(yrs)),
                    stringsAsFactors = FALSE)
  for (v in names(vars)) out[[v]] <- as.vector(res[[vars[[v]]]])
  out
}

add_id <- function(df, id) if (is.null(df) || nrow(df) == 0) df else cbind(species_id = id, df)

#' Catch summed over the BAM fleets assigned to each EwE fleet
#' (config/fleets.csv ewe_fleet): year, ewe_fleet, catch_type, stanza, catch_mt.
catch_by_ewe_fleet <- function(catch_fleet) {
  if (is.null(catch_fleet) || nrow(catch_fleet) == 0) return(NULL)
  out <- stats::aggregate(catch_mt ~ year + ewe_fleet + catch_type + stanza, data = catch_fleet, FUN = sum)
  out[order(out$ewe_fleet, out$catch_type, out$stanza, out$year), ]
}

#' rbind data frames whose columns differ, filling missing columns with NA.
rbind_fill <- function(dfs) {
  dfs <- Filter(function(d) !is.null(d) && nrow(d) > 0, dfs)
  cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    for (cc in setdiff(cols, names(d))) d[[cc]] <- NA
    d[, cols, drop = FALSE]
  })
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL
  out
}

#' Build the combined output tables from a list of process_species() results.
combine_results <- function(results, cfg) {
  bind <- function(f) rbind_fill(lapply(results, f))
  species_info <- rbind_fill(lapply(results, function(r) {
    fl <- cfg$fleets[cfg$fleets$species_id == r$id, ]
    data.frame(species_id = r$id, common_name = r$common_name, file = r$file,
               sedar = fl$sedar[1], first_year = min(as.numeric(rownames(r$B))),
               last_year = max(as.numeric(rownames(r$B))),
               stanzas = paste(colnames(r$B), collapse = ","),
               as.data.frame(r$units, stringsAsFactors = FALSE), stringsAsFactors = FALSE)
  }))
  idx_long <- bind(function(r) r$indices)
  idx_wide <- stats::reshape(idx_long, idvar = "year", timevar = "index", direction = "wide")
  names(idx_wide) <- sub("^obs\\.", "", names(idx_wide))
  idx_wide <- idx_wide[order(idx_wide$year), ]
  rownames(idx_wide) <- NULL

  list(
    species = species_info,
    stanza_timeseries = bind(stanza_long),
    catch_by_fleet = bind(function(r) add_id(r$catch_fleet, r$id)),
    catch_by_ewe_fleet = bind(function(r) add_id(catch_by_ewe_fleet(r$catch_fleet), r$id)),
    ecopath_inputs = bind(function(r) r$ecopath),
    indices = idx_long,
    indices_wide = idx_wide,
    index_summary = bind(function(r) add_id(r$index_summary, r$id)),
    multistanza_params = bind(function(r) add_id(r$multistanza, r$id)),
    stanza_lengths = bind(function(r) add_id(r$stanza_lengths, r$id)),
    alternative_M = bind(function(r) add_id(r$alt_m, r$id)),
    checks = bind(function(r) add_id(r$checks, r$id))
  )
}

#' Write each combined table to <out_dir>/<name>.csv
write_tables <- function(tables, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(tables)) {
    utils::write.csv(tables[[nm]], file.path(out_dir, paste0(nm, ".csv")), row.names = FALSE)
  }
  invisible(file.path(out_dir, paste0(names(tables), ".csv")))
}

#' Per-species workbook with the same sheets as the legacy output, plus
#' stanza lengths, index summary and alternative M.
write_species_xlsx <- function(res, path) {
  wb <- openxlsx::createWorkbook()
  add <- function(sheet, x, row_names = FALSE) {
    openxlsx::addWorksheet(wb, sheet)
    if (row_names) x <- cbind(year = as.numeric(rownames(x)), as.data.frame(x))
    openxlsx::writeData(wb, sheet, x)
  }
  fleet_wide <- function(ct) {
    x <- res$catch_fleet[res$catch_fleet$catch_type == ct, ]
    if (nrow(x) == 0) return(data.frame(year = as.numeric(rownames(res$B))))
    x$col <- paste0(x$bam_fleet, ".", sub("^age", "", x$stanza))
    w <- stats::reshape(x[, c("year", "col", "catch_mt")], idvar = "year", timevar = "col", direction = "wide")
    names(w) <- sub("^catch_mt\\.", "", names(w))
    w
  }
  add("Biomass mt", res$B, TRUE)
  add("F rates", res$F, TRUE)
  add("Z rates", res$Z, TRUE)
  add("M rates", res$M, TRUE)
  add("landings mt", res$landings, TRUE)
  add("discards mt", res$discards, TRUE)
  add("landings mt fleet", fleet_wide("landings"))
  add("discards mt fleet", fleet_wide("discards"))
  add("catch by ewe fleet", catch_by_ewe_fleet(res$catch_fleet))
  add("multistanza", res$multistanza)
  add("stanza lengths", res$stanza_lengths)
  add("index summary", res$index_summary)
  add("alternative M", res$alt_m)
  add("ecopath inputs", res$ecopath)
  add("check", res$checks)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}
