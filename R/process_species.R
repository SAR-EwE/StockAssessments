# Run every step for one species.

#' Process one species from its config row.
#'
#' Returns a list of year x stanza matrices (B, M, F, Z, landings, discards),
#' long-format catch by fleet, indices, index summary, multi-stanza
#' parameters, stanza start lengths, alternative M, and reconciliation checks.
process_species <- function(id, cfg) {
  sp_cfg <- cfg$species[cfg$species$species_id == id, ]
  if (nrow(sp_cfg) != 1) stop("Unknown species_id: ", id, call. = FALSE)
  bam <- read_bam(sp_cfg, cfg$settings$data_dir)

  res <- list(id = id, common_name = sp_cfg$common_name, file = sp_cfg$file,
              ages = bam$ages,
              units = bam$dat$info[grep("^units", names(bam$dat$info))])
  res$B <- stanza_biomass(bam)
  res$M <- stanza_m(bam)
  stanzas <- colnames(res$B)

  fleets <- cfg$fleets[cfg$fleets$species_id == id, ]
  catch <- catch_by_fleet(bam, fleets)
  catch$ewe_fleet <- fleets$ewe_fleet[match(catch$bam_fleet, fleets$bam_fleet)]
  res$catch_fleet <- catch
  res$landings <- catch_matrix(catch, "landings", bam$years, stanzas)
  res$discards <- catch_matrix(catch, "discards", bam$years, stanzas)

  # F = (landings + dead discards) / mid-year biomass; Z = F + M
  res$F <- (res$landings + res$discards) / res$B
  res$Z <- res$F + res$M

  res$indices <- index_series(bam, cfg$indices)
  res$index_summary <- index_summary(bam, cfg$indices, cfg$settings$index_cv)
  res$multistanza <- multistanza_params(bam)
  res$stanza_lengths <- stanza_start_lengths(bam)
  res$alt_m <- lorenzen_alt_m(bam)
  res$checks <- reconciliation_checks(bam, res)
  res$ecopath <- ecopath_inputs(res, cfg$settings$base_year, cfg$settings$model_area_km2)
  res
}
