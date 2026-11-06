# Landings and dead discards (mt, whole weight) by fleet and stanza.

#' Whole:gutted weight ratio at each age (1 when landings are already whole).
gutted_to_whole_ratio <- function(bam) {
  a <- bam$dat$a.series
  method <- bam$cfg$gutted_to_whole
  if (method == "none") return(rep(1, nrow(a)))
  r <- a$weight / a$weight.gut
  switch(method,
    by_age = r,
    # gag: the ratio is < 1 at young ages, so use the mean of the ratios > 1
    mean_ratio_gt1 = rep(mean(r[r > 1]), length(r)),
    stop("Unknown gutted_to_whole method: ", method, call. = FALSE)
  )
}

#' Catch by fleet and stanza in long format.
#'
#' Uses the CLD.est.mats Lw.<fleet> / Dw.<fleet> matrices (1000 lb, in the
#' assessment's landings weight units), converts gutted to whole weight and
#' 1000 lb to mt, and sums ages within stanzas. Returns columns year,
#' bam_fleet, catch_type, stanza, catch_mt.
catch_by_fleet <- function(bam, fleets) {
  ratio <- gutted_to_whole_ratio(bam)
  cld <- bam$dat$CLD.est.mats
  out <- list()
  for (ct in c("landings", "discards")) {
    prefix <- if (ct == "landings") "Lw." else "Dw."
    fl <- fleets_with_catch(fleets[fleets$species_id == bam$id, ], ct)
    for (f in fl$bam_fleet) {
      if (is.null(cld[[paste0(prefix, f)]])) stop(bam$id, ": no ", prefix, f, " in CLD.est.mats", call. = FALSE)
      m <- as.matrix(cld[[paste0(prefix, f)]])
      m <- sweep(m, 2, ratio, "*") * KLB_TO_MT
      s <- sum_by_stanza(m, bam$stanza)
      out[[length(out) + 1]] <- data.frame(
        year = rep(bam$years, times = ncol(s)),
        bam_fleet = f,
        catch_type = ct,
        stanza = rep(colnames(s), each = nrow(s)),
        catch_mt = as.vector(s),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

#' Total catch (all fleets) of one catch type as a year x stanza matrix.
catch_matrix <- function(catch_long, catch_type, years, stanzas) {
  m <- matrix(0, nrow = length(years), ncol = length(stanzas),
              dimnames = list(years, stanzas))
  x <- catch_long[catch_long$catch_type == catch_type, ]
  if (nrow(x) == 0) return(m)
  tot <- tapply(x$catch_mt, list(x$year, x$stanza), sum)
  m[rownames(tot), colnames(tot)] <- tot
  m
}
