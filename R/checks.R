# Reconciliation checks: do the derived totals add back up to the assessment?

#' Year-by-year comparison of assessment totals with the pipeline's sums.
#'
#' Landings and discards are compared in the assessment's own units (1000 lb,
#' before the gutted-to-whole conversion). Biomass compares rowSums(B.age)
#' with Jan-1 numbers x weight (mt); mid-year biomass is reported for reference.
reconciliation_checks <- function(bam, res) {
  d <- bam$dat
  ts <- year_rows(d$t.series, bam$years)
  cld <- d$CLD.est.mats
  sum_cld <- function(prefix) {
    nm <- grep(paste0("^", prefix, "\\."), names(cld), value = TRUE)
    nm <- setdiff(nm, paste0(prefix, ".total"))
    if (length(nm) == 0) return(rep(0, length(bam$years)))
    Reduce(`+`, lapply(nm, function(n) rowSums(as.matrix(cld[[n]]))))
  }
  naa_jan1 <- as.matrix(year_rows(d$N.age, bam$years))
  # B.age is in the assessment's biomass units (mt, or 1000 lb for some)
  b_to_mt <- switch(d$info$units.biomass,
    "metric tons" = 1,
    "1000 lb" = KLB_TO_MT,
    stop(bam$id, ": unexpected units.biomass '", d$info$units.biomass, "'", call. = FALSE)
  )
  out <- data.frame(
    year = bam$years,
    B_age_sum_mt = rowSums(as.matrix(year_rows(d$B.age, bam$years))) * b_to_mt,
    B_jan1_calc_mt = rowSums(sweep(naa_jan1, 2, d$a.series$wgt.mt, "*")),
    B_midyr_calc_mt = rowSums(res$B),
    L_ts_klb = ts$total.L.klb,
    L_fleet_sum_klb = sum_cld("Lw")
  )
  if (!is.null(ts$total.D.klb)) {
    out$D_ts_klb <- ts$total.D.klb
    out$D_fleet_sum_klb <- sum_cld("Dw")
  }
  out
}

#' Largest relative difference between paired check columns.
check_max_rel_diff <- function(chk) {
  rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-6), na.rm = TRUE)
  out <- c(biomass = rel(chk$B_jan1_calc_mt, chk$B_age_sum_mt),
           landings = rel(chk$L_fleet_sum_klb, chk$L_ts_klb))
  if (!is.null(chk$D_ts_klb)) out["discards"] <- rel(chk$D_fleet_sum_klb, chk$D_ts_klb)
  out
}
