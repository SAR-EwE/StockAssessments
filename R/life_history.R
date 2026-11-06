# Multi-stanza parameters and stanza start lengths.

#' vbK, Wmat, Winf, Wmat/Winf and proportion mature in each stanza.
#'
#' Wmat is the weight (lb) at the first age with female maturity >= 0.5 and
#' Winf the largest weight among mature ages, as in the legacy code.
multistanza_params <- function(bam) {
  a <- bam$dat$a.series
  mat <- a[[bam$cfg$maturity_col]]
  w_mature <- a$weight[mat >= 0.5]
  wmat <- min(w_mature)
  winf <- max(w_mature)
  pmat <- round(tapply(mat, bam$stanza, mean), 2)
  data.frame(
    par = c("vbK", "Wmat", "Winf", "Wmat.Winf", paste0("Pmat", names(pmat))),
    val = c(bam$dat$parms$K, wmat, winf, wmat / winf, unname(pmat)),
    stringsAsFactors = FALSE
  )
}

#' Length (mm, floored) at the start of each EwE stanza.
#'
#' The age-0 stanza starts at length 0. Other stanzas start at the smallest
#' mean length at age among the ages in that stanza.
stanza_start_lengths <- function(bam) {
  st <- bam$cfg$stanza_ages[[1]]
  a <- bam$dat$a.series
  size <- vapply(st, function(s) {
    if (s == 0) return(0)
    floor(min(a$length[bam$stanza == s]))
  }, numeric(1))
  data.frame(age_start = st, length_start_mm = size)
}
