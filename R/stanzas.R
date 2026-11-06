# Map assessment ages to EwE stanzas and aggregate age matrices by stanza.

#' EwE stanza (its start age) for each assessment age.
#'
#' Each age goes to the last stanza that starts at or before it. If the
#' assessment starts at age 1, the age-0 stanza simply has no ages, and the
#' stanzas reported are those with at least one age (e.g. "0,1,4" with ages
#' 1-16 gives stanzas 1 and 4).
stanza_of_age <- function(ages, stanza_ages) {
  if (length(stanza_ages) == 0) stop("No EwE stanza ages given", call. = FALSE)
  vapply(ages, function(a) max(stanza_ages[stanza_ages <= a]), numeric(1))
}

stanza_labels <- function(stanza) paste0("age", sort(unique(stanza)))

#' Sum the columns (ages) of a year x age matrix within each stanza.
sum_by_stanza <- function(mat, stanza) {
  out <- t(rowsum(t(as.matrix(mat)), group = stanza))
  colnames(out) <- stanza_labels(stanza)
  out
}

#' Mid-year biomass (mt) by year and stanza: mid-year numbers x weight (mt).
stanza_biomass <- function(bam) {
  naa <- bam$dat$N.age.mdyr
  baa <- sweep(naa, MARGIN = 2, bam$dat$a.series$wgt.mt, "*")
  out <- sum_by_stanza(baa, bam$stanza)
  rownames(out) <- bam$years
  out
}

#' Natural mortality by year and stanza: M at age weighted by mid-year numbers.
stanza_m <- function(bam) {
  naa <- as.matrix(bam$dat$N.age.mdyr)
  m <- bam$dat$a.series$M
  num <- sum_by_stanza(sweep(naa, 2, m, "*"), bam$stanza)
  den <- sum_by_stanza(naa, bam$stanza)
  out <- num / den
  rownames(out) <- bam$years
  out
}
