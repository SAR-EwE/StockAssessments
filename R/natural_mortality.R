# Alternative M by stanza from a Lorenzen (1996) weight-based M curve.
#
# M at monthly ages is 3 * W(g)^-0.288, with W from the von Bertalanffy
# growth curve and the length-weight relationship. A scaled version is
# rescaled so its mean over the fully selected ages equals the assessment M
# (config altM_target). Stanza values are the survivorship-weighted mean, the
# value at the stanza start age, and the value at the stanza mid-age.

#' t0 used for the growth curve: the assessment t0, or -0.5 when t0 < -1.
alt_m_t0 <- function(t0) if (t0 < -1) -0.5 else t0

#' Ages (years) where mean F is at least half of F at the oldest age.
fully_selected_ages <- function(F_age) {
  fbar <- colMeans(F_age)
  ages <- as.numeric(colnames(F_age))
  ages[fbar / utils::tail(fbar, 1) >= 0.5]
}

lorenzen_alt_m <- function(bam) {
  p <- bam$dat$parms
  cfg <- bam$cfg
  n_age_cols <- ncol(bam$dat$N.age)

  age_mo <- 0:(n_age_cols * 12)
  age_yr <- age_mo / 12
  len_mm <- p$Linf * (1 - exp(-p$K * (age_yr - alt_m_t0(p$t0))))
  len <- len_mm / MM_PER_UNIT[[cfg$lw_length_unit]]
  w_kg <- p$wgt.a * len^p$wgt.b * KG_PER_UNIT[[cfg$lw_weight_unit]]

  m_target <- switch(cfg$altM_target,
    M.constant = p$M.constant,
    M.oldest_age = utils::tail(bam$dat$a.series$M, 1)
  )
  m_target <- round(m_target, 3)

  sel_ages <- fully_selected_ages(bam$dat$F.age)
  in_window <- age_yr >= min(sel_ages) & age_yr <= max(sel_ages)

  m_yr <- 3 * (w_kg * 1000)^-0.288
  m_mo <- m_yr / 12
  m_mo_scaled <- (m_target / 12) * m_mo / mean(m_mo[in_window])
  m_yr_scaled <- m_mo_scaled * 12

  survship <- cumprod(c(1, exp(-m_mo[-length(m_mo)])))
  survship_scaled <- cumprod(c(1, exp(-m_mo_scaled[-length(m_mo_scaled)])))

  st <- cfg$stanza_ages[[1]]
  bounds <- c(st, max(age_yr))
  at_age <- function(x, age) x[which.min(abs(age_yr - age))]
  rows <- lapply(seq_along(st), function(i) {
    idx <- age_yr >= bounds[i] & age_yr < bounds[i + 1]
    mid <- bounds[i] + (bounds[i + 1] - bounds[i]) / 2
    data.frame(
      startage_mo = st[i] * 12, startage = st[i], midage = mid,
      M_wtavg_yr = stats::weighted.mean(m_yr[idx], w = survship[idx]),
      M_styr = at_age(m_yr, st[i]),
      M_midyr = at_age(m_yr, mid),
      M_wtavg_yr_scaled = stats::weighted.mean(m_yr_scaled[idx], w = survship_scaled[idx]),
      M_styr_scaled = at_age(m_yr_scaled, st[i]),
      M_midyr_scaled = at_age(m_yr_scaled, mid)
    )
  })
  out <- do.call(rbind, rows)
  attr(out, "M_target") <- m_target
  attr(out, "fully_selected_ages") <- range(sel_ages)
  out
}
