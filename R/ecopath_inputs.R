# Ecopath base-year inputs by stanza.

#' Ecopath inputs for one species in the base year.
#'
#' One row per stanza plus a "total" row: biomass (mt and t/km2), F, M, Z,
#' landings and discards (mt and t/km2), and two biomass-accumulation rates:
#' BA_1yr = B[y+1]/B[y] - 1, and BA_3yr = sum(B[y..y+2]) / sum(B[y-3..y-1]) - 1,
#' with the earlier window shortened near the first model year.
#' For the total row, F = (landings + discards) / biomass; M and Z are NA.
ecopath_inputs <- function(res, base_year, area_km2) {
  yrs <- as.numeric(rownames(res$B))
  i <- match(base_year, yrs)
  if (is.na(i) || i + 2 > length(yrs)) {
    warning(res$id, ": base year ", base_year, " (+2) is outside the model years", call. = FALSE)
    return(NULL)
  }
  before <- (i - min(3, i - 1)):(i - min(1, i - 1))
  after <- i:(i + 2)
  ba <- function(b) {
    b <- unname(b)
    c(BA_1yr = b[i + 1] / b[i] - 1, BA_3yr = sum(b[after]) / sum(b[before]) - 1)
  }

  B <- res$B; L <- res$landings; D <- res$discards
  stz <- colnames(B)
  per_stanza <- data.frame(
    species_id = res$id, stanza = stz,
    biomass_mt = B[i, ], F = res$F[i, ], M = res$M[i, ], Z = res$Z[i, ],
    landings_mt = L[i, ], discards_mt = D[i, ],
    BA_1yr = apply(B, 2, function(b) ba(b)[["BA_1yr"]]),
    BA_3yr = apply(B, 2, function(b) ba(b)[["BA_3yr"]]),
    stringsAsFactors = FALSE
  )
  b_tot <- rowSums(B)
  total <- data.frame(
    species_id = res$id, stanza = "total",
    biomass_mt = b_tot[i], F = (sum(L[i, ]) + sum(D[i, ])) / b_tot[i],
    M = if (length(stz) == 1) res$M[i, 1] else NA_real_,
    Z = if (length(stz) == 1) res$Z[i, 1] else NA_real_,
    landings_mt = sum(L[i, ]), discards_mt = sum(D[i, ]),
    BA_1yr = ba(b_tot)[["BA_1yr"]], BA_3yr = ba(b_tot)[["BA_3yr"]],
    stringsAsFactors = FALSE
  )
  out <- rbind(per_stanza, total)
  out$biomass_t_km2 <- out$biomass_mt / area_km2
  out$landings_t_km2 <- out$landings_mt / area_km2
  out$discards_t_km2 <- out$discards_mt / area_km2
  rownames(out) <- NULL
  out[, c("species_id", "stanza", "biomass_mt", "biomass_t_km2", "F", "M", "Z",
          "landings_mt", "landings_t_km2", "discards_mt", "discards_t_km2", "BA_1yr", "BA_3yr")]
}
