# Unit tests for individual calculations, using small made-up inputs.

test_that("stanza_of_age assigns ages to the last stanza starting at or before them", {
  expect_equal(stanza_of_age(0:6, c(0, 1, 4)), c(0, 1, 1, 1, 4, 4, 4))
  # assessment starting at age 1: no age-0 stanza
  expect_equal(stanza_of_age(1:5, c(0, 1, 4)), c(1, 1, 1, 4, 4))
  # single stanza
  expect_equal(stanza_of_age(1:3, 0), c(0, 0, 0))
  expect_equal(stanza_labels(c(1, 1, 4)), c("age1", "age4"))
})

test_that("sum_by_stanza sums ages within stanzas", {
  m <- matrix(1:12, nrow = 2, dimnames = list(c("2000", "2001"), 0:5))
  s <- sum_by_stanza(m, stanza_of_age(0:5, c(0, 1, 4)))
  expect_equal(colnames(s), c("age0", "age1", "age4"))
  expect_equal(unname(s[1, ]), c(1, 3 + 5 + 7, 9 + 11))
})

test_that("selectivity_summary finds the ascending-limb ages", {
  sel <- c(`1` = 0.1, `2` = 0.5, `3` = 0.9, `4` = 1, `5` = 0.8)
  s <- selectivity_summary(sel)
  expect_equal(c(s$Asel50, s$Asel95, s$Aselmax), c(2, 3, 4))
  expect_equal(c(s$sel50, s$sel95, s$selmax), c(0.5, 0.9, 1))
  # a year x age matrix is averaged over years first
  s2 <- selectivity_summary(rbind(sel, sel))
  expect_equal(s2, s)
})

test_that("alternative-M t0 rule", {
  expect_equal(alt_m_t0(-1.5), -0.5)
  expect_equal(alt_m_t0(-0.8), -0.8)
  expect_equal(alt_m_t0(-1), -1)
})

test_that("fully_selected_ages uses ages from column names", {
  F_age <- matrix(c(0.01, 0.2, 0.6, 1, 1), nrow = 2, ncol = 5, byrow = TRUE, dimnames = list(NULL, 0:4))
  expect_equal(fully_selected_ages(F_age), c(2, 3, 4))
})

test_that("bam_na replaces negative missing-value codes", {
  expect_equal(bam_na(c(1, -99999, NA, 0)), c(1, NA, NA, 0))
})

test_that("Ecopath BA totals use all stanzas", {
  B <- matrix(c(10, 11, 12, 13, 14,
                100, 90, 80, 70, 60), ncol = 2, dimnames = list(2000:2004, c("age1", "age4")))
  zero <- B * 0
  res <- list(id = "X", B = B, F = zero, M = zero + 0.2, Z = zero + 0.2, landings = zero, discards = zero)
  e <- ecopath_inputs(res, base_year = 2001, area_km2 = 10)
  tot <- e[e$stanza == "total", ]
  b <- rowSums(B)
  expect_equal(tot$BA_1yr, b[["2002"]] / b[["2001"]] - 1)
  # window before the base year is shortened near the first model year
  expect_equal(tot$BA_3yr, sum(b[c("2001", "2002", "2003")]) / b[["2000"]] - 1)
  expect_equal(e$BA_1yr[e$stanza == "age1"], 12 / 11 - 1)
  expect_equal(tot$biomass_t_km2, b[["2001"]] / 10)
})

test_that("ecopath_inputs warns when the base year is outside the model years", {
  B <- matrix(1, 3, 1, dimnames = list(2000:2002, "age0"))
  res <- list(id = "X", B = B, F = B, M = B, Z = B, landings = B, discards = B)
  expect_warning(expect_null(ecopath_inputs(res, base_year = 2002, area_km2 = 1)), "outside")
})

test_that("rbind_fill fills missing columns with NA", {
  out <- rbind_fill(list(data.frame(a = 1, b = 2), data.frame(a = 3, c = 4)))
  expect_equal(names(out), c("a", "b", "c"))
  expect_equal(out$c, c(NA, 4))
})
