# Relative abundance indices, their CVs, and selectivity summaries.

#' Index observations in long format: year, index, obs.
index_series <- function(bam, indices) {
  ix <- indices[indices$species_id == bam$id, ]
  ts <- year_rows(bam$dat$t.series, bam$years)
  out <- lapply(ix$bam_index, function(i) {
    obs <- bam_na(ts[[paste0("U.", i, ".ob")]])
    keep <- !is.na(obs)
    data.frame(year = bam$years[keep], index = paste0(bam$id, "_", i), obs = obs[keep],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Summarize a selectivity curve (vector by age, or year x age matrix
#' averaged over years): the value closest to 1 (selmax) and the first age
#' where it occurs, and on the ascending limb the values closest to 0.5 and
#' 0.95 and their ages.
selectivity_summary <- function(sel) {
  if (is.matrix(sel)) sel <- colMeans(sel, na.rm = TRUE)
  ages <- as.numeric(names(sel))
  i_max <- which.min(abs(1 - sel))
  asc <- sel[seq_len(i_max)]
  i50 <- which.min(abs(0.5 - asc))
  i95 <- which.min(abs(0.95 - asc))
  data.frame(Asel50 = ages[i50], Asel95 = ages[i95], Aselmax = ages[i_max],
             sel50 = unname(sel[i50]), sel95 = unname(sel[i95]), selmax = unname(sel[i_max]))
}

#' One row per index: mean CV and selectivity summary.
index_summary <- function(bam, indices, cv_type = "weighted") {
  ix <- indices[indices$species_id == bam$id, ]
  ts <- year_rows(bam$dat$t.series, bam$years)
  cv_prefix <- if (cv_type == "weighted") "cv.U." else "cv.unwgted.U."
  out <- lapply(seq_len(nrow(ix)), function(j) {
    i <- ix$bam_index[j]
    cv <- mean(bam_na(ts[[paste0(cv_prefix, i)]]), na.rm = TRUE)
    sel_name <- ix$selectivity[j]
    sel <- if (is.na(sel_name)) {
      data.frame(Asel50 = NA, Asel95 = NA, Aselmax = NA, sel50 = NA, sel95 = NA, selmax = NA)
    } else {
      selectivity_summary(bam$dat$sel.age[[sel_name]])
    }
    cbind(data.frame(index = paste0(bam$id, "_", i), selectivity = sel_name, cv = cv,
                     stringsAsFactors = FALSE), sel)
  })
  do.call(rbind, out)
}
