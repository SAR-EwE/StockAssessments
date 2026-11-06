# Process BAM stock assessment outputs into Ecopath with Ecosim inputs.
#
# Run from the repo root:
#   Rscript run_bam2ewe.R
#
# Reads config/ and data/bam_rdats/, and writes:
#   output/*.csv              combined tidy tables for all species
#   output/xlsx/<id>.xlsx     one workbook per species

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- load_config(".")

# Stop if the config doesn't match the rdats
problems <- validate_config(cfg)
for (i in seq_len(nrow(problems))) {
  message(sprintf("[%s] %s %s: %s", problems$severity[i], problems$species_id[i],
                  problems$check[i], problems$message[i]))
}
if (any(problems$severity == "error")) stop("Config validation failed; see messages above.")

ids <- cfg$species$species_id[cfg$species$include]
results <- lapply(ids, function(id) {
  message("Processing ", id)
  process_species(id, cfg)
})
names(results) <- ids

# Reconciliation: derived totals must add back up to the assessment
for (r in results) {
  d <- check_max_rel_diff(r$checks)
  if (any(d > 1e-4)) {
    warning(r$id, ": reconciliation differences > 0.01%: ",
            paste(names(d), signif(d, 3), sep = "=", collapse = ", "), call. = FALSE)
  }
}

out_dir <- "output"
tables <- combine_results(results, cfg)
write_tables(tables, out_dir)
dir.create(file.path(out_dir, "xlsx"), showWarnings = FALSE)
for (r in results) write_species_xlsx(r, file.path(out_dir, "xlsx", paste0(r$id, ".xlsx")))

message("Wrote ", length(tables), " tables and ", length(results), " workbooks to ", out_dir, "/")
