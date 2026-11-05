# Runs before the tests: source the pipeline code from the repo root.
# Run all tests from the repo root with:
#   Rscript -e 'testthat::test_dir("tests/testthat")'
repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/")
for (f in list.files(file.path(repo_root, "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}
