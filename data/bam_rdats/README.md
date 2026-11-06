# BAM assessment output files (`.rdat`)

Output files from the Beaufort Assessment Model (BAM) for South Atlantic SEDAR
stock assessments. Each `.rdat` is a text file written with R's `dput()` and is
read with `dget()`.

These are copies of the files in the OneDrive folder
`SAR EwE UF/BAM2EwE/BAM rdats/`, copied 2026-11-05 (checksums verified).
Older versions superseded by updates (SEDAR 56 black sea bass, SEDAR 41 gray
triggerfish, non-FES SEDAR 53 red grouper, SEDAR 73 red snapper, non-FES
SEDAR 55 vermilion snapper) are not in the repo and stay on OneDrive under
`old BAM files/`.

| File | Species (`info$species` in file) | SEDAR | Model years | Ages | Used in EwE |
|---|---|---|---|---|---|
| `BSB_SEDAR76.rdat` | Black sea bass | 76 | 1978–2021 | 0–11 | yes |
| `BTLF_SEDAR50.rdat` | SA Blueline Tilefish | 50 | 1970–2015 | 1–20 | no |
| `GAG_SEDAR71.rdat` | SA Gag Grouper | 71 | 1962–2019 | 1–16 | yes |
| `GAJ_SEDAR59.rdat` | SA Greater Amberjack | 59 | 1980–2017 | 1–10 | yes |
| `GTF_SEDAR82.rdat` | SA gray triggerfish | 82 | 1982–2021 | 1–8 | yes |
| `GTLF_SEDAR66.rdat` | Tilefish (golden tilefish) | 66 | 1972–2018 | 1–20 | yes |
| `RGRP_SEDAR53_FESupdate.rdat` | SA Red Grouper | 53, FES update | 1976–2015 | 1–16 | yes |
| `RPGY_SEDAR60.rdat` | Red Porgy | 60 | 1972–2017 | 1–14 | yes |
| `RSN_SEDAR73_update.rdat` | SA Red Snapper | 73 update | 1950–2023 | 1–20 | yes |
| `SCMP_SEDAR68.rdat` | **SA Red Grouper** (mislabeled; this is scamp) | 68 | 1969–2017 | 1–20 | no |
| `SNWY_SEDAR36_update.rdat` | SA Snowy Grouper | 36 update | 1974–2018 | 1–25 | yes |
| `VSN_SEDAR55_FESupdate.rdat` | SA Vermilion Snapper | 55, FES update | 1946–2016 | 1–12 | yes |

To add or replace an rdat, follow "Updating to a new assessment" in the main
README, then update the table below.

## Notes

- **Scamp:** the species name stored in `SCMP_SEDAR68.rdat` is
  "SA Red Grouper". Code must identify species by file name, never by
  `info$species`.
- **Mutton snapper:** not included. The old SEDAR 15A update file is ASAP3
  output, not BAM. A new Stock Synthesis assessment will replace it.
- **Units differ between assessments** (gutted vs. whole weight, `1000 lb` vs.
  metric tons, and different spawning-output units). See `info$units.*` in each
  file, or `tests/baseline/BAM list.csv`.
- The files are stored byte-for-byte (`*.rdat -text -diff` in
  `.gitattributes`), keeping their original CRLF line endings.

## Citations

TODO: add the full SEDAR report citation (and the update/FES report, where
one applies) for each file.
