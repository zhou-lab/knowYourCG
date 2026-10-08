#!/usr/bin/env Rscript
## scripts/coverage.R -- line coverage of tests/testthat over R/ only.
##
##   Rscript scripts/coverage.R           measure, print per file, write coverage.json
##   Rscript scripts/coverage.R --check   measure and fail if coverage.json is stale
##
## coverage.json is a shields.io endpoint badge, the same scheme as kycg. src/
## is YAME vendored in, with its own suite and badge, so only R/ is counted.
## --check recomputes and fails on drift, so the README never shows a number
## the suite does not back.

library(covr)
args <- commandArgs(trailingOnly = TRUE)
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
pkg <- if (length(self)) normalizePath(file.path(dirname(self), "..")) else "."
badge <- file.path(pkg, "coverage.json")

cov <- package_coverage(pkg, type = "tests", quiet = TRUE)
tally <- tally_coverage(cov, by = "line")
tally <- tally[startsWith(tally$filename, "R/"), ]
pct <- 100 * mean(tally$value > 0)

per <- vapply(split(tally$value > 0, tally$filename), mean, numeric(1)) * 100
per <- sort(per)
print(data.frame(file = names(per), percent = round(per, 1)),
      row.names = FALSE)
cat(sprintf("R coverage: %.1f%% of %d lines\n", pct, nrow(tally)))

color <- if (pct >= 90) "brightgreen" else if (pct >= 75) "green" else
    if (pct >= 60) "yellowgreen" else if (pct >= 40) "yellow" else "orange"
new <- list(schemaVersion = 1L, label = "R coverage",
            message = sprintf("%.1f%%", pct), color = color)

if ("--check" %in% args) {
    old <- jsonlite::read_json(badge)
    old_pct <- as.numeric(sub("%", "", old$message))
    if (abs(old_pct - pct) > 0.5)
        stop(sprintf(paste("coverage.json says %s but the suite gives",
                           "%.1f%%; run Rscript scripts/coverage.R and",
                           "commit coverage.json"), old$message, pct))
    cat("coverage.json is up to date\n")
} else {
    jsonlite::write_json(new, badge, auto_unbox = TRUE, pretty = TRUE)
    cat("wrote", badge, "\n")
}
