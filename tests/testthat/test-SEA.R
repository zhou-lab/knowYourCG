## Offline tests for testEnrichmentSEA (issue #6): sign convention, symmetry
## of the p-value in the sign of the query, edge cases, and the plot helper.
## Synthetic data only, so no ExperimentHub access is needed.

## getDBs() names each set by its "dbname" attribute, not by its list name,
## so synthetic sets carry both, as a loaded knowledgebase does.
as_dbs <- function(dbs) {
    setNames(lapply(names(dbs), function(nm) {
        d <- dbs[[nm]]
        attr(d, "group") <- "synthetic"
        attr(d, "dbname") <- nm
        d
    }), names(dbs))
}

make_sea_data <- function(seed = 1, n = 2000, k = 60, shift = 1) {
    set.seed(seed)
    ids <- sprintf("cg%08d", seq_len(n))
    z <- setNames(rnorm(n), ids)
    up <- ids[seq_len(k)]                 # enriched at large values
    down <- ids[seq(k + 1, 2 * k)]        # enriched at small values
    null <- ids[seq(2 * k + 1, 3 * k)]    # no signal
    z[up] <- z[up] + shift
    z[down] <- z[down] - shift
    dbs <- as_dbs(list(up = up, down = down, null = null))
    list(z = z, dbs = dbs)
}

test_that("estimate sign follows GSEA: positive at large values", {
    d <- make_sea_data()
    set.seed(0)
    res <- testEnrichmentSEA(d$z, d$dbs, platform = "MM285", silent = TRUE)
    expect_s3_class(res, "data.frame")
    expect_setequal(res$dbname, c("up", "down", "null"))
    expect_gt(res$estimate[res$dbname == "up"], 0)
    expect_lt(res$estimate[res$dbname == "down"], 0)
    expect_lt(res$p.value[res$dbname == "up"], 0.01)
    expect_lt(res$p.value[res$dbname == "down"], 0.01)
    expect_gt(res$p.value[res$dbname == "null"], 0.05)
    expect_true(all(c("FDR", "test", "nQ", "nD", "overlap") %in%
                    colnames(res)))
    expect_equal(res$overlap[res$dbname == "up"], 60)
})

test_that("negating the query flips the estimate but keeps the p-value", {
    d <- make_sea_data(seed = 7)
    set.seed(0)
    a <- testEnrichmentSEA(d$z, d$dbs, platform = "MM285", silent = TRUE)
    set.seed(0)
    b <- testEnrichmentSEA(-d$z, d$dbs, platform = "MM285", silent = TRUE)
    a <- a[order(a$dbname), ]; b <- b[order(b$dbname), ]
    expect_equal(a$estimate, -b$estimate, tolerance = 1e-8)
    expect_equal(a$p.value, b$p.value, tolerance = 1e-8)
    expect_equal(a$log10.p.value, b$log10.p.value, tolerance = 1e-8)
})

test_that("p-value is reproducible with the same seed and in (0, 1]", {
    d <- make_sea_data(seed = 3)
    set.seed(11)
    a <- testEnrichmentSEA(d$z, d$dbs, platform = "MM285", silent = TRUE)
    set.seed(11)
    b <- testEnrichmentSEA(d$z, d$dbs, platform = "MM285", silent = TRUE)
    expect_equal(a, b)
    expect_true(all(a$p.value > 0 & a$p.value <= 1))
    expect_true(all(is.finite(a$log10.p.value)))
})

test_that("categorical query against numeric databases is the transpose", {
    d <- make_sea_data(seed = 5)
    set.seed(0)
    res <- testEnrichmentSEA(d$dbs$up, list(score = d$z),
                             platform = "MM285", silent = TRUE)
    expect_equal(nrow(res), 1)
    expect_gt(res$estimate, 0)
    expect_lt(res$p.value, 0.01)
})

test_that("query and databases must be one numeric and one categorical", {
    d <- make_sea_data()
    expect_error(testEnrichmentSEA(d$z, list(score = d$z),
                                   platform = "MM285", silent = TRUE),
                 "numerical and one categorical")
    expect_error(testEnrichmentSEA(d$dbs$up, d$dbs,
                                   platform = "MM285", silent = TRUE),
                 "numerical and one categorical")
    expect_error(testEnrichmentSEA(d$z, NULL, platform = "MM285"))
})

test_that("empty, disjoint and all-covering sets are handled", {
    d <- make_sea_data()
    dbs <- as_dbs(list(up = d$dbs$up,
                       empty = character(0),
                       disjoint = c("cg99999991", "cg99999992"),
                       everything = names(d$z)))
    set.seed(0)
    res <- suppressWarnings(testEnrichmentSEA(
        d$z, dbs, platform = "MM285", silent = TRUE))
    expect_false("empty" %in% res$dbname)          # dropped
    expect_equal(res$estimate[res$dbname == "disjoint"], 0)
    expect_equal(res$p.value[res$dbname == "disjoint"], 1)
    expect_equal(res$overlap[res$dbname == "disjoint"], 0)
    expect_equal(res$estimate[res$dbname == "everything"], 0)
    expect_equal(res$p.value[res$dbname == "everything"], 1)
    ## two warnings: that some set rows are absent, then how many were used
    w <- capture_warnings(testEnrichmentSEA(
        d$z, as_dbs(list(partial = c(d$dbs$up, "cg99999999"))),
        platform = "MM285", silent = TRUE))
    expect_match(w[1], "Not every data")
    expect_match(w[2], "Using 60 in 61")
})

test_that("precise = TRUE refines a significant p-value", {
    d <- make_sea_data(seed = 2, shift = 2)
    set.seed(0)
    res <- testEnrichmentSEA(d$z, d$dbs["up"], platform = "MM285",
                             silent = TRUE, precise = TRUE)
    expect_lt(res$p.value, 0.01)
    expect_gt(res$estimate, 0)
})

test_that("prepPlot returns the walk inputs and the plot builds", {
    d <- make_sea_data()
    set.seed(0)
    res <- testEnrichmentSEA(d$z, d$dbs, platform = "MM285",
                             silent = TRUE, prepPlot = TRUE)
    expect_type(res, "list")
    expect_named(res, c("up", "down", "null"))
    expect_true(all(c("res", "dCont", "dDisc") %in% names(res$up)))
    expect_equal(res$up$dDisc, d$dbs$up, ignore_attr = TRUE)
    p <- KYCG_plotSetEnrichment(res$up, n_sample = 200, n_presence = 20)
    expect_false(is.null(p))
    expect_error(KYCG_plotSetEnrichment(res$up$res), "dDisc")
})

test_that("messages are printed unless silent", {
    d <- make_sea_data()
    set.seed(0)
    expect_message(testEnrichmentSEA(d$z, d$dbs, platform = "MM285"),
                   "Testing against 3 database")
})
