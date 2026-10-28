## Offline tests for the 1.9.6 warnings and the SEA nQ/nD order. Synthetic
## data and the package's own title table only, so no ExperimentHub access.

as_dbs <- function(dbs) {
    setNames(lapply(names(dbs), function(nm) {
        d <- dbs[[nm]]
        attr(d, "group") <- "synthetic"
        attr(d, "dbname") <- nm
        d
    }), names(dbs))
}

test_that("a query with no probe in the universe warns", {
    u <- sprintf("cg%08d", seq_len(1000))
    dbs <- as_dbs(list(a = u[1:100]))
    expect_warning(
        testEnrichment(sprintf("ch%08d", 1:20), dbs, universe = u,
                       platform = "MM285", silent = TRUE),
        "None of the 20 query probes is in the universe")
    expect_no_warning(
        testEnrichment(u[1:20], dbs, universe = u, platform = "MM285",
                       silent = TRUE))
})

test_that("an empty database is dropped with a warning naming it", {
    u <- sprintf("cg%08d", seq_len(1000))
    dbs <- as_dbs(list(a = u[1:100], none = character(0)))
    expect_warning(
        res <- testEnrichment(u[1:20], dbs, universe = u, platform = "MM285",
                              silent = TRUE),
        "Empty database\\(s\\), dropped: none")
    expect_equal(res$dbname, "a")
})

test_that("a database group for another platform is dropped with a warning", {
    expect_warning(
        nms <- knowYourCG:::guess_dbnames("KYCG.MM285.chromHMM.20210210",
                                          platform = "HM450", silent = TRUE),
        "No database group for platform HM450; dropped: KYCG.MM285.chromHMM")
    expect_length(nms, 0)
    ## a name that matches several platforms picks one, silently
    expect_no_warning(knowYourCG:::guess_dbnames(
        "KYCG.MM285.chromHMM.20210210", platform = "MM285", silent = TRUE))
})

test_that("an unknown database group is dropped with a warning", {
    expect_warning(res <- getDBs("no_such_group_xyz", silent = TRUE),
                   "Unknown database group\\(s\\), dropped: no_such_group_xyz")
    expect_null(res)
})

test_that("SEA reports nQ as the query and nD as the database", {
    set.seed(1)
    ids <- sprintf("cg%08d", seq_len(500))
    z <- setNames(rnorm(500), ids)
    dbs <- as_dbs(list(s = ids[1:40]))
    res <- testEnrichmentSEA(z, dbs, platform = "MM285", silent = TRUE)
    expect_equal(res$nQ, 500)       # the numeric query
    expect_equal(res$nD, 40)        # the set
    ## a set as the query against a numeric database: still query first
    r2 <- testEnrichmentSEA(ids[1:40], list(score = z), platform = "MM285",
                            silent = TRUE)
    expect_equal(r2$nQ, 40)
    expect_equal(r2$nD, 500)
})
