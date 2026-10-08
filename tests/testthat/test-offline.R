## Offline tests of the pure-R statistics that need no ExperimentHub data.

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

make_universe <- function(n = 5000) sprintf("cg%08d", seq_len(n))

test_that("testEnrichment runs on a user-supplied universe and databases", {
    u <- make_universe()
    query <- u[1:200]
    dbs <- as_dbs(list(hit = u[1:150], miss = u[1001:1150],
                       half = u[101:250]))
    res <- testEnrichment(query, dbs, universe = u, platform = "MM285",
                          silent = TRUE)
    expect_s3_class(res, "data.frame")
    expect_setequal(res$dbname, names(dbs))
    expect_lt(res$p.value[res$dbname == "hit"], 1e-10)
    expect_gt(res$estimate[res$dbname == "hit"], 0)
    expect_equal(res$overlap[res$dbname == "miss"], 0)
    expect_true(all(c("FDR", "group") %in% colnames(res)))
    expect_true(all(diff(res$log10.p.value) >= 0))       # ordered
    r2 <- testEnrichment(query, dbs, universe = u, platform = "MM285",
                         silent = TRUE, alternative = "two.sided",
                         mtc_by_group = FALSE)
    expect_true(all(r2$FDR >= r2$p.value))
    expect_error(testEnrichment(query, dbs, universe = u, platform = "MM285",
                                alternative = "bogus"))
})

test_that("Fisher helpers agree with fisher.test", {
    nD <- 150; nQ <- 200; nDQ <- 120; nU <- 5000
    res <- testEnrichmentFisherN(nD, nQ, nDQ, nU, alternative = "greater")
    ft <- fisher.test(matrix(c(nDQ, nQ - nDQ, nD - nDQ,
                               nU - nD - nQ + nDQ), 2), alternative = "greater")
    expect_equal(res$p.value, ft$p.value, tolerance = 1e-8)
    expect_true(all(c("estimate", "log10.p.value", "test", "nQ", "nD",
                      "overlap") %in% colnames(res)))
    ## capped log odds ratio when the overlap is complete
    r0 <- testEnrichmentFisherN(nD, nD, nD, nU, alternative = "greater")
    expect_true(is.finite(r0$estimate))
})

test_that("testEnrichmentSpearman matches cor.test and handles no overlap", {
    set.seed(1)
    ids <- make_universe(300)
    q <- setNames(rnorm(300), ids)
    db <- setNames(q + rnorm(300, sd = 0.5), ids)
    res <- testEnrichmentSpearman(q, db)
    ct <- cor.test(q, db, method = "spearman")
    expect_equal(res$estimate, ct$estimate[[1]])
    expect_equal(res$p.value, ct$p.value)
    expect_equal(res$overlap, 300)
    none <- testEnrichmentSpearman(q, setNames(db, paste0("x", ids)))
    expect_equal(none$estimate, 0)
    expect_equal(none$p.value, 1)
    expect_equal(none$overlap, 0)
})

test_that("databases_getMeta collects attributes and tolerates none", {
    a <- c("cg1", "cg2"); attr(a, "group") <- "G"; attr(a, "dbname") <- "a"
    b <- c("cg3"); attr(b, "group") <- "G"; attr(b, "dbname") <- "b"
    m <- databases_getMeta(list(a, b))
    expect_equal(m$dbname, c("a", "b"))
    expect_equal(m$group, c("G", "G"))
    expect_false("hasMeta" %in% colnames(m))
})

test_that("queryCheckPlatform returns the given platform or infers it", {
    expect_equal(queryCheckPlatform("EPIC", NULL), "EPIC")
    expect_error(queryCheckPlatform(NULL, NULL))
})

test_that("set_FDR adjusts within groups and errors without one", {
    res <- data.frame(p.value = c(0.01, 0.02, 0.5, 0.6),
                      group = c("A", "A", "B", "B"))
    out <- set_FDR(res)
    expect_equal(out$FDR[1:2], p.adjust(c(0.01, 0.02), "fdr"))
    expect_equal(out$FDR[3:4], p.adjust(c(0.5, 0.6), "fdr"))
    expect_error(set_FDR(data.frame(p.value = 0.1)), "no 'group'")
})

test_that("aggregateTestEnrichments builds a matrix across results", {
    r1 <- data.frame(dbname = c("a", "b"), estimate = c(1, 2))
    r2 <- data.frame(dbname = c("a", "b"), estimate = c(3, 4))
    m <- aggregateTestEnrichments(list(s1 = r1, s2 = r2))
    expect_equal(dim(m), c(2, 2))
    expect_equal(m["s2", "b"], 4)
    df <- aggregateTestEnrichments(list(s1 = r1, s2 = r2), return_df = TRUE)
    expect_true(is.data.frame(df))
})

test_that("plot helpers accept a minimal enrichment result", {
    set.seed(1)
    df <- data.frame(dbname = paste0("db", 1:30),
                     estimate = rnorm(30), p.value = runif(30) / 10,
                     group = rep(c("G1", "G2"), 15), overlap = 5:34,
                     nD = 100, nQ = 50, test = "Log2(OR)")
    df$log10.p.value <- log10(df$p.value)
    df$FDR <- p.adjust(df$p.value, "fdr")
    expect_s3_class(KYCG_plotBar(df), "ggplot")
    expect_s3_class(KYCG_plotDot(df), "ggplot")
    expect_s3_class(KYCG_plotVolcano(df), "ggplot")
    expect_s3_class(KYCG_plotLollipop(df), "ggplot")
    expect_false(is.null(KYCG_plotWaterfall(df)))
})
