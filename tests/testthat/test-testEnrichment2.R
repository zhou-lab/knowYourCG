## testEnrichment2: the sequencing path through the compiled YAME summary,
## run on the two files the package ships, plus its parsing helpers.

kb <- system.file("extdata", "chromhmm.cm", package = "knowYourCG")
qry <- system.file("extdata", "onecell.cg", package = "knowYourCG")

test_that("testEnrichment2 runs on the shipped .cg and .cm", {
    skip_on_os("windows")
    skip_if(!nzchar(kb) || !nzchar(qry), "extdata not installed")
    ## no warnings: integer counts once overflowed in the odds ratio
    expect_silent(res <- testEnrichment2(qry, kb))
    expect_s3_class(res, "tbl_df")
    expect_gt(nrow(res), 1)
    expect_true(all(c("Mask", "N_mask", "N_query", "N_overlap", "N_univ",
                      "estimate", "p.value", "log10.p.value", "MFile") %in%
                    colnames(res)))
    expect_true(all(diff(res$log10.p.value) >= 0))       # ordered
    expect_true(all(res$N_overlap <= res$N_mask))
    expect_true(all(is.finite(res$estimate)))

    ## the same counts under another alternative, and verbose messages
    expect_message(r2 <- testEnrichment2(qry, kb, alternative = "two.sided",
                                         verbose = TRUE), "Reading query")
    expect_setequal(r2$Mask, res$Mask)
    expect_error(testEnrichment2(qry, kb, alternative = "bogus"))

    ## a minimum overlap nothing meets leaves an empty result, one warning
    big <- max(res$N_overlap) + 1
    w <- capture_warnings(r3 <- testEnrichment2(qry, kb, min_overlap = big))
    expect_length(w, 1)
    expect_match(w, "minimum overlap")
    expect_equal(nrow(r3), 0)
})

test_that("testEnrichment2 validates its inputs", {
    expect_error(testEnrichment2(1, "x"), "single character string")
    expect_error(testEnrichment2(c("a", "b"), "x"), "single character string")
    expect_error(testEnrichment2(qry, 1), "character string or vector")
    expect_error(testEnrichment2(qry, kb, universe_fn = 1),
                 "NULL or a character string")
    expect_error(testEnrichment2("no_such.cg", kb), "Query file not found")
    expect_error(testEnrichment2(qry, "no_such.cm"),
                 "Knowledgebase file\\(s\\) not found")
    expect_error(testEnrichment2(qry, kb, universe_fn = "no_such.cg"),
                 "Universe file not found")
})

test_that("parse_yame_results reads YAME's table and rejects bad output", {
    out <- c("QFile\tQuery\tMFile\tMask\tN_univ\tN_query\tN_mask\tN_overlap",
             "q.cg\tq\tm.cm\tA\t100\t20\t30\t10",
             "q.cg\tq\tm.cm\tB\t100\t20\t40\t2")
    df <- parse_yame_results(out)
    expect_equal(nrow(df), 2)
    expect_equal(df$N_overlap, c(10, 2))
    expect_error(parse_yame_results(NULL), "Empty result")
    expect_error(parse_yame_results(character(0)), "Empty result")
    expect_error(parse_yame_results(c("Mask\tN_mask", "A\t3")),
                 "missing required columns")
    expect_error(parse_yame_results(c("a\tb", "1\t2\t3\t4")),
                 "Failed to parse")
})

test_that("compute_enrichment_stats adds Fisher statistics per mask", {
    df <- tibble::tibble(Mask = c("A", "B"), N_univ = c(100, 100),
                         N_query = c(20, 20), N_mask = c(30, 40),
                         N_overlap = c(10, 2))
    res <- compute_enrichment_stats(df, "greater", 1)
    expect_equal(nrow(res), 2)
    expect_true(all(c("estimate", "p.value", "MFile") %in% colnames(res)))
    expect_equal(res$MFile, c("YAME", "YAME"))
    expect_lt(res$p.value[1], res$p.value[2])
    res2 <- compute_enrichment_stats(df, "greater", 5)
    expect_equal(res2$Mask, "A")
    expect_warning(r0 <- compute_enrichment_stats(df, "greater", 50),
                   "minimum overlap threshold")
    expect_equal(nrow(r0), 0)
})
