## The functions that need ExperimentHub data: one small call each, skipped
## without a network connection. Every record they use is MM285.

mm_probes <- function(n = 50) {
    df <- SummarizedExperiment::rowData(
        sesameData::sesameDataGet("MM285.tissueSignature"))
    head(df$Probe_ID, n)
}

test_that("kycgDataCache caches a record by title", {
    skip_if_offline()
    expect_true(kycgDataCache("KYCG.MM285.chromHMM.20210210"))
    expect_true(kycgDataCache(character(0)))
    expect_error(kycgDataCache("not_a_title"))
})

test_that("buildGeneDBs builds gene sets around a few probes", {
    skip_if_offline()
    dbs <- suppressMessages(buildGeneDBs(mm_probes(20), platform = "MM285",
                                         silent = TRUE))
    expect_type(dbs, "list")
    expect_gt(length(dbs), 0)
    expect_true(grepl("gene", attr(dbs[[1]], "group")))
    expect_false(is.null(attr(dbs[[1]], "gene_name")))
})

test_that("linkProbesToProximalGenes annotates probes with genes", {
    skip_if_offline()
    gr <- suppressMessages(linkProbesToProximalGenes(mm_probes(20),
                                                     platform = "MM285"))
    expect_s4_class(gr, "GRanges")
    expect_equal(convertGeneName("FOXA2"), "Foxa2")
})

test_that("KYCG_plotMeta and KYCG_plotManhattan draw from the manifest", {
    skip_if_offline()
    ids <- mm_probes(2000)
    set.seed(1)
    betas <- setNames(runif(length(ids)), ids)
    expect_s3_class(suppressMessages(KYCG_plotMeta(betas, platform = "MM285")),
                    "ggplot")
    vals <- setNames(rexp(length(ids)), ids)
    expect_s3_class(suppressMessages(KYCG_plotManhattan(
        vals, platform = "MM285", label_min = 1e9)), "ggplot")
})
