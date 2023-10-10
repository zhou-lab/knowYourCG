test_that("testEnrichment works", {
    library(SummarizedExperiment)
    df <- rowData(sesameDataGet('MM285.tissueSignature'))
    query <- df$Probe_ID[df$branch == "B_cell"]
    res <- testEnrichment(query, "chromHMM", platform="MM285")
    expect_type(res, "list")
})
