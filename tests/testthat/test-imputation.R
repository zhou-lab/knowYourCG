test_that("imputation works", {
    se <- sesameDataGet('MM285.467.SE.tissue20Kprobes')
    betas <- assay(se)[1:5000,1:100]
    betas <- imputeMissingProbes(betas)
  expect_equal(sum(is.na(betas)), 0)
})
