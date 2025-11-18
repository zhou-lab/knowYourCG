#' testEnrichmentSpearman uses the Spearman statistical test to estimate
#' the association between two continuous variables.
#'
#' @param num_query named numeric vector of probes of interest 
#' where names are probe IDs (e.g significant probes)
#' @param num_db List of vectors corresponding to the database set of
#' interest with associated meta data as an attribute to each element.
#'
#' @import stats
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name
#' of test for the given results.
testEnrichmentSpearman <- function(num_query, num_db) {
    test <- "Spearman's rho"
    if (length(intersect(names(num_query), names(num_db))) == 0) {
        return(data.frame(
            estimate = 0, p.value = 1,
            log10.p.value = 0, test = test,
            nQ = length(num_query), nD = length(num_db),
            overlap = 0))
    }

    num_db <- num_db[match(names(num_query), names(num_db))]

    res <- cor.test(num_query, num_db, method = "spearman")
    data.frame(
        estimate = res$estimate[[1]], p.value = res$p.value,
        log10.p.value = log10(res$p.value), test = test,
        nQ = length(num_query), nD = length(num_db),
        overlap = length(num_query))
}
