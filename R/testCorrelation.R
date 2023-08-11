#' testEnrichmentSpearman uses the Spearman statistical test to estimate
#' the association between two continuous variables.
#'
#' @param query Vector of probes of interest (e.g., significant probes)
#' @param database List of vectors corresponding to the database set of
#' interest with associated meta data as an attribute to each element.
#'
#' @import stats
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name
#' of test for the given results.
testEnrichmentSpearman <- function(query, database) {
    test <- "Spearman's rho"
    if (length(intersect(names(query), names(database))) == 0) {
        return(data.frame(
            estimate = 0, p.value = 1,
            log10.p.value = 0, test = test,
            nQ = length(query), nD = length(database),
            overlap = 0))
    }

    database <- database[match(names(query), names(database))]

    res <- cor.test(query, database, method = "spearman")
    data.frame(
        estimate = res$estimate[[1]], p.value = res$p.value,
        log10.p.value = log10(res$estimate[[1]]), test = test,
        nQ = length(query), nD = length(database),
        overlap = length(query))
}
