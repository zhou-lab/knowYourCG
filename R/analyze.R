
#' test whether query is significantly enriched in database
testEnrichment = function(setQ, setD, setU) {
    mtx = matrix(c(
        length(intersect(setQ, setD)),
        length(setdiff(setD, setQ)), 
        length(setdiff(setQ, setD)), 
        length(setdiff(setU, union(setD, setQ)))),
        nrow = 2, 
        dimnames = list(
            Query = c("Q_in","Q_out"),
            Database = c("D_in","D_out")))
    list(mtx = mtx, test = fisher.test(mtx))
}
# ****how should the test be set up? test for each category (I/II, each CpG island?) something else?
testCategoricalFisher <- function(database, ProbeIDs, sigProbes) {
    categories <- unique(database)
    categories <- categories[!is.na(categories)]
    result <- sapply(
        categories,
        function(category) {
            setD <- probeIDs[database == category]
            fishertest <- testEnrichment(
                setQ = sigProbes,
                setD = setD,
                setU = probeIDs
            )
            OR <- c(fishertest$test$estimate)
            names(OR) <- category
            return(OR)
        }
    )
    return(result)
}

#' test all databaseSet and return a list ranked by enrichment (odds-ratio)
testEnrichmentAll = function(probeIDs, pVals, databaseSets = NULL, percTop = 25, sig.threshold = 1e-6) {
    # assume same index for results and databases
    # get significant probes. probes and pVals must be in same order
    sigProbes <- probeIDs[pVals <= sig.threshold]
    # if databaseSets is a single vector and not a matrix, turn it into a matrix
    if (is.null(dim(databaseSets))) databaseSets <- matrix(databaseSets, ncol = 1, dimnames = list(probeIDs))
    # calculate number of databases to keep
    nDatabases <- ceiling((percTop/100)*ncol(databaseSets))
	# prioritize database sets by proportion of probes in dataset that have a feature
    # and that overlap with significant hits
    # then take top percTop % database sets
    topDatabaseSelection <- order(
        apply(
            databaseSets, 
            2, 
            function(database) {
              length(intersect(sigProbes, probeIDs[!is.na(database)]))
            }
            ),
        decreasing = TRUE
    )[1:nDatabases]
    
    topDatabases <- databaseSets[, topDatabaseSelection]
    if (is.null(dim(topDatabases))) topDatabases <- matrix(topDatabases, ncol = 1, dimnames = list(probeIDs))
    
    # apply enrichment tests and get test results
    # TODO: need to name database list with database names? Otherwise can't know which is which
    results <- unlist(apply(topDatabases, 2, testCategoricalFisher, ProbeIDs = ProbeIDs, sigProbes = sigProbes))
    
    # sort results by output
    results <- sort(results, decreasing = TRUE)
    
    # TODO: handle continuous value databases
}


