
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

testCategoricalFisher <- function(database, ProbeIDs, sigProbes) {
    categories <- unique(database)
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
            names(OR) <- NULL
            return(OR)
        }
    )
    return(result)
}

#' test all databaseSet and return a list ranked by enrichment (odds-ratio)
testEnrichmentAll = function(probeIDs, pVals, databaseSets = NULL, percTop = 25, sig.threshold = 1e-8) {
    # assume same index for results and databases
    # get significant probes. probes and pVals must be in same order
    sigProbes <- probeIDs[pvals <= sig.threshold]
    # calculate number of databases to keep
    nDatabases <- ceiling((percTop/100)*length(databaseSets))
	# prioritize database sets by proportion of probes in dataset that overlap with significant hits
    # then take top percTop % database sets
    # note: databases read in from tbk files are named vectors
    topDatabases <- databases[
        order(
            lapply(
                names(databases),
                function(dbName) sum(sigProbes %in% databases[[dbName]]$ProbeID)
            ),
            decreasing = TRUE
        )
    ][1:nDatabases]
    
    # apply enrichment tests and get test results
    # TODO: need to name database list with database names? Otherwise can't know which is which
    results <- unlist(lapply(topDatabases, testCategoricalFisher, ProbeIDs = ProbeIDs, sigProbes = sigProbes))
    
    # sort results by output
    results <- sort(results, decreasing = TRUE)
    
    # TODO: handle continuous value databases
}


