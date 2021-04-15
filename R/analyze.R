
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

#' test all databaseSet and return a list ranked by enrichment (odds-ratio)
testEnrichmentAll = function(experiment, databaseSets = NULL, percTop = 25, sig.threshold = 1e-8) {
	# prioritize database sets by proportion of probes in dataset that overlap with significant hits
    
    # take top percTop % database sets
    
    # apply enrichment tests and get scores
    
    # continuous variable stuff
}


