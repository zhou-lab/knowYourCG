
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

calcFoldChange = function(mtx){
	num = mtx[1, 1] / (mtx[1, 1] + mtx[1, 2])
	den = (mtx[1, 1] + mtx[2, 1]) / sum(mtx)
	num / den
}

#' test all databaseSet and return a list ranked by enrichment (odds-ratio)
testEnrichmentAll = function() {
	
}


