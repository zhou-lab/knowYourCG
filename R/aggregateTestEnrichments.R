#' Aggregate test enrichment results
#'
#' @param result_list a list of results from testEnrichment
#' @param column the column name to aggregate (Default: estimate)
#' @param return_df whether to return a merged data frame
#' @return a matrix for all results
#' @importFrom reshape2 melt
#' @examples
#'
#' ## pick some big TFBS-overlapping CpG groups
#' kycgDataCache(data_titles=
#' c("KYCG.MM285.TFBSconsensus.20220116","KYCG.MM285.chromHMM.20210210"))
#'
#' sesameDataCache(data_titles=
#' "probeIDSignature", "MM285.address"))
#' 
#' cg_lists <- getDBs("MM285.TFBS")
#' queries <- cg_lists[(sapply(cg_lists, length) > 40000)]
#' result_list <- lapply(queries, testEnrichment, "MM285.chromHMM")
#' mtx <- aggregateTestEnrichments(result_list)
#'
#' @export
aggregateTestEnrichments <- function(
        result_list, column = "estimate", return_df = FALSE) {
    mtx <- do.call(cbind, lapply(result_list[[1]]$dbname, function(db) {
        vapply(result_list,
               function(x) x$estimate[x$dbname == db], numeric(1))}))
    colnames(mtx) <- result_list[[1]]$dbname
    if (return_df) {
        melt(mtx, value.name = column, varnames = c("query", "db"))
    } else {
        mtx
    }
}
