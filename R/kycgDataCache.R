#' Cache KnowYourCG data
#'
#' @import ExperimentHub
#' @import AnnotationHub
#' @importFrom utils capture.output
#' @param data_titles data to cache, if not given will cache all
#' @return TRUE
#' @examples
#' kycgDataCache("KYCG.HM27.Mask.20220123")
#' ## to cache all data: kycgDataCache()
#' @export
kycgDataCache <- function(data_titles = NULL) {
    setExperimentHubOption(arg="MAX_DOWNLOADS", 100)
    if (is.null(data_titles)) {
        eh_ids <- unique(df_master$EHID) # cache all the records
    } else {
        stopifnot(all(data_titles %in% df_master$Title))
        eh_ids <- df_master$EHID[match(data_titles, df_master$Title)]
    }
    eh_ids <- eh_ids[!is.na(eh_ids)]
    if (length(eh_ids) == 0) return(invisible(TRUE));
    suppressMessages(try({
        eh_ids <- eh_ids[!(eh_ids %in% names(ExperimentHub(localHub=TRUE)))]
    }, silent = TRUE))
    if (length(eh_ids) == 0) return(invisible(TRUE));
    tryCatch({
        ## load meta data
        message(sprintf("Metadata (N=%d):\n", length(eh_ids)))
        suppressMessages(log <- capture.output(
            eh <- ExperimentHub()[eh_ids]))
        
        ## load actual data
        tmp2 <- lapply(seq_along(eh), function(i) {
            message(sprintf(
                "(%d/%d) %s:\n", i, length(eh_ids), eh_ids[i]))
            suppressMessages(log <- capture.output(cache(eh[i])))
        })
    }, error = function(cond) {
        message("ExperimentHub Caching fails:")
        message(cond)
        return(invisible(FALSE))
    })
    invisible(TRUE)
}
