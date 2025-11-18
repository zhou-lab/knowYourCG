#' Get databases by full or partial names of the database group(s)
#'
#' @param group_nms database group names
#' @param db_names name of the database, fetech only the given databases
#' @param platform EPIC, HM450, MM285, ... If given, will restrict to
#' that platform.
#' @param summary return a summary of database instead of db itself
#' @param allow_multi allow multiple groups to be returned for
#' @param type numerical, categorical, default: all
#' @param silent no messages
#' each query.
#' @return a list of databases, return NULL if no database is found
#' @examples
#' kycgDataCache(data_titles=
#' c("KYCG.MM285.chromHMM.20210210","KYCG.MM285.probeType.20210630"))
#' dbs <- getDBs("MM285.chromHMM")
#' dbs <- getDBs(c("MM285.chromHMM", "MM285.probeType"))
#' @export
getDBs <- function(
        group_nms, db_names = NULL, platform = NULL,
        summary = FALSE, allow_multi = FALSE, type = NULL, silent = FALSE) {

    if (!is.character(group_nms)) {
        return(group_nms)
    }

    group_nms <- guess_dbnames(group_nms, platform = platform,
                               allow_multi = TRUE, type = type,
                               silent = silent)
    group_nms <- group_nms[group_nms %in% df_master$Title]
    if (length(group_nms) == 0) {
        return(NULL)
    }
    res <- do.call(c, lapply(unname(group_nms), function(nm) {
        dbs <- kycgDataGet(nm)
        setNames(lapply(seq_along(dbs), function(ii) {
            db <- dbs[[ii]]
            attr(db, "group") <- nm
            attr(db, "dbname") <- names(dbs)[ii]
            db
        }), names(dbs))}))

    if (summary) {
        do.call(bind_rows, lapply(res, attributes))
    } else if (is.null(db_names)) {
        res
    } else {
        stopifnot(all(db_names %in% names(res)))
        res[db_names]
    }
}
