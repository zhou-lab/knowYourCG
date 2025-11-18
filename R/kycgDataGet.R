valid_url <- function(url_in,t=2){
    con <- url(url_in)
    check <- suppressWarnings(try(open.connection(
        con,open="rt",timeout=t),silent=TRUE)[1])
    suppressWarnings(try(close.connection(con),silent=TRUE))
    ifelse(is.null(check),TRUE,FALSE)
}

## access zenodo
download_zenodo <- function(title, u1) {
    if (valid_url(u1)) {
        assign(title, get(load(url(u1))), envir = cacheEnv)
        TRUE
    } else {
        warning(sprintf("Resource %s cannot be retrieved.", title))
        FALSE
    }
    TRUE
}

stopAndCache <- function(title) {
    stop(sprintf('
| File %s either not found or needs to be cached to be
| used in KnowYourCG.
| Please make sure you have updated ExperimentHub and try
| > kycgDataCache("%s")
| or download all data
| > kycgDataCache()
| to retrieve and cache needed knowYourCG data.', title, title))
}

.kycgDataGet <- function(title) {
    idx <- match(title, df_master$Title)
    eh_id <- df_master$EHID[idx]
    stopifnot(length(eh_id) == 1)
    pfx <- df_master$Location_Prefix[idx]
    rdata <- df_master$RDataPath[idx]
    if (is.na(eh_id) && !is.na(pfx)) { # get data from zenodo if no EHID
        ## if title (not ehid) is in the cache
        if (exists(title, envir=cacheEnv, inherits=FALSE)) {
            return(get(title, envir=cacheEnv, inherits=FALSE))
        }

        if (download_zenodo(title, paste0(pfx, "/", rdata))) {
            return(get(title, envir=cacheEnv, inherits=FALSE))
        } else {
            stopAndCache(title)
        }
    } else { # get data from ExperimentHub if EHID
        if (exists(eh_id, envir=cacheEnv, inherits=FALSE)) {
            return(get(eh_id, envir=cacheEnv, inherits=FALSE))
        }
        if (!file.exists(getExperimentHubOption("CACHE"))) {
            stopAndCache(title) }
        tryCatch({
            eh <- ExperimentHub(localHub=TRUE)
        }, error = function(cond) { stopAndCache(title); })
        if (!(eh_id %in% names(eh))) {
            stopAndCache(title); }
        assign(eh_id, eh[[eh_id]], envir = cacheEnv)
        return(get(eh_id, envir=cacheEnv, inherits=FALSE))
    }
}

#' Get KnowYourCG data
#'
#' @param title title of the data
#' @param verbose whether to output ExperimentHub message
#' @return data object
#' @import ExperimentHub
#' @import AnnotationHub
#' @importFrom stringr str_replace
#' @examples
#'
#' kycgDataCache("KYCG.MSA.CGI.20220904")
#' EPIC.1.SigDF <- kycgDataGet('KYCG.MSA.CGI.20220904')
#' @export
kycgDataGet <- function(title, verbose = FALSE) {

    title <- str_replace(title, "MMB", "MM285") # fix potential code discrepancy
    if (verbose) {
        .kycgDataGet(title)
    } else {
        suppressMessages(
            log <- capture.output(obj <- .kycgDataGet(title)))
        obj
    }
}

