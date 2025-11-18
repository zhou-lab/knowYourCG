#' dbStats aggregates methylation of a given betas matrix over specified
#' database set features 
#'
#' @param betas matrix of beta values where probes are on the rows and
#' samples are on the columns
#' @param databases List of vectors corresponding to probe locations for
#' which the features will be extracted
#' @param fun aggregation function, default to mean
#' @param na.rm whether to remove NA
#' @param f_min min fraction of non-NA for aggregation function to apply
#' @param n_min min number of non-NA for aggregation function to apply,
#' overrides f_min
#' @param long produce long-form result
#' @return matrix with samples on the rows and database set on the columns
#' @import methods
#' @examples
#' library(SummarizedExperiment)
#' sesameData::sesameDataCache(data_titles=
#' c("MM285.467.SE.tissue20Kprobes","KYCG.MM285.probeType.20210630"))
#' se <- sesameData::sesameDataGet("MM285.467.SE.tissue20Kprobes")
#' head(dbStats(assay(se), "MM285.probeType")[,1:3])
#' sesameData::sesameDataGet_resetEnv()
#'
#' @export
dbStats <- function(
    betas, databases, fun = mean, na.rm = TRUE, n_min = NULL,
    f_min = 0.1, long = FALSE) {

    if (is(betas, "numeric")) { betas <- cbind(sample = betas); }
    if (is.character(databases)) {
        dbs <- getDBs(databases)
    } else {
        dbs <- databases
    }
    stats <- do.call(cbind, lapply(dbs, function(db) {
        betas1 <- betas[db[db %in% rownames(betas)],,drop=FALSE]
        n_probes <- nrow(betas1)
        if (n_probes == 0) { return(rep(NA, ncol(betas))); }
        nacnt <- colSums(!is.na(betas1), na.rm = TRUE)
        stat1 <- apply(betas1, 2, fun, na.rm = na.rm)
        if(is.null(n_min)) {
            n_min1 <- n_probes * f_min
        } else {
            n_min1 <- n_min
        }
        stat1[nacnt < n_min1] <- NA
        stat1
    }))
    if (!is.null(names(dbs))) {
        colnames(stats) <- names(dbs)
    } else { # use attributes instead of names
        colnames(stats) <- vapply(dbs,
            function(x) attr(x, "dbname"), character(1))
    }
    rownames(stats) <- colnames(betas)
    if (long) {
        stats <- melt(stats, varnames = c("query", "db"), value.name = "value")
    }
    stats
}
