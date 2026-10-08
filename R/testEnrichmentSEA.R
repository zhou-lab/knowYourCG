calcES_Significance <- function(dCont, dDisc, permut=100, precise=FALSE) {

    dCont <- sort(dCont)
    dContName <- names(dCont)
    dDiscN <- length(dDisc)
    dContN <- length(dCont)
    s <- rep(-1/(dContN-dDiscN), dContN)

    ess <- do.call(rbind, lapply(seq_len(permut), function(i) {
        s[sample.int(dContN, dDiscN)] <- 1/dDiscN
        cs <- cumsum(s)
        data.frame(es_max = max(cs), es_min = min(cs))
    }))

    ## es <- calcES(dCont, dDisc)
    presence <- names(dCont) %in% dDisc
    s <- ifelse(presence, 1/sum(presence), -1/sum(!presence))
    cs <- cumsum(s)
    es_max <- max(cs)
    es_min <- min(cs)
    ## Under the null, max(cs) and -min(cs) have the same distribution
    ## (reverse and negate the walk), so pool them into one null for both
    ## tails. Separate nulls made the two tails disagree, and the p-value
    ## depended on the sign of the query (issue #6).
    es_null <- c(ess$es_max, -ess$es_min)
    res <- list(es_small = es_max,
                es_large = -es_min,
                pv_small = mean(es_null >= es_max),
                pv_large = mean(es_null >= -es_min))

    ## if significant, try to be more precise
    if (res$pv_small < 0.01 || res$pv_large < 0.01) {
        if (permut < 1000 && precise) {
            res <- calcES_Significance(dCont, dDisc, permut = 1000)
        } else { # approximated by Gaussian (TODO: also report log.p=TRUE)
            if (res$pv_small == 0) {
                res$pv_small <- pnorm(
                    es_max, mean=mean(es_null),
                    sd=sd(es_null), lower.tail=FALSE) }
            if (res$pv_large == 0) {
                res$pv_large <- pnorm(
                    -es_min, mean=mean(es_null),
                    sd=sd(es_null), lower.tail=FALSE)
            }}}

    res
}


testEnrichmentSEA1 <- function(query, database, precise=FALSE, full=FALSE) {

    test <- "Set Enrichment Score"
    overlap <- intersect(names(query), database)
    if (length(overlap) != length(database)) {
        warning("Not every data in database has query.")
        warning(sprintf("Using %d in %d data for testing.",
                        length(overlap), length(database)))
    }

    if (length(overlap) == 0 || length(overlap) == length(query)) {
        return(data.frame(
            estimate = 0, p.value = 1,
            log10.p.value = 0, test = test,
            nQ = length(database), nD = length(query),
            overlap = length(overlap))) }

    res <- calcES_Significance(query, overlap, precise=precise)

    if (res$es_large > res$es_small) {
        df <- data.frame( ## positive sign: enriched at large values, as GSEA
            estimate = res$es_large, p.value = res$pv_large,
            log10.p.value = log10(res$pv_large), test = test,
            nQ = length(database), nD = length(query),
            overlap = length(overlap))
    } else {
        df <- data.frame( ## negative sign: enriched at small values
            estimate = -res$es_small, p.value = res$pv_small,
            log10.p.value = log10(res$pv_small), test = test,
            nQ = length(database), nD = length(query),
            overlap = length(overlap))
    }

    if (full) {
        list(res = df, dCont = query, dDisc = overlap)
    } else {
        df
    }
}


#' uses the GSEA-like test to estimate the association of a
#' categorical variable against a continuous variable.
#'
#' estimate is the enrichment score. A positive estimate means the
#' categorical set is enriched at the large values of the continuous
#' variable; a negative estimate means it is enriched at the small values.
#' This matches GSEA. Before version 1.9.3 the sign was the reverse.
#'
#' @param query query, if numerical, expect categorical database, if
#' categorical expect numerical database
#' @param databases database, numerical or categorical, but needs to be
#' different from query
#' @param platform EPIC, MM285, ..., infer if not given
#' @param silent suppress message (default: FALSE)
#' @param precise whether to compute precise p-value (up to numerical limit)
#' of interest.
#' @param prepPlot return the raw enrichment scores and presence vectors
#' for plotting
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
#' @examples
#' sesameData::sesameDataCache(data_titles=
#' c("KYCG.MM285.designGroup.20210210","KYCG.MM285.seqContextN.20210630",
#' "probeIDSignature"))
#' query <- getDBs("KYCG.MM285.designGroup")[["TSS"]]
#' res <- testEnrichmentSEA(query, "MM285.seqContextN")
#' @export
testEnrichmentSEA <- function(query, databases,
                              platform = NULL, silent = FALSE,
                              precise = FALSE, prepPlot = FALSE) {

    platform <- queryCheckPlatform(platform, query, silent = silent)
    stopifnot(!is.null(databases))
    if (is.character(databases)) { # infer database from string
        dbs <- getDBs(databases, platform = platform, silent = silent)
    } else {
        dbs <- databases
    }
    dbs <- dbs[vapply(dbs, length, integer(1)) > 0] # no empty db
    if (!silent) {
        message(sprintf("Testing against %d database(s)...", length(dbs))) }

    if (is.character(query) && all(vapply(dbs, is.numeric, logical(1)))) {
        res <- lapply(dbs, function(db) {
            testEnrichmentSEA1(query = db, database = query,
                               precise = precise, full=prepPlot)})
    } else if (
        is.numeric(query) && all(vapply(dbs, is.character, logical(1)))) {
        res <- lapply(dbs, function(db) {
            testEnrichmentSEA1(query = query, database = db,
                               precise = precise, full=prepPlot)})
    } else { stop("query and db must be one numerical and one categorical"); }

    if (prepPlot) { return(res)
    } else { res <- do.call(bind_rows, res) }

    ## adjust p.value after merging
    res$FDR <- p.adjust(res$p.value, method='fdr')
    rownames(res) <- NULL

    ## bind meta data
    res <- cbind(res, databases_getMeta(dbs))
    res[order(res$p.value, -abs(res$estimate)), ]
}

