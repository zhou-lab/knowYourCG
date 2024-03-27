#' testEnrichment tests for the enrichment of a set of probes (query set) in
#' a number of features (database sets).
#'
#' @param probeIDs Vector of probes of interest (e.g., significant probes)
#' @param databases List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element.
#' Optional. (Default: NA)
#' @param universe Vector of probes in the universe set containing all of
#' the probes to be considered in the test. If it is not provided, it will be
#' inferred from the provided platform. (Default: NA).
#' @param alternative "two.sided", "greater", or "less"
#' @param include_genes include gene link enrichment testing
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @param silent output message? (Default: FALSE)
#' @return A data frame containing features corresponding to the test estimate,
#' p-value, and type of test.
#' @importFrom dplyr bind_rows
#' @examples
#'
#' library(SummarizedExperiment)
#' sesameData::sesameDataCache(data_titles=
#' c("KYCG.MM285.tissueSignature.20211211","KYCG.MM285.chromHMM.20210210"))
#' df <- rowData(sesameData::sesameDataGet('MM285.tissueSignature'))
#' probes <- df$Probe_ID[df$branch == "B_cell"]
#' res <- testEnrichment(probes, "chromHMM", platform="MM285")
#' sesameData::sesameDataGet_resetEnv()
#'
#' @export
testEnrichment <- function(
        probeIDs, databases = NULL, universe = NULL, alternative = "greater",
        include_genes = FALSE, platform = NULL, silent = FALSE) {

    platform <- queryCheckPlatform(platform, probeIDs, silent = silent)

    if (is.null(databases)) {
        dbs <- c(getDBs(listDBGroups( # by default, all dbs + gene
            platform, type="categorical")$Title, silent = silent))
    } else if (is.character(databases)) {
        dbs <- getDBs(databases, platform = platform, silent = silent)
    } else {
        dbs <- databases
    }

    if (include_genes) {
        dbs <- c(dbs, buildGeneDBs(probeIDs, platform, silent = silent))
    }

    ## there shouldn't be empty databases, but just in case
    dbs <- dbs[vapply(dbs, length, integer(1)) > 0]
    if (!silent) {
        message(sprintf("Testing against %d database(s)...", length(dbs)))
    }

    if (is.null(universe)) {
        universe <- inferUniverse(platform)
    } else { # subset the dbs by universe
        dbs <- subsetDBs(dbs, universe) }

    res <- do.call(bind_rows, lapply(dbs, function(db) {
        testEnrichmentFisher(query = probeIDs, database = db,
                             universe = universe, alternative = alternative)}))

    ## adjust p.value after merging
    res$FDR <- p.adjust(res$p.value, method='fdr')
    rownames(res) <- NULL

    ## bind meta data
    res <- cbind(res, databases_getMeta(dbs))
    res[order(res$log10.p.value, -abs(res$estimate)), ]
}



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
#' sesameData::sesameDataCache(data_titles=
#' c("KYCG.MM285.TFBSconsensus.20220116","KYCG.MM285.chromHMM.20210210"))
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

#' testEnrichmentFisher uses Fisher's exact test to estimate the association
#' between two categorical variables.
#'
#' Estimates log2 Odds ratio
#'
#' @param query Vector of probes of interest (e.g., significant probes)
#' @param database Vectors corresponding to the database set of
#' interest with associated meta data as an attribute to each element.
#' @param universe Vector of probes in the universe set containing all of
#' @param alternative greater or two.sided (default: greater)
#' the probes to be considered in the test. (Default: NULL)
#'
#' @import stats
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
testEnrichmentFisher <- function(query, database, universe,
                                 alternative = "greater") {

    nD <- length(database)
    nQ <- length(query)
    nDQ <- length(intersect(query, database))
    nU <- length(universe)

    testEnrichmentFisherN(nD, nQ, nDQ, nU, alternative = alternative)
}


testEnrichmentFisherN <- function(
        nD, nQ, nDQ, nU, alternative = "greater") {

    nDmQ <- nD - nDQ
    nQmD <- nQ - nDQ
    nUmDQ <- nU - nQ - nD + nDQ

    if (alternative == "two.sided") {
        pvg <- phyper(
            nDQ-1, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
            lower.tail = FALSE, log.p = TRUE) / log(10)
        pvl <- phyper(
            nDQ, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
            lower.tail = TRUE, log.p = TRUE) / log(10)
        log10.p.value <- pmin(pmin(pvg, pvl) + log(2), 0) / log(10)
        ## log10.p.value <- log10(fisher.test(matrix(c(
        ##     nDQ, nDmQ, nQmD, nUmDQ), nrow = 2))$p.value)
    } else if (alternative == "greater") {
        log10.p.value <- phyper(
            nDQ-1, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
            lower.tail = FALSE, log.p = TRUE) / log(10)
    } else if (alternative == "less") {
        log10.p.value <- phyper(
            nDQ, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
            lower.tail = TRUE, log.p = TRUE) / log(10)
    } else {
        stop("alternative must be either greater, less or two-sided.")
    }

    odds_ratio <- nDQ / nQmD / nDmQ * nUmDQ # can be NaN if 0
    odds_ratio[odds_ratio == Inf] <- .Machine$double.xmax
    odds_ratio[odds_ratio == 0] <- .Machine$double.xmin
    data.frame(
        estimate = log2(odds_ratio),
        p.value = 10**(log10.p.value),
        log10.p.value = log10.p.value,
        test = "Log2(OR)",
        nQ = nQ, nD = nD, overlap = nDQ,
        cf_Jaccard = nDQ / (nD + nQmD),
        cf_overlap = nDQ / pmin(nD, nQ), # Szymkiewicz–Simpson
        cf_NPMI = (log2(nD)+log2(nQ)-2*log2(nU))/(log2(nDQ)-log2(nU))-1,
        cf_SorensenDice = 2 * nDQ/(nD + nQ))
}



#' build gene-probe association database
#'
#' @param probeIDs the query probe list. If NULL, use all the probes
#' on the platform
#' @param platform HM450, EPIC, MM285, Mammal40, will infer from
#' query if not given
#' @param genome hg38, mm10, ..., will infer if not given.
#' @param max_distance probe-gene distance for association
#' @param silent suppress messages
#' @return gene databases
#' @import sesameData
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @examples
#' sesameData::sesameDataCache(data_titles=
#' c("EPIC.address","genomeInfo.hg38","probeIDSignature"))
#' query <- c("cg04707299", "cg13380562", "cg00480749")
#' dbs <- buildGeneDBs(query, platform = "EPIC")
#' testEnrichment(query, dbs, platform = "EPIC")
#' @export
buildGeneDBs <- function(
    probeIDs = NULL, platform = NULL,
        genome = NULL, max_distance = 10000, silent = FALSE) {

    platform <- queryCheckPlatform(platform, probeIDs, silent = silent)
    genes <- sesameData_txnToGeneGRanges(
        sesameData_getTxnGRanges(
            sesameData_check_genome(NULL, platform)))
    all_probes <- sesameData_getManifestGRanges(platform, genome = genome)
    if (!is.null(probeIDs)) {
        probes <- all_probes[names(all_probes) %in% probeIDs] }

    ## skip non-overlapping genes, strand always ignored
    genes <- subsetByOverlaps(
        genes, probes + max_distance, ignore.strand = TRUE)
    hits <- findOverlaps(
        genes, all_probes + max_distance, ignore.strand = TRUE)
    dbs <- split(names(all_probes)[subjectHits(hits)],
                 names(genes)[queryHits(hits)])
    gene_names <- genes[names(dbs)]$gene_name
    res <- lapply(seq_along(dbs), function(i) {
        d1 <- dbs[[i]];
        attr(d1, "group") <- sprintf("KYCG.%s.gene.00000000", platform);
        attr(d1, "dbname") <- names(dbs)[i];
        attr(d1, "gene_name") <- gene_names[i];
        d1;})
    names(res) <- names(dbs)
    message(sprintf("Building %d gene DBs for %s...", length(res), platform))
    res
}
