#' Test for enrichment of query in knowledgebase sets
#'
#' @param query For array input, a vector of probes of interest
#' (e.g., significant differential methylated probes). For sequencing data
#' input, the file name for YAME-compressed CG sets.
#' @param databases List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element.
#' If NULL, all available databases for the platform are used. (Default: NULL)
#' @param universe Vector of probes in the universe set containing all
#' probes to be considered in the test. If NULL, will be inferred from the
#' provided platform. (Default: NULL)
#' @param alternative Test alternative: "two.sided", "greater", or "less".
#' (Default: "greater")
#' @param include_genes Include gene link enrichment testing. (Default: FALSE)
#' @param platform String corresponding to the type of platform to use:
#' MM285, EPIC, HM450, or HM27. If NULL, will be inferred from query set
#' probe IDs. (Default: NULL)
#' @param silent Suppress output messages? (Default: FALSE)
#' @param mtc_by_group Perform multiple testing correction within 
#' knowledgebase groups. (Default: TRUE)
#' @param mtc_method Method for multiple test correction. (Default: "fdr")
#' @return A data frame containing features corresponding to the test estimate,
#' p-value, and type of test, ordered by significance.
#' @importFrom dplyr bind_rows
#' @examples
#'
#' library(SummarizedExperiment)
#' library(sesameData)
#' library(knowYourCG)
#' kycgDataCache(data_titles = "KYCG.MM285.chromHMM.20210210")
#' sesameDataCache("MM285.tissueSignature")
#' df <- rowData(sesameDataGet("MM285.tissueSignature"))
#' probes <- df$Probe_ID[df$branch == "B_cell"]
#' res <- testEnrichment(probes, "chromHMM", platform = "MM285")
#' 
#' \donttest{
#' # Define temporary directory and file URLs
#' temp_dir <- tempdir()
#' knowledgebase <- file.path(temp_dir, "ChromHMM.20220414.cm")
#' query <- file.path(temp_dir, "mm10_f3_10cells.cg")
#' 
#' # URLs for the knowledgebase and query files
#' knowledgebase_url <- paste0(
#'   "https://zenodo.org/records/18175656/files/",
#'   "ChromHMM.20220414.cm"
#' )
#' query_url <- paste0(
#'   "https://zenodo.org/records/18176004/files/",
#'   "mm10_f3_10cells.cg"
#' )
#' 
#' # Download the files
#' download.file(knowledgebase_url, destfile = knowledgebase)
#' download.file(query_url, destfile = query)
#' 
#' # Confirm file download
#' list.files(temp_dir)
#' res <- testEnrichment(query, knowledgebase)
#' }
#' @export
testEnrichment <- function(
    query, 
    databases = NULL, 
    universe = NULL, 
    alternative = "greater",
    include_genes = FALSE, 
    platform = NULL, 
    silent = FALSE,
    mtc_by_group = TRUE, 
    mtc_method = "fdr") {
    
    ## Validate alternative parameter
    alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
    
    ## Check if query is a file path for sequencing data
    is_seq_data <- length(query) == 1 && 
                   !grepl("^c[gh]", query) &&
                   !grepl("^rs", query) && 
                   is.null(platform)
    
    if (is_seq_data) {
        res <- testEnrichment2(
            query, 
            databases, 
            universe_fn = universe,
            alternative = alternative
        )
    } else {
        ## Array-based enrichment testing
        res <- perform_array_enrichment(
            query = query,
            databases = databases,
            universe = universe,
            alternative = alternative,
            include_genes = include_genes,
            platform = platform,
            silent = silent
        )
    }
    
    ## Apply multiple testing correction
    res <- set_FDR(res, mtc_by_group = mtc_by_group, mtc_method = mtc_method)
    
    ## Order results by significance
    res[order(res$log10.p.value, -abs(res$estimate)), ]
}

## Helper function to perform array-based enrichment testing
perform_array_enrichment <- function(
    query, databases, universe, alternative, 
    include_genes, platform, silent) {
    
    ## Validate platform
    platform <- queryCheckPlatform(platform, query, silent = silent)
    
    ## Get databases
    dbs <- get_databases(databases, platform, silent)
    
    ## Add gene databases if requested
    if (include_genes) {
        gene_dbs <- buildGeneDBs(query, platform, silent = silent)
        dbs <- c(dbs, gene_dbs)
    }
    
    ## Remove empty databases
    dbs <- Filter(function(x) length(x) > 0, dbs)
    
    if (!silent) {
        message(sprintf("Testing against %d database(s)...", length(dbs)))
    }
    
    ## Get or validate universe
    universe <- get_universe(universe, platform, dbs)
    
    ## Perform enrichment tests
    res_list <- lapply(dbs, function(db) {
        nD <- length(db)
        nQ <- length(query)
        nDQ <- length(intersect(query, db))
        nU <- length(universe)
        
        testEnrichmentFisherN(nD, nQ, nDQ, nU, alternative = alternative)
    })
    
    ## Combine results
    res <- dplyr::bind_rows(res_list)
    rownames(res) <- NULL
    
    ## Add metadata
    cbind(res, databases_getMeta(dbs))
}

## Helper function to get databases
get_databases <- function(databases, platform, silent) {
    if (is.null(databases)) {
        ## Get all available databases by default
        db_groups <- listDBGroups(platform)$Title
        getDBs(db_groups, silent = silent, type = "categorical")
    } else if (is.character(databases)) {
        getDBs(databases, platform = platform, silent = silent)
    } else {
        databases
    }
}

## Helper function to get or validate universe
get_universe <- function(universe, platform, dbs) {
    if (is.null(universe)) {
        address_data <- sesameDataGet(paste0(platform, ".address"))
        address_data$ordering$Probe_ID
    } else {
        ## Subset databases by universe
        subsetDBs(dbs, universe)
        universe
    }
}

## Apply FDR correction to enrichment results
set_FDR <- function(res, mtc_by_group = TRUE, mtc_method = "fdr") {
    
    if (!mtc_by_group) {
        res$FDR <- p.adjust(res$p.value, method = mtc_method)
        return(res)
    }
    
    ## Determine grouping variable
    group <- determine_group(res)
    
    ## Calculate FDR by group
    grp_ind <- split(seq_len(nrow(res)), group)
    grp_fdr <- lapply(grp_ind, function(idx) {
        p.adjust(res$p.value[idx], method = mtc_method)
    })
    
    ## Assign FDR values back to results
    res$FDR[unlist(grp_ind)] <- unname(unlist(grp_fdr))
    res
}

## Determine grouping variable for FDR correction
determine_group <- function(res) {
    if (!is.null(res$group)) {
        ## Array data
        res$group
    } else if (!is.null(res$MFile)) {
        ## Sequencing data
        res$MFile
    } else {
        stop(
            "Cannot adjust p-values by group: ",
            "no 'group' or 'MFile' column found.")
    }
}

## Calculate Fisher's exact test from counts
testEnrichmentFisherN <- function(
    nD, nQ, nDQ, nU, alternative = "greater") {
    
    ## Calculate contingency table values
    nDmQ <- nD - nDQ      ## In database but not query
    nQmD <- nQ - nDQ      ## In query but not database
    nUmDQ <- nU - nQ - nD + nDQ  ## In neither
    
    ## Calculate p-value based on alternative hypothesis
    log10.p.value <- calculate_fisher_pvalue(
        nDQ, nQmD, nUmDQ, nDmQ, alternative
    )
    
    ## Calculate odds ratio with safeguards
    odds_ratio <- calculate_odds_ratio(nDQ, nQmD, nDmQ, nUmDQ)
    
    ## Calculate additional effect size metrics
    cf_metrics <- calculate_effect_sizes(nD, nQ, nDQ, nU, nQmD, nDmQ, nUmDQ)
    
    data.frame(
        estimate = log2(odds_ratio),
        p.value = 10^log10.p.value,
        log10.p.value = log10.p.value,
        test = "Log2(OR)",
        nU = nU, 
        nQ = nQ, 
        nD = nD, 
        overlap = nDQ,
        cf_Jaccard = cf_metrics$jaccard,
        cf_MCC = cf_metrics$mcc,
        cf_overlap = cf_metrics$overlap,
        cf_NPMI = cf_metrics$npmi,
        cf_SorensenDice = cf_metrics$sorensen_dice
    )
}

## Calculate Fisher's exact test p-value
calculate_fisher_pvalue <- function(nDQ, nQmD, nUmDQ, nDmQ, alternative) {
    m <- nDQ + nQmD
    n <- nUmDQ + nDmQ
    k <- nDmQ + nDQ
    
    if (alternative == "two.sided") {
        pvg <- phyper(nDQ - 1, m, n, k, lower.tail = FALSE,
            log.p = TRUE) / log(10)
        pvl <- phyper(nDQ, m, n, k, lower.tail = TRUE,
            log.p = TRUE) / log(10)
        pmin(pmin(pvg, pvl) + log10(2), 0)
    } else if (alternative == "greater") {
        phyper(nDQ - 1, m, n, k, lower.tail = FALSE, log.p = TRUE) / log(10)
    } else if (alternative == "less") {
        phyper(nDQ, m, n, k, lower.tail = TRUE, log.p = TRUE) / log(10)
    } else {
        stop("alternative must be 'greater', 'less', or 'two.sided'.")
    }
}

## Calculate odds ratio with safeguards for extreme values
calculate_odds_ratio <- function(nDQ, nQmD, nDmQ, nUmDQ) {
    odds_ratio <- (nDQ * nUmDQ) / (nQmD * nDmQ)
    
    ## Handle extreme values
    odds_ratio[is.infinite(odds_ratio)] <- .Machine$double.xmax
    odds_ratio[odds_ratio == 0] <- .Machine$double.xmin
    odds_ratio[is.nan(odds_ratio)] <- NA_real_
    
    odds_ratio
}

## Calculate effect size metrics
calculate_effect_sizes <- function(nD, nQ, nDQ, nU, nQmD, nDmQ, nUmDQ) {
    list(
        jaccard = nDQ / (nD + nQmD),
        mcc = calculate_mcc(nDQ, nUmDQ, nQmD, nDmQ, nD, nU, nQ),
        overlap = nDQ / pmin(nD, nQ),  ## Szymkiewicz-Simpson coefficient
        npmi = calculate_npmi(nD, nQ, nDQ, nU),
        sorensen_dice = (2 * nDQ) / (nD + nQ)
    )
}

## Calculate Matthews correlation coefficient
calculate_mcc <- function(nDQ, nUmDQ, nQmD, nDmQ, nD, nU, nQ) {
    numerator <- as.numeric(nDQ) * as.numeric(nUmDQ) - 
                 as.numeric(nQmD) * as.numeric(nDmQ)
    denominator <- sqrt(
        as.numeric(nD) * (nU - nD) * nQ * (nU - nQ)
    )
    numerator / denominator
}

## Calculate Normalized Pointwise Mutual Information
calculate_npmi <- function(nD, nQ, nDQ, nU) {
    (log2(nD) + log2(nQ) - 2 * log2(nU)) / 
    (log2(nDQ) - log2(nU)) - 1
}
