#' Test enrichment from YAME-compressed CG sets
#' 
#' Tests for enrichment of genomic regions in YAME-compressed CG sets using
#' Fisher's exact test. This function is optimized for sequencing data and
#' uses compiled C code for efficient processing.
#' 
#' @param query_fn Character string specifying the file path to the query
#' CG set file (YAME-compressed format).
#' @param knowledge_fn Character string or vector specifying the file path(s)
#' to the knowledgebase file(s) (YAME-compressed format). Can be a single
#' file or multiple files.
#' @param universe_fn Optional character string specifying the file path to
#' the universe CG set file. If NULL, universe will be inferred from the
#' knowledgebase. (Default: NULL)
#' @param alternative Character string specifying the alternative hypothesis:
#' "greater" (enrichment), "less" (depletion), or "two.sided". 
#' (Default: "greater")
#' @param min_overlap Minimum number of overlapping CGs required for a test
#' to be included in results. (Default: 1)
#' @param verbose Logical indicating whether to print progress messages.
#' (Default: FALSE)
#' 
#' @return A tibble containing enrichment test results with the following columns:
#' \describe{
#'   \item{Mask}{Name/identifier of the knowledgebase mask}
#'   \item{N_mask}{Number of CGs in the mask}
#'   \item{N_query}{Number of CGs in the query}
#'   \item{N_overlap}{Number of overlapping CGs}
#'   \item{N_univ}{Total number of CGs in universe}
#'   \item{estimate}{Log2 odds ratio}
#'   \item{p.value}{P-value from Fisher's exact test}
#'   \item{log10.p.value}{Log10-transformed p-value}
#'   \item{test}{Type of test performed}
#'   \item{Additional effect size metrics}{Jaccard, MCC, etc.}
#' }
#' 
#' @useDynLib knowYourCG, .registration = TRUE
#' @importFrom tibble as_tibble
#' @importFrom utils read.table
#' 
#' @export
#' 
#' @examples
#' if (.Platform$OS.type != "windows") {
#'   kfn <- system.file("extdata", "chromhmm.cm", package = "knowYourCG")
#'   qfn <- system.file("extdata", "onecell.cg", package = "knowYourCG")
#'   res <- testEnrichment2(qfn, kfn)
#'   head(res)
#' }
testEnrichment2 <- function(
    query_fn, 
    knowledge_fn, 
    universe_fn = NULL, 
    alternative = "greater",
    min_overlap = 1,
    verbose = FALSE) {
    
    ## Validate inputs
    validate_inputs(query_fn, knowledge_fn, universe_fn)
    
    ## Check platform compatibility
    if (.Platform$OS.type == "windows") {
        stop(
            "Testing sequencing data is not currently supported on Windows. ",
            "This feature requires compiled C code that is not available on Windows.",
            call. = FALSE
        )
    }
    
    ## Validate alternative parameter
    alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
    
    if (verbose) {
        message("Reading query file: ", query_fn)
        message("Reading knowledgebase file(s): ", 
                paste(knowledge_fn, collapse = ", "))
    }
    
    ## Call compiled C function to compute overlaps
    yame_result <- tryCatch(
        .Call("yame_summary_cfunc", query_fn, knowledge_fn),
        error = function(e) {
            stop(
                "Failed to process YAME files: ", e$message,
                "\nPlease ensure files are in valid YAME-compressed format.",
                call. = FALSE
            )
        }
    )
    
    ## Parse results into data frame
    df <- parse_yame_results(yame_result)
    
    if (nrow(df) == 0) {
        warning("No results returned from YAME processing.", call. = FALSE)
        return(tibble::tibble())
    }
    
    if (verbose) {
        message(sprintf("Processing %d mask(s)...", nrow(df)))
    }
    
    ## Perform Fisher's exact test for each mask
    res <- compute_enrichment_stats(df, alternative, min_overlap)
    
    ## Remove rows with missing mask names
    res <- res[!is.na(res$Mask) & res$Mask != "", ]
    
    if (nrow(res) == 0) {
        warning("No valid results after filtering.", call. = FALSE)
        return(tibble::tibble())
    }
    
    ## Order by significance
    res <- res[order(res$log10.p.value, -abs(res$estimate)), ]
    
    if (verbose) {
        message(sprintf(
            "Completed enrichment testing for %d mask(s).", 
            nrow(res)
        ))
    }
    
    res
}

## Validate input file paths
validate_inputs <- function(query_fn, knowledge_fn, universe_fn) {
    ## Check input types
    if (!is.character(query_fn) || length(query_fn) != 1) {
        stop("'query_fn' must be a single character string.", call. = FALSE)
    }
    
    if (!is.character(knowledge_fn)) {
        stop("'knowledge_fn' must be a character string or vector.", call. = FALSE)
    }
    
    if (!is.null(universe_fn) && !is.character(universe_fn)) {
        stop("'universe_fn' must be NULL or a character string.", call. = FALSE)
    }
    
    ## Check file existence
    if (!file.exists(query_fn)) {
        stop("Query file not found: ", query_fn, call. = FALSE)
    }
    
    missing_kb <- knowledge_fn[!file.exists(knowledge_fn)]
    if (length(missing_kb) > 0) {
        stop(
            "Knowledgebase file(s) not found: ",
            paste(missing_kb, collapse = ", "),
            call. = FALSE
        )
    }
    
    if (!is.null(universe_fn) && !file.exists(universe_fn)) {
        stop("Universe file not found: ", universe_fn, call. = FALSE)
    }
}

## Parse YAME results into a tibble
parse_yame_results <- function(yame_result) {
    if (is.null(yame_result) || length(yame_result) == 0) {
        stop("Empty result from YAME processing.", call. = FALSE)
    }
    
    ## Convert to tibble
    df <- tryCatch(
    {
        text_data <- paste(yame_result, collapse = "\n")
        temp_df <- read.table(
            text = text_data, 
            header = TRUE,
            stringsAsFactors = FALSE,
            comment.char = ""
        )
        tibble::as_tibble(temp_df)
    },
    error = function(e) {
        stop(
            "Failed to parse YAME output: ", e$message,
            "\nOutput may be malformed.",
            call. = FALSE
        )
    }
    )
    
    ## Validate expected columns
    required_cols <- c("Mask", "N_mask", "N_query", "N_overlap", "N_univ")
    missing_cols <- setdiff(required_cols, names(df))
    
    if (length(missing_cols) > 0) {
        stop(
            "YAME output missing required columns: ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
        )
    }
    
    df
}

## Compute enrichment statistics for each mask
compute_enrichment_stats <- function(df, alternative, min_overlap) {
    ## Filter by minimum overlap if specified
    if (min_overlap > 1) {
        df <- df[df$N_overlap >= min_overlap, ]
    }
    
    if (nrow(df) == 0) {
        warning(
            "No masks meet minimum overlap threshold (", 
            min_overlap, ").",
            call. = FALSE
        )
        return(tibble::tibble())
    }
    
    ## Available from the main enrichment file
    ## It calculates Fisher's exact test statistics from count data
    enrichment_stats <- testEnrichmentFisherN(
        nD = df$N_mask,
        nQ = df$N_query,
        nDQ = df$N_overlap,
        nU = df$N_univ,
        alternative = alternative
    )
    
    ## Combine with original data
    res <- cbind(df, enrichment_stats)
    
    ## Add MFile column to track source if needed
    if (!"MFile" %in% names(res)) {
        res$MFile <- "YAME"
    }
    
    tibble::as_tibble(res)
}
