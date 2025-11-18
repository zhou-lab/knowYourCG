
#' Test enrichment from YAME-compressed CG sets
#' 
#' @param query_fn File path to query
#' @param knowledge_fn File path to knowledgebase
#' @param universe_fn optional file path to universe
#' @param alternative greater, less
#' @return A single concatenated string.
#' @useDynLib knowYourCG, .registration = TRUE
#' @export
#' @examples
#' if (.Platform$OS.type!="windows") {
#' kfn = system.file("extdata", "chromhmm.cm", package = "knowYourCG")
#' qfn = system.file("extdata", "onecell.cg", package = "knowYourCG")
#' testEnrichment2(qfn, kfn)
#' }
testEnrichment2 <- function(
    query_fn, knowledge_fn, universe_fn=NULL, alternative="greater") {

    stopifnot(is.character(query_fn), is.character(knowledge_fn))
    if (.Platform$OS.type == "windows") {
        stop("Testing sequencing data does not support Windows.")
    }
    yame_result <- .Call("yame_summary_cfunc", query_fn, knowledge_fn)
    
    ## args <- c("summary", "-m", knowledgebase_fn, query_fn)
    ## yame_result <- system2("yame", args, stdout = TRUE, stderr = TRUE)
    
    ## command <- paste("yame summary -m", knowledgebase_fn, query_fn)
    ## yame_result <- system(command, intern = TRUE)
    df <- tibble::as_tibble(read.table(
        text = paste(yame_result, collapse = "\n"), header = TRUE))
    res <- cbind(df, testEnrichmentFisherN(
        nD = df$N_mask, nQ = df$N_query, nDQ = df$N_overlap, nU = df$N_univ))
    res <- res[!is.na(res$Mask),]
    res
}
