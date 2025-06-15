convertGeneName <- function(gene) {
    first <- substr(gene, 1, 1)
    rest <- substr(gene, 2, nchar(gene))
    paste(toupper(first), tolower(rest), sep = "")
}

#' find genes in genomic proximity to given Infinium probes
#'
#' This is a convenient function that uses sesameData_getGenomeInfo()
#' to retrieve stored gene models.
#' 
#' For finer control, such as taking only genes by their promoters,
#' please use sesameData_getTxnGRanges followed by
#' sesameData_annoProbes(). See code of this convenient function for details.
#' 
#' @param probeIDs character vector of probe IDs
#' @param platform HM450, EPIC, EPICv2, MM285, MSA, ..., will infer
#' from probe ID if not given
#' @param genome hg38, hg19, mm10, this is usually inferred from platform. 
#' @return a data frame annotate gene list linked to each given probes
#' @examples
#' library(SummarizedExperiment)
#' probes = rowData(sesameData::sesameDataGet('MM285.tissueSignature'))$Probe_ID[1:10]
#' linkProbesToProximalGenes(probes, platform = "MM285")
#' @export
linkProbesToProximalGenes <- function(probeIDs, platform = NULL, genome = NULL) {
    platform <- queryCheckPlatform(platform, probeIDs, silent = FALSE)
    if (is.null(genome)) {
        genome <- ifelse(platform == "MM285","mm10","hg38")
    }
    genes <- sesameData_getTxnGRanges(genome, merge2gene = TRUE)
    sesameData_annoProbes(
        probeIDs, genes, return_ov_features=TRUE)
}

## #' tests Gene Ontology of genes overlapping CpG query
## #'
## #' estimate represent enrichment score and negative estimate indicate a
## #' test for depletion
## #'
## #' @param probeIDs Vector of CpG probes IDs or a data frame with 
## #' gene_name column, usually the output of testEnrichment() function
## #' @param platform EPIC, MM285, ..., infer if not given
## #' @param organism The organism corresponding to the CpG platform
## #' or genes in gene_name column
## #' @param gene_name If query is data frame from testEnrichment output,
## #' whether to use the gene_name column. If set to FALSE,
## #' TFBS will be used (default: FALSE)
## #' @param ... Additional arguments to sesameData_getGenesByProbes and gost()
## #' @return A list of enriched terms and meta data from gprofiler2 output
## #' @examples
## #' library(SummarizedExperiment)
## #' sesameData::sesameDataCache(data_titles=
## #' c("MM285.tissueSignature","probeIDSignature",
## #' "MM285.address","genomeInfo.mm10"))
## #'
## #' df <- rowData(sesameData::sesameDataGet('MM285.tissueSignature'))
## #' query_probes <- df$Probe_ID[df$branch == "fetal_liver" & df$type == "Hypo"]
## #' genes <- linkProbesToGenes(query_probes, platform = "MM285")
## #' res <- testGO(query,platform="MM285")
## #' @export
## testGO <- function(
##     probeIDs, platform=NULL, organism="hsapiens",gene_name=TRUE,...) {
  
##     if (!requireNamespace("gprofiler2",quietly = TRUE)) {
##       stop("Install 'gprofiler2' to use this function")
##     }
  
##     platform <- queryCheckPlatform(platform, probeIDs, silent = FALSE)

##     if (is.character(probeIDs)) {
##         query <- linkProbesToGenes(probeIDs, platform)$gene_name
##     }

##     if (is.data.frame(query)) { # output of testEnrichment
##         if (gene_name) {
##             query <- query[["gene_name"]]
##         } else {
##             ind <- grepl("TFBS",query[["group"]])
##             query <- query[["dbname"]][ind]
##             if (organism == "mmusculus") {
##                 query <- convertGeneName(query)
##             }
##         }
##     }

##     gostres <- gprofiler2::gost(query, organism = organism,...)
##     gostres$result <- gostres$result[order(gostres$result$p_value),]
##     gostres
## }

