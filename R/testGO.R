convertGeneName <- function(gene) {
    first <- substr(gene, 1, 1)
    rest <- substr(gene, 2, nchar(gene))
    paste(toupper(first), tolower(rest), sep = "")
}


#' tests Gene Ontology of genes overlapping CpG query
#'
#' estimate represent enrichment score and negative estimate indicate a
#' test for depletion
#'
#' @param probeIDs Vector of CpG probes IDs or a data frame with 
#' gene_name column, usually the output of testEnrichment() function
#' @param platform EPIC, MM285, ..., infer if not given
#' @param organism The organism corresponding to the CpG platform
#' or genes in gene_name column
#' @param gene_name If query is data frame from testEnrichment output,
#' whether to use the gene_name column. If set to FALSE,
#' TFBS will be used (default: FALSE)
#' @param ... Additional arguments to sesameData_getGenesByProbes and gost()
#' @return A list of enriched terms and meta data from gprofiler2 output
#' @examples
#' library(SummarizedExperiment)
#' sesameData::sesameDataCache(data_titles=
#' c("KYCG.MM285.tissueSignature.20211211","probeIDSignature",
#' "MM285.address","genomeInfo.mm10"))
#' df <- rowData(sesameData::sesameDataGet('MM285.tissueSignature'))
#' query <- df$Probe_ID[df$branch == "fetal_liver" & df$type == "Hypo"]
#' res <- testGO(query,platform="MM285")
#' @export
testGO <- function(
    probeIDs, platform=NULL, organism="hsapiens",gene_name=TRUE,...) {
  
    if (!requireNamespace("gprofiler2",quietly = TRUE)) {
      stop("Install 'gprofiler2' to use this function")
    }
  
    platform <- queryCheckPlatform(platform, probeIDs, silent = FALSE)

    if (is.character(probeIDs)) {
        genome <- ifelse(platform == "MM285","mm10","hg38")
        genes <- sesameData_getTxnGRanges(genome,merge2gene = TRUE)
        query <- sesameData_annoProbes(
          probeIDs, genes, return_ov_features = TRUE
        )
        query <- query$gene_name
    }

    if (is.data.frame(query)) {
        if (gene_name) {
            query <- query[["gene_name"]]
        } else {
            ind <- grepl("TFBS",query[["group"]])
            query <- query[["dbname"]][ind]
            if (organism == "mmusculus") {
                    query <- convertGeneName(query)
            }
        }
    }

    gostres <- gprofiler2::gost(query, organism = organism,...)
    gostres$result <- gostres$result[order(gostres$result$p_value),]
    gostres
}

