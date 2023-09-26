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
#' @param query Vector of CpG probes IDs or a data frame with gene_name column,
#' usually the output of testEnrichment() function
#' @param platform EPIC, MM285, ..., infer if not given
#' @param organism The organism corresponding to the CpG platform
#' or genes in gene_name column
#' @param gene_name If query is data frame from testEnrichment output,
#' whether to use the gene_name column. If set to FALSE,
#' TFBS will be used (default: FALSE)
#' @param ... Additional arguments to sesameData_getGenesByProbes and gost()
#' @importFrom gprofiler2 gost
#' @return A list of enriched terms and meta data from gprofiler2 output
#' @examples
#' library(SummarizedExperiment)
#' df <- rowData(sesameDataGet('MM285.tissueSignature'))
#' query <- df$Probe_ID[df$branch == "fetal_liver" & df$type == "Hypo"]
#' res <- testGO(query)
#' @export
testGO <- function(
        query, platform=NULL, organism="hsapiens",gene_name=TRUE,...) {

    if (is.character(query)) {
        query <- sesameData_getGenesByProbes(query, platform=platform,...)
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

    gostres <- gost(query, organism = organism,...)
    gostres$result <- gostres$result[order(gostres$result$p_value),]
    gostres
}

