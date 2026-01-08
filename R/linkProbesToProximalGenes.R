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
#' probes = rowData(
#'     sesameData::sesameDataGet('MM285.tissueSignature'))$Probe_ID[1:10]
#' linkProbesToProximalGenes(probes, platform = "MM285")
#' @export
linkProbesToProximalGenes <- function(
    probeIDs, platform = NULL, genome = NULL) {
    platform <- queryCheckPlatform(platform, probeIDs, silent = FALSE)
    if (is.null(genome)) {
        genome <- ifelse(platform == "MM285","mm10","hg38")
    }
    genes <- sesameData_getTxnGRanges(genome, merge2gene = TRUE)
    sesameData_annoProbes(
        probeIDs, genes, return_ov_features=TRUE)
}
