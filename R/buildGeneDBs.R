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
