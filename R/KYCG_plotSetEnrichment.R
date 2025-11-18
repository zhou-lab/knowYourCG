#' Plot Set Enrichment
#'
#' @param result result object as returned from an element of the list of
#' testEnrichmentSEA(..., prepPlot=TRUE)
#' @param n_sample number of CpGs to sample
#' @param n_presence number of overlap to sample for the plot
#' @return grid object for plot
#' @importFrom wheatmap WGG
#' @importFrom wheatmap Beneath
#' @examples
#' query <- getDBs("KYCG.MM285.designGroup")[["VMR"]]
#' db <- getDBs("MM285.seqContextN", "distToTSS")
#' res <- testEnrichmentSEA(query, db, prepPlot = TRUE)
#' KYCG_plotSetEnrichment(res[[1]])
#' 
#' @export
KYCG_plotSetEnrichment <- function(
    result, n_sample = 1000, n_presence = 200) {

    stopifnot("dDisc" %in% names(result))
    dCont <- sort(result$dCont)
    dDisc <- result$dDisc
    presence <- names(dCont) %in% dDisc
    cs <- cumsum(ifelse(presence, 1/sum(presence), -1/sum(!presence)))
    index <- as.integer(seq(1, length(cs), length.out=n_sample))

    pos <- which(presence)
    if(length(pos) > n_presence) { pos <- sample(pos, n_presence) }

    WGG(ggplot(data.frame(index=index, cs=cs[index])) +
        geom_segment(data=data.frame(pos=pos),
            aes(x = .data[["pos"]], xend = .data[["pos"]],
                y = -0.02, yend = 0.02),
            color="grey50") +
        geom_line(aes(x=.data[["index"]], y=.data[["cs"]]), color="darkred") +
        xlab("") + ylab("ES(S)")) +
    WGG(ggplot(data.frame(index=index, var=dCont[index]),
        aes(x=.data[["index"]], y=.data[["var"]])) +
        geom_area() +
        xlab("CpGs") + ylab("Phenotype Var"), Beneath(height=0.5))
}
