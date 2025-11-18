#' Plot meta gene or other meta genomic features
#'
#' @param betas a named numeric vector or a matrix
#' (row: probes; column: samples)
#' @param platform if not given and x is a SigDF, will be inferred
#' the meta features
#' @importFrom reshape2 melt
#' @importFrom sesameData sesameDataGet
#' @return a grid plot object
#' @examples
#' library(sesameData)
#' library(sesame)
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' KYCG_plotMeta(getBetas(sdf))
#' @export
KYCG_plotMeta <- function(betas, platform = NULL) {

    if (!is.matrix(betas)) {
        betas <- cbind(sample=betas)
    }
    if (is.null(platform)) {
        platform <- inferPlatformFromProbeIDs(rownames(betas))
    }
    stopifnot(!is.null(platform))

    dbs <- getDBs(sprintf("%s.metagene", platform))
    df <- dbStats(betas, dbs, long=TRUE)
    dflabel <- data.frame(
        ord = as.integer(names(dbs)),
        reg = vapply(dbs, function(x) attr(x, "label"), character(1)))
    
    ggplot(df) +
        annotate("rect", xmin = -1, xmax = 10, ymin = -Inf,
            ymax = Inf, fill = "grey80", alpha = .5, color = NA) +
        geom_line(aes(.data$db, .data$value, group=.data$query)) +
        scale_x_continuous(breaks=dflabel$ord, labels=dflabel$reg) +
        ylab("Mean DNA Methylation Level") + xlab("") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
