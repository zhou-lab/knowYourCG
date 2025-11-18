#' Bar plot to show most enriched CG groups from testEnrichment
#'
#' The input data frame should have an "estimate" and
#' a "FDR" columns.
#' 
#' Top CG groups are determined by estimate (descending order).
#'
#' @param df KYCG result data frame
#' @param y the column to be plotted on y-axis
#' @param n number of CG groups to plot
#' @param order_by the column by which CG groups are ordered
#' @param label whether to label significant bars
#' @return grid plot object
#'
#' @import ggplot2
#' @examples
#' KYCG_plotBar(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=10,
#'   overlap=as.integer(runif(10,0,30)), group="g", dbname=seq_len(10)))
#' @export
KYCG_plotBar <- function(df, y = "-log10(FDR)",
    n = 20, order_by = "FDR", label = FALSE) {

    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))
    df1 <- preparePlotDF(df, n, order_by)
    if (y == "-log10(FDR)") {
        df1[["-log10(FDR)"]] <- -log10(df1$FDR)
    }

    p <- ggplot(df1, aes(.data[["db1"]], .data[[y]])) +
        geom_bar(stat="identity") +
        coord_flip() + ylab(y) + xlab("CpG Group")
    
    if (label) {
        ## only significant ones are labeled
        df1_label <- df1[df1$FDR < 0.05,]
        df1_label$pos_label <- df1_label[[y]]/2
        df1_label$label <- sprintf("N=%d", df1_label$overlap)
        p <- p + geom_label(aes(
            x=.data[["db1"]], y=.data[["pos_label"]],
            label = .data[["label"]]),
            data = df1_label, alpha=0.6, hjust=0.5)
    }
    p
}
