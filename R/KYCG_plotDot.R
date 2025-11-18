#' Dot plot to show most enriched CG groups from testEnrichment
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
#' @param size_by the column by which CG group size plot
#' @param color_by the column by which CG groups are colored
#' @param label_by the column for label
#' @param short_label omit group in label
#' @param title plot title
#' @return grid plot object (by ggplot)
#'
#' @import ggplot2
#' @examples
#' KYCG_plotDot(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=runif(10,10,20),
#'   overlap=as.integer(runif(10,0,30)), group="g", dbname=seq_len(10)))
#' @export
KYCG_plotDot <- function(df, y = "-log10(FDR)",
    n = 20, order_by = "FDR", title = "Enriched Knowledgebases",
    label_by = "dbname", size_by = "overlap", color_by = "estimate",
    short_label = FALSE) {

    df1 <- preparePlotDF(
        df, n, order_by, short_label = short_label, label_by = label_by)

    if (y == "-log10(FDR)") {
        df1[["-log10(FDR)"]] <- -log10(df1$FDR)
    }
    ggplot(df1) +
        geom_point(aes(.data[["db1"]], .data[[y]],
            size = .data[[size_by]], color = .data[[color_by]])) +
        coord_flip() + ggtitle(title) +
        scale_color_gradient(low="blue",high="red") +
        ylab(y) + xlab("")
}
