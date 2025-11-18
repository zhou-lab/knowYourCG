#' create a waterfall plot of log(estimate) given test enrichment
#'
#' @param df data frame where each row is a database with test
#' enrichment result
#' @param order_by the column by which CG groups are ordered
#' @param size_by the column by which CG group size plot
#' @param n_label number of datapoints to label
#' @param label_by column in df to be used as the label (default: dbname)
#' @return grid
#' @importFrom utils head
#' @import ggplot2
#' @import ggrepel
#' @importFrom sesameData sesameDataGet
#' @examples
#'
#' library(SummarizedExperiment)
#' library(sesameData)
#' df <- rowData(sesameDataGet('MM285.tissueSignature'))
#' query <- df$Probe_ID[df$branch == "fetal_brain" & df$type == "Hypo"]
#' results <- testEnrichment(query, "TFBS", platform="MM285")
#' KYCG_plotWaterfall(results)
#' 
#' @export
KYCG_plotWaterfall <- function(df,
    order_by="Log2(OR)", size_by="-log10(FDR)",
    label_by="dbname", n_label=10) {

    df$label <- df[[label_by]]
    if (size_by == "-log10(FDR)" ||
        order_by == "-log10(FDR)" ||
        label_by == "-log10(FDR)") {
        df[["-log10(FDR)"]] <- -log10(df$FDR)
    }
    if (df$test[[1]] == "Log2(OR)" && (
        size_by == "Log2(OR)" || order_by == "Log2(OR)" ||
        label_by == "Log2(OR)")) {
        df[["Log2(OR)"]] <- df$estimate
        message(sprintf("%d extremes are capped.",
            sum(abs(df[["Log2(OR)"]]) > 1000)))
        ## cap extremes
        df[["Log2(OR)"]][df[["Log2(OR)"]] > 1000] <- 1000
        df[["Log2(OR)"]][df[["Log2(OR)"]] < -1000] <- -1000
        ## df <- df[abs(df$estimate) < 1000,] # skip extremes
    }

    df <- df[order(df[[order_by]]),]
    df$index <- seq_len(nrow(df))
    
    requireNamespace("ggrepel")
    ggplot(df, aes(.data[["index"]], .data[[order_by]])) +
        geom_point(aes(size=.data[[size_by]]), alpha=0.6) +
        geom_hline(yintercept=0, linetype="dashed", color="grey60") +
        theme_minimal() + ylab(order_by) + xlab("Databases") +
        ggrepel::geom_text_repel(
            data = df[head(order(df$log10.p.value),
                n = min(n_label, nrow(df)*0.5)),],
            aes_string(label="label"), nudge_x=-nrow(df)/10,
            max.overlaps=999)
}

