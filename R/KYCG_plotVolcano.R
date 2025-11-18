#' creates a volcano plot of -log2(p.value) and log(estimate)
#' given data with fields estimate and p.value.
#'
#' @param df DataFrame where each field is a database name with two fields
#' for the estimate and p.value.
#' @param label_by column in df to be used as the label (default: dbname)
#' @param alpha Float representing the cut-off alpha value for the plot. 
#' Optional. (Default: 0.05)
#' @return ggplot volcano plot
#' @import ggplot2
#' @import ggrepel
#' @examples
#' 
#' KYCG_plotVolcano(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=runif(10,10,20),
#'   overlap=as.integer(runif(10,0,30)), group="g", dbname=seq_len(10)))
#'
#' @export
KYCG_plotVolcano <- function(df, label_by="dbname", alpha=0.05) {
    ## suppress R CMD CHECK no visible binding warning
    estimate <- FDR <- label <- NULL

    ## volcano plot cannot plot extreme effect size
    df <- df[abs(df$estimate) < 1000,]
    df[["-log10(FDR)"]] <- -log10(df$FDR)
    df$Significance <- ifelse(
        df$FDR < alpha, "Significant", "Not significant")
    ## TODO: replace with column specifying sig vs non sig
    g <- ggplot(data = df,
        aes(x = .data[["estimate"]], y = .data[["-log10(FDR)"]],
            color = .data[["Significance"]]))
    g <- g + geom_point() + xlab("log2(OR)")
    g <- g + ylab("-log10 FDR") +
        scale_colour_manual(
            name = sprintf("Significance (q < %s)", alpha),
            values = c("Significant" = "red", "Not significant" = "black"))
    requireNamespace("ggrepel")
    g <- g + ggrepel::geom_text_repel(
        data = df[df$FDR < alpha & df$estimate > 0,],
        aes(label = .data[[label_by]]), size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines"),
        show.legend = FALSE)
    g
}
