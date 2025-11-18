#' creates a lollipop plot of log(estimate) given data with
#' fields estimate.
#'
#' @param df DataFrame where each row is a database name with its
#' estimate.
#' @param label_column column in df to be used as the label (default: dbname)
#' @param n Integer representing the number of top enrichments to report.
#' Optional. (Default: 10)
#' @return ggplot lollipop plot
#' @importFrom utils head
#' @import ggplot2
#' @examples
#' 
#' KYCG_plotLollipop(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=runif(10,10,20),
#'   overlap=as.integer(runif(10,0,30)), group="g",
#'   dbname=as.character(seq_len(10))))
#' 
#' @export
KYCG_plotLollipop <- function(df, label_column="dbname", n=20) {
    ## suppress R CMD CHECK no visible binding warning
    estimate <- label <- NULL

    df$label <- df[[label_column]]
    df <- head(df[order(df$estimate, decreasing=TRUE), ], n=n)

    allest <- df$estimate[!is.infinite(df$estimate)]
    cap <- max(allest) * 1.4
    cap_line <- max(allest) * 1.2
    df$estimate[df$estimate == Inf] <- cap
    
    ggplot(df, aes(x = .data[["label"]],
        y = .data[["estimate"]], label = .data[["label"]])) +
        geom_hline(yintercept = 0) +
        geom_segment(aes(
            x=reorder(label, -estimate), y=0,
            yend=estimate, xend=label), color='black') +
        geom_point(fill="black", stat='identity', size=15,
            alpha=0.95, shape=21) +
        scale_fill_gradientn(name='Log2(OR)',
            colours=c('#2166ac','#333333','#b2182b')) +
        geom_text(color='white', size=3) +
        ylab("Log2(OR)") +
        geom_hline(yintercept = cap_line, linetype = "dashed") +
        ylim(min(min(allest)*1.3,0), max(max(allest)*1.5,0)) +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
}

