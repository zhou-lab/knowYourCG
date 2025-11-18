#' Plot meta gene or other meta genomic features
#'
#' @param result_list one or a list of testEnrichment
#' @return a grid plot object
#' @examples
#' cg_lists <- getDBs("MM285.TFBS")
#' queries <- cg_lists[(sapply(cg_lists, length) > 40000)]
#' result_list <- lapply(queries, testEnrichment,
#'     "MM285.metagene", silent=TRUE, platform="MM285")
#' 
#' KYCG_plotMetaEnrichment(result_list)
#' @export
KYCG_plotMetaEnrichment <- function(result_list) {

    if (is.data.frame(result_list)) { # 1 testEnrichment result
        result_list <- list(result_list)
    }

    stopifnot(all(c("dbname", "label") %in% colnames(result_list[[1]])))
    df <- aggregateTestEnrichments(result_list, return_df = TRUE)
    
    ggplot(df) +
        annotate("rect", xmin = -1, xmax = 10, ymin = -Inf,
            ymax = Inf, fill = "grey80", alpha = .5, color = NA) +
        geom_line(aes_string("db", "estimate", color="query")) +
        scale_x_continuous(breaks=as.integer(result_list[[1]]$db),
            labels=result_list[[1]]$label) +
        annotate("text", x=min(as.integer(result_list[[1]]$db)),
            y=0.05, label="Enrichment", hjust = 0, vjust = 0) +
        annotate("text", x=min(as.integer(result_list[[1]]$db)),
            y=-0.05, label="Depletion", hjust = 0, vjust = 1) +
        geom_hline(yintercept = 0.0, linetype="dashed") +
        ylab("Log2 Fold Enrichment") + xlab("") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
