#' Plot point range for a list of enrichment testing results
#' against the same set of databases
#'
#' @param result_list a list of testEnrichment resultsx
#' @return grid plot object
#' @importFrom reshape2 melt
#' @importFrom dplyr summarize
#' @examples
#'
#' ## pick some big TFBS-overlapping CpG groups
#' cg_lists <- getDBs("MM285.TFBS")
#' queries <- cg_lists[(sapply(cg_lists, length) > 40000)]

#' result_list <- lapply(queries, testEnrichment,
#'     "MM285.chromHMM", platform="MM285")
#' KYCG_plotPointRange(result_list)
#' 
#' @export
KYCG_plotPointRange <- function(result_list) {

    ord <- mean_betas <- state <- est <- NULL
    
    mtx <- aggregateTestEnrichments(result_list)
    df <- melt(mtx, varnames = c("sample","state"), value.name = "est")
    df <- summarize(group_by(df, state), 
        ave = mean(pmax(-4, est),na.rm=TRUE),
        sd = sd(pmax(-10,est),na.rm=TRUE))
    df$ymin <- df$ave - df$sd
    df$ymax <- df$ave + df$sd
    
    df$state <- factor(df$state, levels = df$state[order(df$ave)])
    
    ggplot(df) +
        geom_pointrange(aes_string("state", "ave", ymin="ymin", ymax="ymax")) +
        geom_hline(yintercept=0, linetype='dashed') +
        ylab("Log2 Fold Enrichment") + xlab("") +
        scale_y_continuous(position="right") +
        annotate("text", x=0.5, y=0.5,
            label="Enrichment", angle=-90, hjust=1) +
        annotate("text", x=0.5, y=-0.5,
            label="Depletion", angle=90, hjust=0) +
        coord_flip()
}

