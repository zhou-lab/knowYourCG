#' plot enrichment test result
#'
#' @param df test enrichment result data frame
#' @param fdr_max maximum fdr for capping
#' @param n_label number of database to label
#' @param min_estimate minimum estimate
#' @return grid object
#' @import utils
#' @importFrom stringr str_replace
#' @importFrom tibble rownames_to_column
#' @import ggplot2
#' @examples
#' query <- KYCG_getDBs("MM285.designGroup")[["PGCMeth"]]
#' res <- testEnrichment(query, platform="MM285")
#' KYCG_plotEnrichAll(res)
#'
#' @export
KYCG_plotEnrichAll <- function(
        df, fdr_max = 25, n_label = 15, min_estimate = 0) {

    gp_size <- sort(table(df$group))
    gp_width <- log(2+gp_size)
    e1 <- df[order(factor(df$group, levels=names(gp_size)), df$dbname),]
    e1$inc <- (gp_width / gp_size)[e1$group]
    e1$inc1 <- c(0,ifelse(e1$group[-1] != e1$group[-nrow(e1)], 1, 0))
    e1$inc2 <- cumsum(e1$inc + e1$inc1)

    e1$group <- str_replace(e1$group,"KYCG.","")
    e1$group <- vapply(strsplit(e1$group, "\\."),
                       function(x) paste0(x[2:(length(x)-1)], collapse="."), character(1))
    if ("gene_name" %in% colnames(e1)) {
        e1$dbname[e1$group == "gene"] <- e1$gene_name[e1$group == "gene"] }

    e2 <- e1[e1$estimate > min_estimate & e1$FDR < 0.01 ,]
    e2$FDR[e2$FDR < 10**-fdr_max] <- 10**-(fdr_max*1.1)

    e3 <- rownames_to_column(as.data.frame(do.call(rbind, lapply(
        split(e1$inc2, e1$group), function(x)
            c(beg=min(x), middle=mean(x), end=max(x))))), "group")

    inc2 <- FDR <- estimate <- group <- dbname <- beg <- middle <- NULL
    requireNamespace("ggrepel")
    ggplot(e2, aes(inc2, -log10(FDR))) +
        geom_point(aes(size=estimate, color=group), alpha=0.5) +
        ggrepel::geom_text_repel(data = e2[head(order(e2$FDR), n = n_label),],
                                 aes(label=dbname, color=group), size = 3,
                                 ## box.padding = unit(0.35, "lines"),
                                 ## point.padding = unit(0.3, "lines"),
                                 direction="y", nudge_y=0.2, max.overlaps=100) +
        annotate("text", -1, fdr_max*0.96,
                 label="Values above this line are capped.",
                 hjust=0, vjust=1, color="grey60") +
        geom_hline(yintercept = fdr_max, linetype="dotted", color="grey60") +
        geom_segment(aes(x = beg, y = 0, xend = end, yend = 0, color=group),
                     size=3, data=e3) +
        geom_text(data=e3,aes(middle, -1, label=group, color=group),
                  vjust=1, hjust=1, angle=30) + scale_color_discrete(guide="none") +
        ylim(-6, fdr_max*1.2) + xlab("") +
        scale_size_continuous(guide=guide_legend(title="log2(OR)")) +
        coord_cartesian(clip="off") + theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.minor.x = element_blank())
}




