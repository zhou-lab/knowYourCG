#' plotVolcano creates a volcano plot of -log2(p.value) and log(estimate)
#' given data with fields estimate and p.value.
#'
#' @param data DataFrame where each field is a database name with two fields
#' for the estimate and p.value.
#' @param title String representing the title label. Optional. (Default: NA)
#' @param subtitle String representing the subtitle label. Optional. (Default:
#' NA)
#'
#' @return ggplot volcano plot
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @examples
#' data=data.frame(estimate=c(runif(10)), p.value=c(runif(10)))
#' plotVolcano(data)
#'
#' @export
plotVolcano = function(data, title=NA, subtitle=NA) {
    data["gsm"] = rownames(data)
    # TODO: create mapping between GSM and gene name
    #
    if (is.na(title)) {
        title = "Volcano plot"
    }
    title = gsub('(.{1,80})(\\s|$)', '\\1\n', title)

    if (is.na(subtitle)) {
        subtitle = ''
    }

    subtitle = gsub('(.{1,80})(\\s|$)', '\\1\n', subtitle)

    if (any(data$p.value <= 0.05)) {
        g = ggplot(data=data,
                   aes(
                       x=log2(estimate),
                       y=-log10(p.value),
                       color = cut(p.value, c(-Inf, 0.05))))
    } else {
        g = ggplot(data=data, aes(x=estimate, y=p.value))
    }
    g = g + geom_point() +
    xlab("log2 Fold Change") +
    ylab("-log10 p-value") +
    labs(title = title,
         subtitle = subtitle,
         legend= "HELLO",
         fill = "pvalue") +
    theme(plot.title = element_text(size=16, face = "bold"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)) +
    scale_colour_discrete(name = "Significance (alpha = 0.05)",
        labels=c("Significant", "Not Significant")) +
    geom_text_repel(
        data = subset(data, p.value < 0.05),
        aes(label = gsm),
        size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
        )
    g
}

#' plotLollipop creates a lollipop plot of log(estimate) given data with fields
#' estimate.
#'
#' @param data DataFrame where each field is a database name with its estimate.
#' @param n Integer representing the number of top enrichments to report.
#' Optional. (Default: 10)
#' @param title String representing the title label. Optional. (Default: NA)
#' @param subtitle String representing the subtitle label. Optional. (Default:
#' NA)
#'
#' @return ggplot lollipop plot
#'
#' @import ggplot2
#'
#' @examples
#' data=data.frame(estimate=c(runif(10, 0, 10)))
#' plotLollipop(data)
#'
#' @export
plotLollipop = function(data, n=10, title=NA, subtitle=NA) {
    data$x = rownames(data)
    data = head(data[order(data$estimate, decreasing=TRUE), ], n=n)

    if (is.na(title)) {
        title = 'Lollipop Plot'
    }

    if (is.na(subtitle)) {
        subtitle = ''
    }

    ggplot(data, aes(x=x, y=log2(estimate), label=sprintf('%.2f',log2(estimate)))) +
        geom_hline(yintercept=0) +
        geom_segment(aes(y=0, x=reorder(x,- estimate), yend=log2(estimate), xend=x), color='black') +
        geom_point(aes(fill=pmax(-1.5,log2(estimate))), stat='identity', size=10, alpha=0.95, shape=21) +
        scale_fill_gradientn(name='Fold Change',
                             colours=c('#2166ac','#333333','#b2182b'),
                             limits=c(0,
                                      max(log2(data$estimate + 1)))) +
        geom_text(color='white', size=3) +
        labs(title=title, subtitle=subtitle) +
        geom_label(aes(x=x,
                       y=ifelse(estimate>1,
                                log2(estimate) + 0.8,
                                log2(estimate) - 0.5),
                       label=x),
                   alpha=0.8) +
        # new_scale("fill") +
        # scale_fill_manual(values=hmm_colors) +
        ylab("Log2 Enrichment") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
}

