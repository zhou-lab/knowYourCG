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
    options(ggrepel.max.overlaps = 10)

    if ("Target" %in% colnames(data))
        data["label"] = unlist(data[["Target"]])
    else
        data["label"] = rownames(data)

    if (is.na(title)) {
        title = "Volcano plot"
    }
    title = gsub('(.{1,80})(\\s|$)', '\\1\n', title)

    if (is.na(subtitle)) {
        subtitle = ''
    }

    subtitle = gsub('(.{1,80})(\\s|$)', '\\1\n', subtitle)

    if (any(data$p.value <= 0.05)) {
        g = ggplot(data=data, aes(x=log2(estimate), y=-log10(p.value),
                                  color = cut(p.value, c(-Inf, 0.05))))
    } else {
        g = ggplot(data=data, aes(x=estimate, y=p.value))
    }
    g = g + geom_point() +
    xlab("log2 Fold Change") +
    ylab("-log10 p-value") +
    labs(
        title = title,
         subtitle = subtitle,
         legend= "HELLO",
         fill = "pvalue"
        ) +
    theme(
        plot.title = element_text(size=16, face = "bold"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)
        ) +
    scale_colour_discrete(
        name = "Significance (alpha = 0.05)",
        labels=c("Significant", "Not Significant")
        ) +
    geom_text_repel(
        data = subset(data, p.value < 0.05),
        aes(label = label),
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
    if ("Target" %in% colnames(data))
        data["label"] = unlist(data[["Target"]])
    else
        data["label"] = rownames(data)

    data = head(data[order(data$estimate, decreasing=TRUE), ], n=n)

    if (is.na(title)) {
        title = 'Lollipop Plot'
    }

    if (is.na(subtitle)) {
        subtitle = ''
    }

    ggplot(data, aes(x=label, y=log2(estimate), label=sprintf('%.2f',log2(estimate)))) +
        geom_hline(yintercept=0) +
        geom_segment(aes(y=0, x=reorder(label, -estimate), yend=log2(estimate), xend=label), color='black') +
        geom_point(aes(fill=pmax(-1.5,log2(estimate))), stat='identity', size=10, alpha=0.95, shape=21) +
        scale_fill_gradientn(name='Fold Change',
                             colours=c('#2166ac','#333333','#b2182b'),
                             limits=c(0,
                                      max(log2(data$estimate + 1)))) +
        geom_text(color='white', size=3) +
        labs(title=title, subtitle=subtitle) +
        geom_label(aes(x=label,
                       y=ifelse(estimate>1,
                                log2(estimate) + 0.8,
                                log2(estimate) - 0.5),
                       label=label),
                   alpha=0.8) +
        # new_scale("fill") +
        # scale_fill_manual(values=hmm_colors) +
        ylab("Log2 Enrichment") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
}

#' createGeneNetwork creates databaseSet network using the given similarity
#' metric.
#'
#' @param databaseSet Vector of probes corresponding to a single database set
#' of interest.
#' @param metric String representing the similarity score to use. Optional.
#' (Default: "Jaccard").
#'
#' @return ggplot lollipop plot
#'
#' @import RCy3
#' @import reshape2
#'
#' @examples
#' createDatabaseSetNetwork(list(a=c("a", "b"), b=c("a", "e", "f"), c=c("q", "a")))
#'
#' @export
createDatabaseSetNetwork = function(databaseSets, title="Database Interaction Network", collection="DatabaseSets") {
    m = getDatabaseSetPairwiseDistance(databaseSets, metric=metric)
    saveRDS(m, "/Users/ethanmoyer/Dropbox/Ongoing_knowYourCpG/data/databaseSetNetwork.rds")

    m_ = m
    m = m_[1:50, 1:50]

    m_melted = melt(m); colnames(m_melted) = c("gene1", "gene2", "metric")
    m_melted = m_melted[m_melted$metric != 0, ]

    # Used for additional attributes like color, size, name. This is for GSM
    nodes <- data.frame(id=colnames(m),
                        # group=c("A","A","B","B"), # categorical strings
                        # score=as.integer(c(20,10,15,5)), # integers
                        stringsAsFactors=FALSE)
    # This is for Target
    edges <- data.frame(source=m_melted$gene1,
                        target=m_melted$gene2,
                        # interaction=NULL, Maybe for positive/negative assocation
                        weight=m_melted$metric, # numeric
                        stringsAsFactors=FALSE)

    # Return nodes and edges
    return(list(nodes=nodes, edges=edges))

    # createNetworkFromDataFrames(nodes, edges, title=title, collection=collection)
    #
    # filename = file.path("images","text.png")
    # if(file.exists(filename)){
    #     file.remove(filename)
    # }
    #
    # layoutNetwork('force-directed defaultSpringLength=70 defaultSpringCoefficient=0.000003')
    #
    # #export the network
    # exportImage(filename, type = "png")

}

