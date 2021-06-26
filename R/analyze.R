#' Test for enriched categories of CpGs
#' testEnrichmentAll tests a default set of categories (database sets) for enrichment
#' in the provided vector of CpGs (sigProbes).
#'
#' @param sigProbes Vector of probes, or a named vector of ranks (e.g. coefficient), with names being
#' Illumina Probe IDs
#' @param database Vector of probes, or a named vector of numeric values with names being
#' Probe IDs. sigProbes is tested against database for enrichment.
#'
#' @return A list of two tables, one each for database sets of categorical vs continuous
#' values. The tables rank database sets by odds ratio (Fisher's exact test,
#' categorical sigProbes and categorical database set), or P-value (FGSEA or Spearman, categorical vs
#'  continuous or continuous vs continuous respectively).
#'
#' @import dplyr
#' @import fgsea
#' @import sesameData
#'
#' @export
testEnrichment1 = function(sigProbes, databaseName) {

    # probe IDs: mouseMethylation285_probeIDs

    # TODO: load database sets based on mapping (excel sheet?)
    # for now, just load the default one
    database <- readRDS(url("http://zhouserver.research.chop.edu/kyCG/20210601_MM285_TFBS_ENCODE.rds"))[1:100]

    if (all(sapply(database, is.numeric))) { # numeric db
        if(is.numeric(sigProbes)) { # a named vector of continuous value
            # Spearman correlation
            spearmanResult <- bind_rows(lappy(
                names(database),
                function(db) {
                    df <- testEnrichmentSpearman(
                        sigProbes = sigProbes,
                        database = database,
                        dbName = databaseName
                    )
                    df$instanceName <- db
                    df$databaseLength <- length(db)

                    return(df)
                }
            ))

            return(spearmanResult)
        } else { # categorical query
            # FGSEA
            print('FGSEA for continuous database sets')

            ## continuous database set
            fgseaResult <- bind_rows(lapply(
                names(database),
                function(db) {
                    cbind(fgsea(
                        pathways = list(sigProbes = sigProbes),
                        stats = database[[db]]
                    ), data.frame(database = databaseName, instanceName = db))
                }
            ))

            return(fgseaResult)
        }
    } else if (all(sapply(database, is.character))) { # categorical database
        if(is.numeric(sigProbes)) { # a named vector of continuous value
            ## do fgsea(switched arguments)
            fgseaResult <- bind_rows(lapply(
                names(database),
                function(db) {
                    df <- testEnrichmentFGSEA(
                        probeIDs = probeIDs,
                        database = database,
                        sigProbes = sigProbes,
                        dbName = databaseName,
                        sigProbesRank = sigProbesRank
                    )
                    df$instanceName <- db
                    df$databaseLength <- length(database[[db]])

                    return(df)
                }
            ))

            # TODO: sort by p value, and join with database info table
            return(fgseaResult)

        } else { # categorical db
            print("Fisher's exact test for categorical database sets")
            ## if only significant probes provided, Fisher (categorical db set) and fgsea (continuous db set)
            # categorical database sets first
            fisherResult <- bind_rows(lapply(
                names(database),
                function(db) {
                    df <- testEnrichmentFisher(
                        probeIDs = mouseMethylation285_probeIDs,
                        categoryProbes = database[[db]],
                        sigProbes = sigProbes,
                        dbName = databaseName
                    )
                    df$instanceName <- db
                    df$databaseLength <- length(database[[db]])
                    df$test <- 'Fisher'

                    df
                }
            ))
            fisherResult <- fisherResult[order(fisherResult, decreasing = TRUE), ]

            # TODO: join results with database info table
            return(fisherResult)
        }
    } else {
        # This may go away depending on the final organization of database sets
        # for example, if we handle one vector at a time, there is no need for this
        stop('Database contains mixture of categorical and continuous')
    }
}



testEnrichmentAll = function(query) {
    if (!is.null(databaseSets)) {
        # TODO: support custom lists of databaseSets
    }
    # TODO: decide on what the organization of database sets will be
    # right now: vector of database set keys
    # each group of database sets is a list of vectors
    # right now, this list has only one group of database sets
    defaultDatabaseSets = c(MM285_TFBS = "http://zhouserver.research.chop.edu/kyCG/20210601_MM285_TFBS_ENCODE.rds")
    lapply(names(defaultDatabaseSets), function(db) testEnrichment1(query=query, db=db))
    ## apply multi-test correction, BH, FDR
}


# testEnrichment = function(setQ, setD, setU) {
testEnrichmentFisher = function(probeIDs, categoryProbes, sigProbes, dbName, instanceName) {
    # setD = categoryProbes
    # setQ = sigProbes
    # setU = probeIDs

    mtx = matrix(c(
        length(intersect(sigProbes, categoryProbes)),
        length(setdiff(categoryProbes, sigProbes)),
        length(setdiff(sigProbes, categoryProbes)),
        length(setdiff(probeIDs, union(categoryProbes, sigProbes)))),
        nrow = 2,
        dimnames = list(
            Query = c("Q_in","Q_out"),
            Database = c("D_in","D_out")))

    test <- fisher.test(mtx)

    result <- data.frame(
        databaseName = dbName,
        OddsRatio = test$estimate,
        P.Value = test$p.value
    )
    return(result)
}

#' calculate fold change given a matrix
calcFoldChange = function(mtx){
	num = mtx[1, 1] / (mtx[1, 1] + mtx[1, 2])
	den = (mtx[1, 1] + mtx[2, 1]) / sum(mtx)
	num / den
}

#' create a volcano plot given a filename, fold change (log2), pvalue (-log10), and title (optional), xlabel (optional), ylabel (optional)
plotVolcano = function(filename, fc, pvalue, title="", xlabel=NA, ylabel=NA) {

    data = na.omit(data.frame(fc=fc,
        pvalue=pvalue))
    rownames(data) = 1:nrow(data)
    colnames(data) = c("log2fc", "pvalue")

    title = gsub('(.{1,80})(\\s|$)', '\\1\n', title)

    if (is.na(xlabel)) {
        xlabel = "log2 fold change"
    }

    if (is.na(ylabel)) {
        ylabel = "-log10 pvalue"
    }

    pdf(filename, height=50, width=50, onefile=FALSE)
    ggplot(data=data, aes(x=fc, y=pvalue, color = cut(pvalue, c(2.995732, Inf) ))) + geom_point() +
    xlab(xlabel) +
    ylab(ylabel) +
    labs(title = title, fill = "pvalue") +
    theme(plot.title = element_text(size=80, face = "bold"),
        axis.text=element_text(size=48),
        axis.title=element_text(size=60),
        legend.title = element_text(size=48),
        legend.text=element_text(size=40)) # +
    # scale_color_manual(name = "pvalue", values = c("[2.995732, Inf)" = "red"), labels = c("> 2.995732"))
    dev.off(0)
}


testEnrichmentFGSEA <- function(category, probeIDs, database, sigProbes, dbName) {
    pathway <- probeIDs[database == category]

    names(sigProbesRank) <- sigProbes
    result <- fgsea(pathways = list(category = pathway), stats = sigProbes)
    result$pathway <- db
    return(result)
}

testEnrichmentSpearman <- function(sigProbes, database, dbName) {
    sharedProbes <- intersect(names(sigProbes), names(database))
    test <- cor.test(
        sigProbes,
        database,
        method = 'spearman'
    )
    result <- data.frame(
        DatabaseAccession = dbName,
        rho = test$estimate,
        P.Value = test$p.value
    )
    return(result)
}
