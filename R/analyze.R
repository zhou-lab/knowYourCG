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
testEnrichment1 = function(sigProbes, database) { # sigProbesRank = NULL, databaseSets = NULL) {

    # probe IDs: mouseMethylation285_probeIDs

    # default database sets: test_defaultDatabaseSets
    # TODO: finalize database sets

    if (is.numeric(sigProbes)) { # a named vector of continuous value
        if(is.numeric(database)) { # numeric db
            ## do fgsea(switched arguments)
            ## correlation-based
            ## spearman
        } else {
            print('FGSEA for continuous database sets')

            ## continuous database set
            bind_rows(lapply(
                names(test_defaultDatabaseSets$continuous),
                function(db) {
                    cbind(fgsea(
                        pathways = list(sigProbes = sigProbes),
                        stats = test_defaultDatabaseSets$continuous[[db]]
                    ), pathway = db)
                }
            ))
        }
    } else { # categorical query
        if(is.numeric(database)) { # numeric db
            ## do fgsea(switched arguments)

            print("Fisher's exact test for categorical database sets")
            database = sesameDataGet("data/www/html/kyCG") # instead of test_defaultDatabaseSets
            ## if only significant probes provided, Fisher (categorical db set) and fgsea (continuous db set)
                                        # categorical database sets first
            bind_rows(lapply(
                names(test_defaultDatabaseSets$categorical),
                function(db) {
                    testEnrichmentFisher(
                        probeIDs = mouseMethylation285_probeIDs,
                        categoryProbes = test_defaultDatabaseSets$categorical[[db]],
                        sigProbes = sigProbes,
                        dbName = db
                    )
                    ## add size of the db
                }
            ))

                                        # sort results by output
            results$categorical <- results$categorical[order(results$categorical$OddsRatio, decreasing = TRUE), ]
            results$categorical <- left_join(results$categorical, databaseSetManifest, by = c('DatabaseAccession' = 'FileAccession'))
        } else { # categorical db
        }
    }
    ## else { # if significant probes AND ranking provided, fgsea (categorical db set) and spearman (continuous db set)
    ##     results$categorical <- lapply(
    ##         categoricalDatabases,
    ##         function(db) {
    ##             print(db)
    ##             database = sesameDataGet("kycg/20210601_TFBS_ENCODE")
    ##             ## database <- tbk_data(
    ##             ##     idx_fname = '~/Dropbox/Ongoing_knowYourCpG/TBK_INDICES/MM285.idx.gz',
    ##             ##     tbk_fnames = file.path('~/Dropbox/Ongoing_knowYourCpG/DATABASE_SETS/MM285/', paste0(db, '.tbk'))
    ##             ## )

    ##             categories <- getCategories(database)

    ##             fgseaResult <- lapply(
    ##                 categories,
    ##                 testEnrichmentFGSEA,
    ##                 probeIDs = probeIDs,
    ##                 database = database,
    ##                 sigProbes = sigProbes,
    ##                 dbName = db,
    ##                 sigProbesRank = sigProbesRank
    ##             ) %>%
    ##                 bind_rows()
    ##             return(fgseaResult)
    ##         }
    ##     ) %>%
    ##         bind_rows()
    ##     # sort results by output
    ##     results$categorical <- results$categorical[order(results$categorical$p.value), ]

    ##     results$continuous <- lapply(
    ##         continuousDatabases,
    ##         function(db) {
    ##             print(db)
    ##             database <- tbk_data(
    ##                 idx_fname = '~/Dropbox/Ongoing_knowYourCpG/TBK_INDICES/MM285.idx.gz',
    ##                 tbk_fnames = file.path('~/Dropbox/Ongoing_knowYourCpG/DATABASE_SETS/MM285/', paste0(db, '.tbk'))
    ##             )
    ##             result <- testEnrichmentSpearman(
    ##                 sigProbes = sigProbes,
    ##                 sigProbesRank = sigProbesRank,
    ##                 database = database,
    ##                 dbName = db
    ##             )
    ##             return(result)
    ##         }
    ##     ) %>%
    ##         bind_rows()
    ## }

    return(results)
}



testEnrichmentAll = function(query) {
    if (!is.null(databaseSets)) {
        # TODO: support custom lists of databaseSets
    }
    databases = sesameDataGet("kycg/20210601_TFBS_ENCODE") # instead of test_defaultDatabaseSets
    lapply(databases, function(db) testEnrichment1(query=query, db=db))
    ## apply multi-test correction, BH, FDR
}


# testEnrichment = function(setQ, setD, setU) {
testEnrichmentFisher = function(probeIDs, categoryProbes, sigProbes, dbName) {
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
        DatabaseAccession = dbName,
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


testEnrichmentFGSEA <- function(category, probeIDs, database, sigProbes, dbName, sigProbesRank) {
    pathway <- probeIDs[database == category]

    names(sigProbesRank) <- sigProbes
    result <- fgsea(pathways = list(category = pathway), stats = sigProbesRank)
    result$pathway <- db
    return(result)
}

testEnrichmentSpearman <- function(sigProbes, sigProbesRank, database, dbName) {
    sharedProbes <- intersect(sigProbes, names(database))
    names(sigProbesRank) <- sigProbes
    test <- cor.test(
        sigProbesRank[names(sigProbesRank) %in% sharedProbes],
        database[names(database) %in% sharedProbes],
        method = 'spearman'
    )
    result <- data.frame(
        DatabaseAccession = dbName,
        rho = test$estimate,
        P.Value = test$p.value
    )
    return(result)
}
