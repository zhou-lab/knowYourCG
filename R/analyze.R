#' Test for enriched categories of CpGs
#' testEnrichmentAll tests a default set of categories (database sets) for enrichment
#' in the provided vector of CpGs (sigProbes).
#'
#' @param sigProbes Vector of probes of interested (e.g., significant probes)
#' @param sigProbesRank Numerical ranking for sigProbes such as P-values. Optional. (Default: NULL)
#' @param databaseSets Vector of database sets of interest. Optional. (Default: NULL)
#'
#' @return A list of two tables, one each for database sets of categorical vs continuous
#' values. The tables rank database sets by odds ratio (Fisher's exact test,
#' categorical sigProbes and categorical database set), or P-value (FGSEA or Spearman, categorical vs
#'  continuous or continuous vs continuous respectively).
#'
#' @import dplyr
#' @import fgsea
#'
#' @export
testEnrichmentAll = function(sigProbes, sigProbesRank = NULL, databaseSets = NULL) {

    # probe IDs: mouseMethylation285_probeIDs

    if (!is.null(databaseSets)) {
        # TODO: support custom lists of databaseSets
    }

    # default database sets: test_defaultDatabaseSets
    # TODO: finalize database sets

    results <- list(categorical = data.frame(), continuous = data.frame())

    # if only significant probes provided, Fisher (categorical db set) and fgsea (continuous db set)
    if (is.null(sigProbesRank)) {
        # categorical database sets first
        print("Fisher's exact test for categorical database sets")
        results$categorical <- lapply(
            names(test_defaultDatabaseSets$categorical),
            function(db) {
                fisherResult <- testEnrichmentFisher(
                    probeIDs = mouseMethylation285_probeIDs,
                    categoryProbes = test_defaultDatabaseSets$categorical[[db]],
                    sigProbes = sigProbes,
                    dbName = db
                )

                return(fisherResult)
            }
        ) %>%
            bind_rows()
        # sort results by output
        results$categorical <- results$categorical[order(results$categorical$OddsRatio), ]

        # continuous database set
        print('FGSEA for continuous database sets')
        results$continuous <- lapply(
            names(test_defaultDatabaseSets$continuous),
            function(db) {

                result <- fgsea(
                    pathways = list(sigProbes = sigProbes),
                    stats = test_defaultDatabaseSets$continuous[[db]]
                )
                result$pathway <- db

                return(result)
            }
        ) %>%
            bind_rows()
    } else { # if significant probes AND ranking provided, fgsea (categorical db set) and spearman (continuous db set)
        results$categorical <- lapply(
            categoricalDatabases,
            function(db) {
                print(db)
                database <- tbk_data(
                    idx_fname = '~/Dropbox/Ongoing_knowYourCpG/TBK_INDICES/MM285.idx.gz',
                    tbk_fnames = file.path('~/Dropbox/Ongoing_knowYourCpG/DATABASE_SETS/MM285/', paste0(db, '.tbk'))
                )

                categories <- getCategories(database)

                fgseaResult <- lapply(
                    categories,
                    testEnrichmentFGSEA,
                    probeIDs = probeIDs,
                    database = database,
                    sigProbes = sigProbes,
                    dbName = db,
                    sigProbesRank = sigProbesRank
                ) %>%
                    bind_rows()
                return(fgseaResult)
            }
        ) %>%
            bind_rows()
        # sort results by output
        results$categorical <- results$categorical[order(results$categorical$p.value), ]

        results$continuous <- lapply(
            continuousDatabases,
            function(db) {
                print(db)
                database <- tbk_data(
                    idx_fname = '~/Dropbox/Ongoing_knowYourCpG/TBK_INDICES/MM285.idx.gz',
                    tbk_fnames = file.path('~/Dropbox/Ongoing_knowYourCpG/DATABASE_SETS/MM285/', paste0(db, '.tbk'))
                )
                result <- testEnrichmentSpearman(
                    sigProbes = sigProbes,
                    sigProbesRank = sigProbesRank,
                    database = database,
                    dbName = db
                )
                return(result)
            }
        ) %>%
            bind_rows()
    }

    return(results)
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
