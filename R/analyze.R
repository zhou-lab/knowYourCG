
#' test whether query is significantly enriched in database
testEnrichment = function(setQ, setD, setU) {
    mtx = matrix(c(
        length(intersect(setQ, setD)),
        length(setdiff(setD, setQ)), 
        length(setdiff(setQ, setD)), 
        length(setdiff(setU, union(setD, setQ)))),
        nrow = 2, 
        dimnames = list(
            Query = c("Q_in","Q_out"),
            Database = c("D_in","D_out")))
    list(mtx = mtx, test = fisher.test(mtx))
}

testCategoricalFisher <- function(database, probeIDs, sigProbes, dbName) {
    # function to get unique database categories, filter by overlap with query probes, and call testEnrichment
    out <- list()
    # get unique categories for database. Filter out ones with no overlap with significant probes
    categories <- unique(database[!is.na(database)])
    print(paste0('number of categories: ', length(categories)))
    # count number of overlaps with sig probes per category
    categories <- sapply(
        categories, 
        function(category) {
            length(intersect(sigProbes, names(database)[database == category]))
        }
    )
    # keep only categories with non-zero overlap
    categories <- names(categories)[categories > 0]
    if (length(categories) == 0) stop('no categories with overlapping probes')
    
    # perform test for each category
    out <- lapply(
        categories,
        function(category) {
            setD <- probeIDs[database == category]
            fishertest <- testEnrichment(
                setQ = sigProbes,
                setD = setD,
                setU = probeIDs
            )
            result <- list(
                DatabaseAccession = dbName,
                Category = category,
                TestMatrix = fishertest$mtx,
                Fisher.Test = fishertest$test,
                OddsRatio = fishertest$test$estimate,
                P.Value = fishertest$test$p.value
            )
            return(result)
        }
    )
    
    names(out) <- categories
    return(out)
}

#' test all databaseSet and return a list ranked by enrichment (odds-ratio)
testEnrichmentAll = function(probeIDs, sigProbes, databaseSets = NULL, sig.threshold = .05) {
    # probeIDs: list of all probes
    # sigProbes: list of significant probes
    # assume same index for results and databases
    # master list of databaseSets. Probably shouldn't load this in the function
    suppressMessages(library(readxl))
    suppressMessages(library(dplyr))
    suppressMessages(library(fgsea))
    databaseMaster <- read_xlsx('~/Dropbox/Ongoing_knowYourCpG/20210329_MM285_databaseSets.xlsx')
    # TODO: get list of core database sets if not none given. For now, use temporary list
    if (is.null(databaseSets)) {
        # coredatabaseSets <- getCoreDatabaseSets()
        coreDatabaseSets <- c('20210409_Infinium_Type', '20210416_cpg_density')
        databaseInfo <- databaseMaster[databaseMaster$FileAccession %in% coreDatabaseSets, ]
    }
    
    # apply enrichment tests and get test results
    results <- list(categorical = list(), continuous = list())
    # categorical first
    for (db in databaseInfo$FileAccession[databaseInfo$Format == 'Categorical']) {
        print(db)
        database <- tbk_data(
            idx_fname = '~/Dropbox/Ongoing_knowYourCpG/TBK_INDICES/MM285.idx.gz',
            tbk_fnames = file.path('~/Dropbox/Ongoing_knowYourCpG/DATABASE_SETS/MM285/', paste0(db, '.tbk'))
        )
        results$categorical[[db]] <- testCategoricalFisher(database = database, probeIDs = probeIDs, 
                                                           sigProbes = sigProbes, dbName = db)
        rm(database)
    }
    
    # sort results by output
    # results <- sort(results, decreasing = TRUE)
    
    # continuous value databases
    results$continuous <- lapply(
        databaseInfo$FileAccession[databaseInfo$Format == 'Continuous'],
        function(db) {
            print(db)
            database <- tbk_data(
                idx_fname = '~/Dropbox/Ongoing_knowYourCpG/TBK_INDICES/MM285.idx.gz',
                tbk_fnames = file.path('~/Dropbox/Ongoing_knowYourCpG/DATABASE_SETS/MM285/', paste0(db, '.tbk'))
            )
            result <- fgsea(pathways = list(sigProbes = sigProbes), stats = database)
            result$pathway <- db
            return(result)
        }
    ) %>%
        bind_rows()
    
    return(results)
}


