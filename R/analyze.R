
#' test whether query is significantly enriched in database
# testEnrichment = function(setQ, setD, setU) {
testEnrichmentFisher = function(category, probeIDs, database, sigProbes, dbName) {
    setD <- probeIDs[database == category]
    setQ = sigProbes
    setU = probeIDs
    
    mtx = matrix(c(
        length(intersect(setQ, setD)),
        length(setdiff(setD, setQ)), 
        length(setdiff(setQ, setD)), 
        length(setdiff(setU, union(setD, setQ)))),
        nrow = 2, 
        dimnames = list(
            Query = c("Q_in","Q_out"),
            Database = c("D_in","D_out")))
    
    test <- fisher.test(mtx)
    
    result <- data.frame(
        DatabaseAccession = dbName,
        Category = category,
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

getCategories <- function(database) {
    suppressMessages(library(dplyr))
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
    return(categories)


#' test all databaseSet and return a list ranked by enrichment (odds-ratio)
testEnrichmentAll = function(probeIDs, sigProbes, sigProbesRank = NULL, databaseSets = NULL) {
    # probeIDs: list of all probes
    # sigProbes: list of significant probes
    # sigProbesRank: corresponding ranking or other metric of interest for sigProbes, e.g. p values
    # databaseSets: vector of databaseSets of interest. Overrides default core database sets

    # assume same index for results and databases
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
    
    results <- list(categorical = data.frame(), continuous = data.frame())
    
    categoricalDatabases <- databaseInfo$FileAccession[databaseInfo$Format == 'Categorical']
    continuousDatabases <- databaseInfo$FileAccession[databaseInfo$Format == 'Continuous']
    
    # if only significant probes provided, Fisher (categorical db set) and fgsea (continuous db set)
    if (is.null(sigProbesRank)) {
        results$categorical <- lapply(
            categoricalDatabases,
            function(db) {
                print(db)
                database <- tbk_data(
                    idx_fname = '~/Dropbox/Ongoing_knowYourCpG/TBK_INDICES/MM285.idx.gz',
                    tbk_fnames = file.path('~/Dropbox/Ongoing_knowYourCpG/DATABASE_SETS/MM285/', paste0(db, '.tbk'))
                )
                
                categories <- getCategories(database)
                
                fisherResult <- lapply(
                    categories,
                    testEnrichmentFisher,
                    probeIDs = probeIDs,
                    database = database,
                    sigProbes = sigProbes,
                    dbName = db
                ) %>%
                    bind_rows()
                return(fisherResult)
            }
        ) %>%
            bind_rows()
        # sort results by output
        results$categorical <- results$categorical[order(results$categorical$OddsRatio), ]
        
        results$continuous <- lapply(
            continuousDatabases,
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


