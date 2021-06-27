#' Test for the enrichment of set of probes (query set) in a given feature (database set)
#' testEnrichmentAll1 tests a default set of categories (database sets) for enrichment
#' in the provided vector of CpGs (querySet).
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param sigProbesRank Numerical ranking for querySet such as P-values. Optional. (Default: NULL)
#' @param databaseSet Vector of database sets of interest. Optional. (Default: NULL)
#'
#' @return A list of two tables, one each for database sets of categorical vs continuous
#' values. The tables rank database sets by odds ratio (Fisher's exact test,
#' categorical querySet and categorical database set), or P-value (FGSEA or Spearman, categorical vs
#'  continuous or continuous vs continuous respectively).
#'
#' @import dplyr
#' @import fgsea
#' @import sesameData
#'
#' @export
testEnrichment1 = function(querySet, databaseSet, universeSet, verbose=FALSE) {

    if (is.numeric(querySet)) { # a named vector of continuous value
        if(is.numeric(databaseSet)) { # numeric db
            if (verbose) {
                cat("Query set: Continuous\tDatabase set: Continuous\t[Spearman test]\n")
            }
            if (length(intersect(names(querySet), names(databaseSet))) == 0) {
                return(NULL)
            }
            results = testEnrichmentSpearman(
                querySet=querySet,
                databaseSet=databaseSet)
        } else {
            if (verbose) {
                cat("Query set: Continuous\tDatabase set: Discrete\t[FGSEA test]\n")
            }
            if (length(intersect(names(querySet), databaseSet)) == 0) {
                return(NULL)
            }
            results = testEnrichmentFGSEA(
                querySet=querySet,
                databaseSet=databaseSet)
        }
    } else { # categorical query
        if(is.numeric(databaseSet)) { # numeric db
            ## do fgsea(switched arguments)
            if (verbose) {
                cat("Query set: Discrete\tDatabase set: Continuous\t[FGSEA test]\n")
            }
            if (length(intersect(querySet, names(databaseSet))) == 0) {
                return(NULL)
            }
            results = testEnrichmentFGSEA(
                querySet=databaseSet,
                databaseSet=querySet)
        } else { # categorical db
            if (verbose) {
                cat("Query set: Discrete\tDatabase set: Discrete\t\t[Fisher exact test]\n")
            }
            if (length(intersect(querySet, databaseSet)) == 0) {
                return(NULL)
            }
            results = testEnrichmentFisher(
                querySet=querySet, 
                databaseSet=databaseSet, 
                universeSet=universeSet)
        }
    }
    return(results)
}


testEnrichmentAll = function(querySet, databaseSets=NULL, platform = c("MM285", "EPIC", "HM450", "HM27"), verbose=FALSE) {

    if (platform == "MM285") {
        # TODO: upload MM285 manifest to sever... use tmp file for now
        # universeSet = readRDS(url("http://zhouserver.research.chop.edu/InfiniumAnnotation/current/MM285/MM285.mm10.manifest.rds"))
        universeSet = readRDS(url("http://zhouserver.research.chop.edu/moyerej/InfiniumAnnotation/MM285.mm10.manifest.rds"))$probeID

    } else if (platform == "EPIC") {
        universeSet = readRDS(url("http://zhouserver.research.chop.edu/InfiniumAnnotation/current/EPIC/EPIC.hg19.manifest.rds"))$probeID
    } else if (platform == "HM450") {
        universeSet = readRDS(url("http://zhouserver.research.chop.edu/InfiniumAnnotation/current/HM450/HM450.hg19.manifest.rds"))$probeID
    } else if (platform == "HM27") {
        universeSet = readRDS(url("http://zhouserver.research.chop.edu/InfiniumAnnotation/current/HM27/HM27.hg19.manifest.rds"))$probeID
    }

    if (is.null(databaseSets)) {
        print("Database set was not defined. Loading in default database sets of transcription factor binding sites.")
        databaseSets = readRDS(url("http://zhouserver.research.chop.edu/kyCG/20210601_MM285_TFBS_ENCODE.rds"))
    }
   
    results = do.call(rbind, lapply(databaseSets, function(databaseSet) testEnrichment1(querySet=querySet, databaseSet=databaseSet, universeSet=universeSet, verbose=verbose)))
    results = results[order(results$p.value, decreasing=F), ]
    return(results)
    ## apply multi-test correction, BH, FDR
}


# testEnrichment = function(setQ, setD, setU) {
testEnrichmentFisher = function(querySet, databaseSet, universeSet) {
    # setD = categoryProbes
    # setQ = querySet
    # setU = universeSet

    mtx = matrix(c(
        length(intersect(querySet, databaseSet)),
        length(setdiff(databaseSet, querySet)),
        length(setdiff(querySet, databaseSet)),
        length(setdiff(universeSet, union(databaseSet, querySet)))),
        nrow = 2,
        dimnames = list(
            querySet = c("Q_in","Q_out"),
            databaseSet = c("D_in","D_out")))

    test <- fisher.test(mtx)

    result <- data.frame(
        estimate = calcFoldChange(mtx),
        p.value = test$p.value,
        test = "fisher"
    )
    return(result)
}


#' calculate fold change given a matrix
calcFoldChange = function(mtx){
    num = mtx[1, 1] / (mtx[1, 1] + mtx[1, 2])
    den = (mtx[1, 1] + mtx[2, 1]) / sum(mtx)
    num / den
}


testEnrichmentFGSEA <- function(querySet, databaseSet) {
    test <- fgsea(pathways=list(pathway=databaseSet), stats=querySet)
    result <- data.frame(
        estimate = test$ES[[1]],
        p.value = test$pval,
        test = "fgsea"
    )
    return(result)
}


testEnrichmentSpearman <- function(querySet, databaseSet) {
    test <- cor.test(
        querySet,
        querySet,
        method = 'spearman'
    )
    result <- data.frame(
        estimate = test$estimate[[1]],
        p.value = test$p.value,
        test = "spearman"
    )
    return(result)
}





