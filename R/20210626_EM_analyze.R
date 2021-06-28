#' testEnrichment1 tests for the enrichment of set of probes (query set) in a 
#' single given feature (database set)
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet Vector of probes corresponding to a single database set of interest.
#' @param universeSet Vector of probes in the universe set containing all of 
#' the probes to be considered in the test.
#' @param verbose Logical value indicating whether to display intermediate 
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return One list containing features corresponding the test estimate, 
#' p-value, and type of test.
#'
#' @import fgsea
#'
#' @export
testEnrichment1 = function(querySet, databaseSet, universeSet, estimate.type="ES", verbose=FALSE) {
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
                databaseSet=databaseSet, estimate.type=estimate.type)
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
                databaseSet=querySet, estimate.type=estimate.type)
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


#' testEnrichmentAll tests for the enrichment of set of probes (query set) in a number of features (database sets).
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSets List of vectors of probes corresponding to multiple 
#' database sets. Optional. (Default: NULL)
#' @param platform String corresponding to the type of platform to use. Either 
#' MM285, EPIC, HM450, or HM27.
#' @param verbose Logical value indicating whether to display intermediate 
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return One list containing features corresponding the test estimate, p-value, and type of test.
#'
#' @export
testEnrichmentAll = function(querySet, databaseSets=NULL, platform = c("MM285", "EPIC", "HM450", "HM27"), estimate.type="ES", verbose=FALSE) {
    # options(timeout=1000)
    if (platform == "MM285") {
        universeSet = get(load(url("http://zhouserver.research.chop.edu/sesameData/MM285.mm10.manifest.rda")))$probeID
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
   
    results = do.call(rbind, lapply(databaseSets, function(databaseSet) testEnrichment1(querySet=querySet, databaseSet=databaseSet, universeSet=universeSet, estimate.type=estimate.type, verbose=verbose)))
    results = results[order(results$p.value, decreasing=F), ]
    return(results)
    ## apply multi-test correction, BH, FDR
}


#' testEnrichmentFisher uses Fisher's exact test to estimate the association 
#' between two categorical variables.
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet Vector of probes corresponding to a single database set 
#' of interest.
#' @param universeSet Vector of probes in the universe set containing all of 
#' the probes to be considered in the test. (Default: NULL)
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test 
#' for the given results.
testEnrichmentFisher = function(querySet, databaseSet, universeSet) {
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


#' calcFoldChange calculates fold change given a 2x2 matrix of counts.
#'
#' @param mtx 2x2 matrix of values corresponding to overlapping counts between 
#' two sets of a categorical variable.
#'
#' @return A numerical value corresponding to the fold change enrichment,
calcFoldChange = function(mtx){
    num = mtx[1, 1] / (mtx[1, 1] + mtx[1, 2])
    den = (mtx[1, 1] + mtx[2, 1]) / sum(mtx)
    num / den
}

#' testEnrichmentFGSEA uses the FGSEA test to estimate the association of a 
#' categorical variable against a continuous variable.
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet Vector of probes corresponding to a single database set 
#' of interest.
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test 
#' for the given results.
testEnrichmentFGSEA <- function(querySet, databaseSet, p.pvalue.adj=FALSE, estimate.type="ES") {
    test <- fgsea(pathways=list(pathway=databaseSet), stats=querySet)

    if (p.pvalue.adj) {
        p.value = test$padj
    } else {
        p.value = test$pval
    }    

    if (estimate.type == "log2err") {
        estimate = test$log2err
    } else if (estimate.type == "NES") {
        estimate = test$NES
    } else if (estimate.type == "leadingEdge") {
        estimate = test$leadingEdge
    } else if (estimate.type == "ES") {
        estimate = test$ES
    } else {
        print(sprintf("Incorrect estimate.type: [", estimate.type, "].", sep=""))
        return(NULL)
    }

    result <- data.frame(
        estimate = estimate,
        p.value = p.value,
        test = "fgsea"
    )
    return(result)
}

#' testEnrichmentSpearman uses the spearman test to estimate the association 
#' between two continuous variables.
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet Vector of probes corresponding to a single database set 
#' of interest.
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test 
#' for the given results.
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





