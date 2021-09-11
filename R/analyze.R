baseurl = "http://zhouserver.research.chop.edu"
# source("R/data.R")

#' getDatabaseSetPairwiseDistance tests for the pariwise overlap between given
#' list of database sets using a distance metric.
#'
#' @param databaseSets List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element. Optional.
#' (Default: NA)
#' @param metric String representing the similarity metric to use. Optional.
#' (Default: "Jaccard").
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return An upper triangular matrix containing a metric (Jaccard) comparing
#' the pairwise distances between database sets.
#'
#' @examples
#' databaseSets = list(a=c("a", "b"), b=c("a", "e", "f"), c=c("q", "a"))
#' getDatabaseSetPairwiseDistance(databaseSets)
#'
#' @export
getDatabaseSetPairwiseDistance = function(databaseSets=NA,
                                          metric="Jaccard",
                                          verbose=FALSE) {
    ndatabaseSets = length(databaseSets)
    names = names(databaseSets)
    m = matrix(0, nrow=ndatabaseSets, ncol=ndatabaseSets)
    colnames(m) = names
    rownames(m) = names
    for (i in seq(ndatabaseSets - 1)) { 
        for (j in seq(i + 1, ndatabaseSets)) {
            cat(i, j, '\n')
            m[i, j] = length(intersect(databaseSets[[i]], databaseSets[[j]])) /
                length(union(databaseSets[[i]], databaseSets[[j]]))
        }
    }
    return(m)
}

#' getDatabaseSetOverlap tests for the overlap of set of probes (querySet) in a
#' single given feature (database set)
#'
#' @param querySet Vector of probes corresponding to a single database set
#' of interest.
#' @param databaseSets List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element. Optional.
#' (Default: NA)
#' @param release Integer indicating the release number of the database set. 
#' Optional. (Default: NA).
#' @param array String corresponding to the type of array to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @param dev Logical value indicating whether to use development version.
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return A sparse data.frame containing all of the meta data from all database
#' sets.
#'
#' @examples
#' querySet=c("cg29176188_TC21", "cg29176794_TC21")
#' getDatabaseSetOverlap(querySet, release=2)
#'
#' @export
getDatabaseSetOverlap = function(querySet,
                           databaseSets=NA,
                           release=NA,
                           array=NA,
                           dev=TRUE, 
                           verbose=TRUE) {
    if (all(is.na(databaseSets))) {
        databaseSets = getDatabaseSets(release=release, dev=dev)
    }
    if (all(is.na(databaseSets))) {
        return(NULL)
    }

    metadata = as.data.frame(
        do.call(rbind,
                lapply(databaseSets,
                       function(databaseSet) {
                           rowmeta = attr(databaseSet, "meta")
                           if (!is.null(rowmeta)) 
                               rowmeta = c(meta=TRUE, rowmeta)
                            else
                                rowmeta = c(meta=FALSE, rowmeta)
                           nQ = length(querySet)
                           nD = length(databaseSet)
                           overlap = length(intersect(querySet, databaseSet))
                            rowmeta = c(rowmeta, nQ=nQ, nD=nD, overlap=overlap)
                           return(rowmeta)
                       })
        ))
    
    metadata$meta = as.logical(metadata$meta)

    return(metadata)
}

#' testEnrichment1 tests for the enrichment of set of probes (query set) in a
#' single given feature (database set)
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet Vector corresponding to the database sets of
#' interest with associated meta data as an attribute to each element.
#' @param universeSet Vector of probes in the universe set containing all of
#' the probes to be considered in the test.
#' @param estimate.type String indicating the estimate to report. (Default:
#' "ES")
#' @param p.value.adj Logical value indicating whether to report the adjusted
#' p-value. (Default: FALSE)
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @import utils
#'
#' @return One list containing features corresponding the test estimate,
#' p-value, and type of test.
testEnrichment1 = function(querySet, databaseSet, universeSet,
                           estimate.type="ES", p.value.adj=FALSE,
                           verbose=FALSE) {
    if (is.numeric(querySet)) { # a named vector of continuous value
        if(is.numeric(databaseSet)) { # numeric db
            if (verbose) {
                cat("Query set: Continuous\t",
                    "Database set: Continuous\t[Spearman test]\n")
            }
            results = testEnrichmentSpearman(
                querySet=querySet,
                databaseSet=databaseSet)
        } else {
            if (verbose) {
                cat("Query set: Continuous\t",
                    "Database set: Discrete\t\t[FGSEA test]\n")
            }
            results = testEnrichmentFGSEA(
                querySet=querySet,
                databaseSet=databaseSet,
                p.value.adj=p.value.adj,
                estimate.type=estimate.type)
        }
    } else { # categorical query
        if(is.numeric(databaseSet)) { # numeric db
            ## do fgsea(switched arguments)
            if (verbose) {
                cat("Query set: Discrete\t", 
                    "Database set: Continuous\t[FGSEA test]\n")
            }
            results = testEnrichmentFGSEA(
                querySet=databaseSet,
                databaseSet=querySet,
                p.value.adj=p.value.adj,
                estimate.type=estimate.type)
        } else { # categorical db
            if (verbose) {
                cat("Query set: Discrete\t",
                    "Database set: Discrete\t\t[Fisher exact test]\n")
            }
            results = testEnrichmentFisher(
                querySet=querySet,
                databaseSet=databaseSet,
                universeSet=universeSet)
        }
    }
    return(results)
}


#' testEnrichmentAll tests for the enrichment of set of probes (query set) in
#' a number of features (database sets).
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSets List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element. Optional.
#' (Default: NA)
#' @param universeSet Vector of probes in the universe set containing all of
#' the probes to be considered in the test. If it is not provided, it will be
#' inferred from the provided array. (Default: NA).
#' @param array String corresponding to the type of array to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @param estimate.type String indicating the estimate to report. (Default:
#' "ES")
#' @param p.value.adj Logical value indicating whether to report the adjusted
#' p-value. (Default: FALSE).
#' @param n.fdr Integer corresponding to the number of comparisons made. 
#' Optional. (Default: NA).
#' @param release Integer indicating the release number of the database set. 
#' Optional. (Default: NA).
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE).
#'
#' @return One list containing features corresponding the test estimate,
#' p-value, and type of test.
#'
#' @examples
#' testEnrichmentAll(c("cg0000029"))
#'
#' @export
testEnrichmentAll = function(querySet, databaseSets=NA, universeSet=NA,
                             array=NA, estimate.type="ES", p.value.adj=FALSE,
                             n.fdr=NA, release=2, verbose=FALSE) {
    options(timeout=1000)
    if (all(is.na(universeSet))) {
        if (verbose) {
            cat("The universeSet was not defined.", 
                  "Loading in universeSet based on array.")
        }
        if (is.na(array)) {
            if (verbose) {
                cat("The array was not defined.",
                      "Inferring array from probeIDs.")
            }
            array = inferarrayFromProbeIDs(querySet)
        }
        universeSet = getUniverseSet(array)
    }

    # Pull from release instead

    if (all(is.na(databaseSets))) {
        if (!is.na(release)) {
            if (verbose)
                print("Database set was not defined. Loading in database
                  sets from the given release.")
            databaseSets = getDatabaseSets(release=release)
        } else {
            if (verbose)
                print("Database set was not defined. Loading in database
                  sets from the most recent release.")
            databaseSets = getDatabaseSets(release=2)
        }
    }

    results = data.frame(
        do.call(rbind,
                lapply(databaseSets,
                       function(databaseSet) testEnrichment1(
                           querySet=querySet,
                           databaseSet=databaseSet,
                           universeSet=universeSet,
                           p.value.adj=p.value.adj,
                           estimate.type=estimate.type,
                           verbose=verbose)
                       )
                ))
    
    if (is.na(n.fdr))
        results$p.adjust.fdr = p.adjust(results$p.value, method='fdr')
    else
        results$p.adjust.fdr = p.adjust(results$p.value, method='fdr', n=n.fdr)
    
    metadata = data.frame(
        do.call(rbind,
                lapply(databaseSets,
                       function(databaseSet) {
                           output = attr(databaseSet, "meta")
                           if (!is.null(output)) 
                               return(append(c(meta=TRUE), output))
                           return(c(meta=FALSE))
                       })
                )
        )
    
    rank = list()
    rank$estimate.rank = rank(-results$estimate, ties.method='first')
    rank$p.value.rank = rank(results$p.value, ties.method='first')
    rank$overlap.rank = rank(results$overlap, ties.method='first')
    rank$max.rank = apply(data.frame(rank), 1, max)
    rank$mean.rank = apply(data.frame(rank), 1, mean)
    rank = data.frame(rank, row.names=row.names(results))
    output = cbind(results, rank, metadata)

    output = output[order(output$p.value, decreasing=FALSE), ]
    return(output)
    ## apply multi-test correction, BH, FDR
}


#' testEnrichmentGene tests for the enrichment of set of probes
#' (querySet) in gene regions.
#'
#' @param querySet Vector of probes of interest (e.g., probes belonging to a
#' given array)
#' @param array String corresponding to the type of array to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set querySet (Default: NA)
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return One list containing features corresponding the test estimate,
#' p-value, and type of test.
#'
#' @examples
#' testEnrichmentGene(c("cg0000029"), array="EPIC")
#'
#' @export
testEnrichmentGene = function(querySet, array=NA, verbose=FALSE) {
    if (is.na(array)) {
        if (verbose) {
            print("The array was not defined. Inferring array from probeIDs.")
        }
        array = inferarrayFromProbeIDs(querySet)
    }
    probeID2gene = getProbeID2Gene(array)

    databaseSetNames = probeID2gene$genesUniq[match(querySet, 
                                                    probeID2gene$probeID)]

    databaseSetNames = na.omit(unique(
        unlist(lapply(databaseSetNames,
                      function(databaseSetName) {
                          strsplit(databaseSetName, ";")
                          }))))

    if (length(databaseSetNames) == 0) return(NULL)

    databaseSets = getDatabaseSets(group="Gene", array=array)

    n = length(databaseSets)
    
    databaseSets = databaseSets[names(databaseSets) %in% databaseSetNames]
    
    return(testEnrichmentAll(querySet, databaseSets, array=array, n.fdr=n))
}

#' inferarrayFromProbeIDs infers the Infinium MicroArray array using the
#' given probeIDs
#'
#' @param probeIDs Vector of probes of interest (e.g., probes belonging to a
#' given array)
#'
#' @return String corresponding to the inferred array
inferarrayFromProbeIDs = function(probeIDs) {
    sig = get(load(url(sprintf("%s/sesameData/probeIDSignature.rda", baseurl))))
    names(which.max(vapply(
        sig, function(x) sum(probeIDs %in% x), integer(1))))
}

#' testEnrichmentFisher uses Fisher's exact test to estimate the association
#' between two categorical variables.
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet Vectors corresponding to the database set of
#' interest with associated meta data as an attribute to each element.
#' @param universeSet Vector of probes in the universe set containing all of
#' the probes to be considered in the test. (Default: NULL)
#'
#' @import stats
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
testEnrichmentFisher = function(querySet, databaseSet, universeSet) {
    test = "fisher"
    if (length(intersect(querySet, databaseSet)) == 0) {
        return(data.frame(estimate=0,
                    p.value=1,
                    test=test,
                    nQ=length(querySet),
                    nD = length(databaseSet),
                    overlap=0
                    ))
    }

    mtx = matrix(c(
        length(intersect(querySet, databaseSet)),
        length(setdiff(databaseSet, querySet)),
        length(setdiff(querySet, databaseSet)),
        length(setdiff(universeSet, union(databaseSet, querySet)))),
        nrow = 2,
        dimnames = list(
            querySet = c("Q_in","Q_out"),
            databaseSet = c("D_in","D_out")))

    res = fisher.test(mtx)

    result = data.frame(
        estimate = calcFoldChange(mtx),
        p.value = res$p.value,
        test = test,
        nQ = length(querySet),
        nD = length(databaseSet),
        overlap = length(intersect(querySet, databaseSet))
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
#' @param p.value.adj Logical value indicating whether to report the adjusted
#' p-value. (Default: FALSE).
#' @param estimate.type String indicating the estimate to report. Optional.
#' (Default: "ES").
#'
#'
#' @import fgsea
#' @import BiocGenerics
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
testEnrichmentFGSEA = function(querySet, databaseSet, p.value.adj=FALSE,
                               estimate.type="ES") {
    test="fgsea"
    overlap = length((intersect(names(querySet), databaseSet)))
    if (overlap == 0) {
        return(data.frame(estimate=0,
                    p.value=1,
                    test=test,
                    nQ=length(databaseSet),
                    nD=length(querySet),
                    overlap=overlap
                    ))
    }
    res = suppressWarnings(fgsea(pathways=list(pathway=databaseSet), 
                                 stats=querySet))

    if (p.value.adj) {
        p.value = res$padj
    } else {
        p.value = res$pval
    }

    if (estimate.type == "log2err") {
        estimate = res$log2err
    } else if (estimate.type == "NES") {
        estimate = res$NES
    } else if (estimate.type == "leadingEdge") {
        estimate = res$leadingEdge
    } else if (estimate.type == "ES") {
        estimate = res$ES
    } else {
        print(sprintf("Incorrect estimate.type: [%s].", estimate.type))
        return(NULL)
    }

    result = data.frame(
        estimate = estimate,
        p.value = p.value,
        test = test,
        nQ=length(databaseSet),
        nD=length(querySet),
        overlap=overlap
    )
    return(result)
}


#' testEnrichmentSpearman uses the Spearman statistical test to estimate the 
#' association between two continuous variables.
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet List of vectors corresponding to the database set of
#' interest with associated meta data as an attribute to each element.
#'
#' @import stats
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
testEnrichmentSpearman = function(querySet, databaseSet) {
    test = "spearman"
    if (length(intersect(names(querySet), names(databaseSet))) == 0) {
        return(data.frame(estimate=0,
                    p.value=1,
                    test=test,
                    nQ=length(querySet),
                    nD=length(databaseSet),
                    overlap=0
                    ))
    }

    databaseSet = databaseSet[match(names(querySet), names(databaseSet))]

    res = cor.test(
        querySet,
        databaseSet,
        method = test
    )
    result = data.frame(
        estimate = res$estimate[[1]],
        p.value = res$p.value,
        test = test
    )
    return(result)
}
