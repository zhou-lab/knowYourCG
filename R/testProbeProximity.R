

getDistance <- function(v) {
    if (length(v) == 1) return(NA)
    vapply(seq(v), function(x) v[x + 1] - v[x], numeric(1))
}

getPairwiseDistance <- function(gr,q) #takes in GRanges
{
    df <- as.data.frame(sort(gr[q,])) %>%
        dplyr::group_by(seqnames) %>%
        dplyr::mutate(distance=getDistance(start))
}


#' testProbeProximity tests if a query set of probes share closer genomic proximity than
#' if randomly distributed
#'
#' @param query Vector of probes of interest (e.g., significant probes)
#' @param gr GRanges to draw samples and compute genomic distances
#' @param iterations Number of random samples to generate null distribution
#' (Default: 100).
#' @param bin_size the poisson interval size for computing neighboring hits
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @return list containing a dataframe for the poisson statistics and a
#' data frame for the probes in close proximity
#' @importFrom dplyr mutate group_by %>%
#' @examples
#'
#' library(SummarizedExperiment)
#' df <- rowData(sesameDataGet('MM285.tissueSignature'))
#' query <- df$Probe_ID[df$branch == "B_cell"]
#' res <- testProbeProximity(query,platform="MM285")
#' sesameDataGet_resetEnv()
#'
#' @export
testProbeProximity <- function (query,gr=NULL,platform=NULL,iterations=100,
                                bin_size=1500)
{
    if (is.null(gr)) {
        if (is.null(platform)) {
            platform = inferPlatformFromProbeIDs(query)
        }
        gr <- sesameData_getManifestGRanges(platform)
    }
    gr_q <- getPairwiseDistance(gr,query)
    qd <- gr_q[["distance"]]
    qd <- qd[!is.na(qd)]
    q_hits <- sum(abs(qd) <= bin_size)

    if (q_hits == 0) {
        stats <- data.frame(
            num_query=length(query),
            hits_query=q_hits,
            lambda=NA,
            p.val=1
        )
        return(list(Stats=stats,Clusters=NA))
    }

    null_dist <- do.call(c, lapply(1:iterations, function(x) {
        query2 <- sample(names(gr), length(query))
        gr_q2 <- getPairwiseDistance(gr,query2)
        q2d <- gr_q2[["distance"]]
        q2d <- q2d[!is.na(q2d)]
        null_hits <- sum(abs(q2d) <= bin_size)
    }))

    lambda <- mean(null_dist)
    if (lambda == 0) {
        lambda <- 0.001
        warning("Forcing lambda to .001")
    }

    pval <- ppois(
        q = q_hits,
        lambda = lambda,
        lower.tail=FALSE
    )

    stats <- data.frame(
        num_query=length(query),
        hits_query=q_hits,
        lambda=lambda,
        p.val=pval
    )

    ind <- which(abs(gr_q[["distance"]]) <= bin_size)
    ind <- unique(c(ind, ind + 1)) %>% sort()

    clusters <- gr_q[ind,] %>%
        as.data.frame() %>%
        dplyr::select(seqnames,start,end,distance)

    list(Stats=stats,Clusters=clusters)

}




