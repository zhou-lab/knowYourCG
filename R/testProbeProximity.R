getDistance <- function(v) {
    if (length(v) == 1) return(NA)
    vapply(seq(v), function(x) v[x + 1] - v[x], numeric(1))
}

getPairwiseDistance <- function(gr,q)
{
    df <- as.data.frame(sort(gr[q,])) |>
        dplyr::group_by(.data$seqnames) |>
        dplyr::mutate(distance=getDistance(.data$start))
}


#' testProbeProximity tests if a query set of probes share closer
#' genomic proximity than if randomly distributed
#'
#' @param probeIDs Vector of probes of interest (e.g., significant probes)
#' @param gr GRanges to draw samples and compute genomic distances
#' @param iterations Number of random samples to generate null distribution
#' (Default: 100).
#' @param bin_size the poisson interval size for computing neighboring hits
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @return list containing a dataframe for the poisson statistics and a
#' data frame for the probes in close proximity
#' @importFrom dplyr mutate group_by 
#' @importFrom rlang .data
#' @examples
#' sesameData::sesameDataCache(data_titles=
#' c("KYCG.MM285.tissueSignature.20211211","MM285.address","probeIDSignature"))
#' library(SummarizedExperiment)
#' df <- rowData(sesameData::sesameDataGet('MM285.tissueSignature'))
#' probes <- df$Probe_ID[df$branch == "B_cell"]
#' res <- testProbeProximity(probeIDs=probes,platform="MM285")
#' sesameData::sesameDataGet_resetEnv()
#'
#' @export
testProbeProximity <- function (probeIDs,gr=NULL,platform=NULL,iterations=100,
                                bin_size=1500)
{
    if (is.null(gr)) {
        if (is.null(platform)) {
            platform <- inferPlatformFromProbeIDs(probeIDs)
        }
        gr <- sesameData_getManifestGRanges(platform)
    }
    gr_q <- getPairwiseDistance(gr,probeIDs)
    qd <- gr_q[["distance"]]
    qd <- qd[!is.na(qd)]
    q_hits <- sum(abs(qd) <= bin_size)

    if (q_hits == 0) {
        stats <- data.frame(
            num_query=length(probeIDs),
            hits_query=q_hits,
            lambda=NA,
            p.val=1
        )
        return(list(Stats=stats,Clusters=NA))
    }

    null_dist <- do.call(c, lapply(seq(iterations), function(x) {
        query2 <- sample(names(gr), length(probeIDs))
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
        nQ=length(probeIDs),
        Hits=q_hits,
        Lambda=lambda,
        P.val=pval
    )

    ind <- which(abs(gr_q[["distance"]]) <= bin_size)
    ind <- unique(c(ind, ind + 1)) |> sort()

    clusters <- gr_q[ind,] |>
        as.data.frame() |>
        dplyr::select(.data$seqnames,.data$start,.data$end,.data$distance)

    list(Stats=stats,Clusters=clusters)

}
